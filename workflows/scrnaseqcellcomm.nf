/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { paramsSummaryMap       } from 'plugin/nf-schema'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_scrnaseqcellcomm_pipeline'

include { PREP_DATA              } from '../subworkflows/local/prep_data'
include { PREP_DATA_ALT          } from '../subworkflows/local/prep_data_alt'
include { CELL2CELL              } from '../subworkflows/local/cell2cell'
include { CELLCHAT               } from '../subworkflows/local/cellchat'
include { CPDB                   } from '../subworkflows/local/cpdb'
include { LIANA                  } from '../subworkflows/local/liana'
include { CONSENSUS              } from '../subworkflows/local/consensus'
include { AGGREGATION            } from '../subworkflows/local/aggregation'
include { UTILS_SAVE_AS_XLSX     } from '../modules/local/utils/save_as_xlsx'
include { RUN_CCI                } from '../subworkflows/local/run_cci/main.nf'
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow SCRNASEQCELLCOMM {

    main:

    // Create channels for inputs
    ch_input_file = !params.input_file ? channel.empty() : channel.fromPath(params.input_file)
    ch_metadata_csv = !params.metadata_csv ? channel.empty() : channel.fromPath(params.metadata_csv)
    ch_metadata_rds = !params.metadata_rds ? channel.empty() : channel.fromPath(params.metadata_rds)
    ch_ranked_cci_rds = channel.empty()
    ch_mvoted_rds = channel.empty()
    // Create channels
    ch_versions = channel.empty()

    // Create channels for scrnaseqcellcomm_modules
    ch_cci = channel.empty()

    def scrnaseqcellcomm_modules = params.scrnaseqcellcomm_modules ? params.scrnaseqcellcomm_modules.split(',').collect { it.trim().toLowerCase() } : []
    def avail_cci_tools = ("cell2cell,cellchat,cellphonedb,liana").split(",").collect { it.trim().toLowerCase() }
    def cci_tools = params.cci_tools ? params.cci_tools.split(',').collect { it.trim().toLowerCase() } : []

    if (scrnaseqcellcomm_modules.contains("prep_data")) {

        if (!params.seurat_obj_dir) {
            PREP_DATA(
                ch_input_file,
                params.sample_var,
                params.annot,
                params.min_cells,
                params.skip_reduction,
            )

            ch_metadata_csv = PREP_DATA.out.metadata_csv
            ch_metadata_rds = PREP_DATA.out.metadata_rds
            mtx_dir_prepped = PREP_DATA.out.mtx_dir
            seurat_obj_prepped = PREP_DATA.out.seurat_obj
            ch_versions = ch_versions.mix(PREP_DATA.out.versions)
        }
        else {
            PREP_DATA_ALT(
                params.seurat_obj_dir,
                params.sample_var,
                params.annot,
                params.min_cells,
            )
            ch_metadata_csv = PREP_DATA_ALT.out.metadata_csv
            ch_metadata_rds = PREP_DATA_ALT.out.metadata_rds
            mtx_dir_prepped = PREP_DATA_ALT.out.mtx_dir
            seurat_obj_prepped = PREP_DATA_ALT.out.seurat_obj
            ch_versions = ch_versions.mix(PREP_DATA_ALT.out.versions)
        }
    }

    if (scrnaseqcellcomm_modules.contains("run_cci")) {
        RUN_CCI(
            mtx_dir_prepped,
            seurat_obj_prepped,
            ch_metadata_csv,
        )
        ch_cci = RUN_CCI.out.cci
        ch_versions = ch_versions.mix(RUN_CCI.out.versions)
    }

    // All CCI tools need to be run for the consensus and aggregation
    if (cci_tools.intersect(avail_cci_tools).size() == 4) {
        if (!scrnaseqcellcomm_modules.contains("run_cci") && scrnaseqcellcomm_modules.contains("consensus")) {
            ch_cci = channel.fromPath(params.sample_sheet)
                .splitCsv(header: true)
                .map { row ->
                    def sid = row.sample_id
                    tuple(
                        [sample_id: sid],
                        file("${params.interactions}/02_run_cci/01_cell2cell/02_formatted/cell2cell__${sid}__postproc.rds"),
                        file("${params.interactions}/02_run_cci/02_cellchat/02_formatted/cellchat__${sid}__postproc.rds"),
                        file("${params.interactions}/02_run_cci/03_cellphonedb/02_formatted/cpdb__${sid}__postproc.rds"),
                        file("${params.interactions}/02_run_cci/04_liana/02_formatted/liana__${sid}__postproc.rds"),
                    )
                }
        }

        if (scrnaseqcellcomm_modules.contains("consensus")) {
            CONSENSUS(
                ch_cci,
                ch_metadata_rds,
                params.alpha,
                params.n_perm,
                params.condition_var,
                params.sample_var,
                params.patient_var,
            )
            ch_versions = ch_versions.mix(CONSENSUS.out.versions)
            ch_ranked_cci_rds = CONSENSUS.out.ranked_rds
            ch_mvoted_rds = CONSENSUS.out.mvoted_rds
        }



        if (scrnaseqcellcomm_modules.contains("aggregation")) {
            // Required by aggregation, either generated above or supplied as params
            if (!scrnaseqcellcomm_modules.contains("consensus")) {
                ch_ranked_cci_rds = !params.ranked_cci_rds ? channel.empty() : channel.fromPath(params.ranked_cci_rds)
                ch_mvoted_rds = !params.mvoted_rds ? channel.empty() : channel.fromPath(params.mvoted_rds)
            }


            AGGREGATION(
                ch_ranked_cci_rds,
                ch_mvoted_rds,
                params.min_patients,
            )
            ch_versions = ch_versions.mix(AGGREGATION.out.versions)
        }
    }
    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'scrnaseqcellcomm_software_' + 'versions.yml',
            sort: true,
            newLine: true,
        )
        .set { ch_collated_versions }

    emit:
    versions = ch_versions // channel: [ path(versions.yml) ]
}
