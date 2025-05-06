/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { paramsSummaryMap       } from 'plugin/nf-schema'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_scrnaseqcellcomm_pipeline'

include { PREP_DATA              } from '../subworkflows/local/prep_data'
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
    ch_input_file = Channel.fromPath(params.input_file)
    ch_metadata_csv = !params.metadata_csv ? Channel.empty() : Channel.fromPath(params.metadata_csv)
    ch_metadata_rds = !params.metadata_rds ? Channel.empty() : Channel.fromPath(params.metadata_rds)

    // Create channels
    ch_versions = Channel.empty()

    // Create channels for scrnaseqcellcomm_modules
    ch_cci = channel.empty()
    ch_consensus = Channel.empty()

    def scrnaseqcellcomm_modules = params.scrnaseqcellcomm_modules ? params.scrnaseqcellcomm_modules.split(',').collect { it.trim().toLowerCase() } : []
    def avail_cci_tools = ("cell2cell,cellchat,cellphonedb,liana").split(",").collect { it.trim().toLowerCase() }
    def cci_tools = params.cci_tools ? params.cci_tools.split(',').collect { it.trim().toLowerCase() } : []


    if (scrnaseqcellcomm_modules.contains("prep_data")) {
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
            ch_consensus = CONSENSUS.out.rds
        }

        if (scrnaseqcellcomm_modules.contains("aggregation")) {
            AGGREGATION(
                ch_consensus,
                params.condition_var,
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
