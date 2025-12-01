include { UTILS_CREATE_SAMPLESHEET  } from '../../../modules/local/utils/create_samplesheet'
include { UTILS_EXTRACT_METADATA    } from '../../../modules/local/utils/extract_metadata'
include { SEURAT_REDUCE_OBJECT_SIZE } from '../../../modules/local/seurat/reduce_object_size'
include { SEURAT_PREPROCESS_SAMPLE  } from '../../../modules/local/seurat/preprocess_sample'
include { SEURAT_EXTRACT_SAMPLE     } from '../../../modules/local/seurat/extract_sample'
include { SEURAT_SUBSET_OBJECT      } from '../../../modules/local/seurat/subset_object'
include { SEURAT_CONVERT_TO_MTX     } from '../../../modules/local/seurat/convert_seurat_to_mtx'

workflow PREP_DATA {
    take:
    input_file
    sample_var
    annot
    min_cells
    skip_reduction

    main:
    ch_seurat_obj = Channel.empty()
    ch_versions = Channel.empty()

    UTILS_EXTRACT_METADATA(input_file)
    ch_versions = ch_versions.mix(UTILS_EXTRACT_METADATA.out.versions)

    UTILS_CREATE_SAMPLESHEET(
        UTILS_EXTRACT_METADATA.out.rds,
        sample_var,
        annot,
        min_cells,
    )
    ch_versions = ch_versions.mix(UTILS_CREATE_SAMPLESHEET.out.versions)

    sample_sheet = UTILS_CREATE_SAMPLESHEET.out.csv
        | splitCsv(header: true)
        | map { row -> [sample_id: row.sample_id] }

    if (!skip_reduction) {
        SEURAT_REDUCE_OBJECT_SIZE(input_file)
        ch_versions = ch_versions.mix(SEURAT_REDUCE_OBJECT_SIZE.out.versions)

        SEURAT_SUBSET_OBJECT(SEURAT_REDUCE_OBJECT_SIZE.out.rds, sample_var, UTILS_CREATE_SAMPLESHEET.out.csv)
        ch_seurat_obj = SEURAT_SUBSET_OBJECT.out.rds
        ch_versions = ch_versions.mix(SEURAT_SUBSET_OBJECT.out.versions)
    }
    else {
        ch_seurat_obj = input_file
    }

    // Only preprocess samples that have at least 2 cell types each having at least min_cells
    // By using combine, for each sample_id you couple the seurat object (tuple) -> returns list of tuples
    sample_sheet.combine(ch_seurat_obj).set { ch_samples }

    SEURAT_EXTRACT_SAMPLE(ch_samples, sample_var)
    ch_versions = ch_versions.mix(SEURAT_EXTRACT_SAMPLE.out.versions.first())

    SEURAT_PREPROCESS_SAMPLE(
        SEURAT_EXTRACT_SAMPLE.out.rds,
        annot,
        min_cells,
    )
    ch_versions = ch_versions.mix(SEURAT_PREPROCESS_SAMPLE.out.versions.first())

    SEURAT_CONVERT_TO_MTX(SEURAT_PREPROCESS_SAMPLE.out.rds)
    ch_versions = ch_versions.mix(SEURAT_CONVERT_TO_MTX.out.versions.first())

    emit:
    metadata_csv = UTILS_EXTRACT_METADATA.out.csv
    metadata_rds = UTILS_EXTRACT_METADATA.out.rds
    mtx_dir      = SEURAT_CONVERT_TO_MTX.out.tsv.join(SEURAT_CONVERT_TO_MTX.out.mtx)
    seurat_obj   = SEURAT_PREPROCESS_SAMPLE.out.rds
    versions     = ch_versions
}
