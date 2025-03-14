include { UTILS_CREATE_SAMPLESHEET         } from '../../../modules/local/utils/create_samplesheet'
include { UTILS_EXTRACT_METADATA           } from '../../../modules/local/utils/extract_metadata'
include { SEURAT_REDUCE_OBJECT_SIZE        } from '../../../modules/local/seurat/reduce_object_size'
include { SEURAT_SPLIT_OBJECT_INTO_SAMPLES } from '../../../modules/local/seurat/split_object_into_samples'
include { SEURAT_PREPROCESS_SAMPLE         } from '../../../modules/local/seurat/preprocess_sample'
include { SEURAT_EXTRACT_SAMPLE            } from '../../../modules/local/seurat/extract_sample'
include { SEURAT_SUBSET_OBJECT             } from '../../../modules/local/seurat/subset_object'
workflow PREP_DATA {
    take:
    input_file
    sample_var
    annot
    min_cells
    skip_reduction

    main:
    seurat_obj = Channel.empty()
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
        | map { row -> [sample_id: row.Sample] }

    if (!skip_reduction) {
        SEURAT_REDUCE_OBJECT_SIZE(input_file)
        ch_versions = ch_versions.mix(SEURAT_REDUCE_OBJECT_SIZE.out.versions)

        SEURAT_SUBSET_OBJECT(SEURAT_REDUCE_OBJECT_SIZE.out.rds, sample_var, UTILS_CREATE_SAMPLESHEET.out.csv)
        seurat_obj = SEURAT_SUBSET_OBJECT.out.rds
        ch_versions = ch_versions.mix(SEURAT_SUBSET_OBJECT.out.versions)
    }
    else {
        seurat_obj = input_file
    }


    // Only preprocess samples that have at least 2 cell types each having at least min_cells
    sample_sheet.combine(seurat_obj).set { samples }

    SEURAT_EXTRACT_SAMPLE(samples, sample_var)
    ch_versions = ch_versions.mix(SEURAT_EXTRACT_SAMPLE.out.versions.first())

    SEURAT_PREPROCESS_SAMPLE(
        SEURAT_EXTRACT_SAMPLE.out.rds,
        annot,
        min_cells,
    )
    ch_versions = ch_versions.mix(SEURAT_PREPROCESS_SAMPLE.out.versions.first())

    emit:
    metadata_csv = UTILS_EXTRACT_METADATA.out.csv
    metadata_rds = UTILS_EXTRACT_METADATA.out.rds
    mtx_dir      = SEURAT_PREPROCESS_SAMPLE.out.tsv.join(SEURAT_PREPROCESS_SAMPLE.out.mtx)
    seurat_obj   = SEURAT_PREPROCESS_SAMPLE.out.rds
    versions     = ch_versions
}
