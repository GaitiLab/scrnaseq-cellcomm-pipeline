include { UTILS_CREATE_SAMPLESHEET       } from '../../../modules/local/utils/create_samplesheet'
include { UTILS_COMBINE_SAMPLES_METADATA } from '../../../modules/local/utils/combine_samples_metadata'
include { SEURAT_PREPROCESS_SAMPLE       } from '../../../modules/local/seurat/preprocess_sample'
include { SEURAT_CONVERT_TO_MTX          } from '../../../modules/local/seurat/convert_seurat_to_mtx'

workflow PREP_DATA_ALT {
    take:
    seurat_obj_dir
    sample_var
    annot
    min_cells

    main:
    ch_seurat_obj = channel.fromPath("${seurat_obj_dir}/*.rds")
    ch_versions = channel.empty()

    UTILS_COMBINE_SAMPLES_METADATA(ch_seurat_obj.collect())
    ch_versions = ch_versions.mix(UTILS_COMBINE_SAMPLES_METADATA.out.versions)

    UTILS_CREATE_SAMPLESHEET(
        UTILS_COMBINE_SAMPLES_METADATA.out.rds,
        sample_var,
        annot,
        min_cells,
    )
    ch_versions = ch_versions.mix(UTILS_CREATE_SAMPLESHEET.out.versions)

    sampleIds = UTILS_CREATE_SAMPLESHEET.out.csv | splitCsv(header: true) | map { row -> row.sample_id }.collect()

    ch_samples_all = ch_seurat_obj
        .combine(sampleIds)
        .map { file, ref_sample_id -> [[ref_id: ref_sample_id, sample_id: file.simpleName], file] }

    // Only preprocess samples that have at least 2 cell types each having at least min_cells
    // By using combine, for each sample_id you couple the seurat object (tuple) -> returns list of tuples
    ch_samples_all
        .filter { meta, _file -> meta.ref_id == meta.sample_id }
        .set { ch_samples }

    SEURAT_PREPROCESS_SAMPLE(
        ch_samples,
        annot,
        min_cells,
    )
    ch_versions = ch_versions.mix(SEURAT_PREPROCESS_SAMPLE.out.versions.first())

    SEURAT_CONVERT_TO_MTX(SEURAT_PREPROCESS_SAMPLE.out.rds)
    ch_versions = ch_versions.mix(SEURAT_CONVERT_TO_MTX.out.versions.first())

    emit:
    metadata_csv = UTILS_COMBINE_SAMPLES_METADATA.out.csv
    metadata_rds = UTILS_COMBINE_SAMPLES_METADATA.out.rds
    mtx_dir      = SEURAT_CONVERT_TO_MTX.out.tsv.join(SEURAT_CONVERT_TO_MTX.out.mtx)
    seurat_obj   = SEURAT_PREPROCESS_SAMPLE.out.rds
    versions     = ch_versions
}
