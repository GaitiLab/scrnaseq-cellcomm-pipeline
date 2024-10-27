// include { EXTRACT_METADATA; REDUCE_SEURAT_OBJECT_SIZE; SAMPLE_PREPROCESSING; split_seurat_object_into_samples } from "../nf-modules/prep_data.nf"
include { EXTRACT_METADATA } from '../../../modules/local/extractmetadata.nf'
include { REDUCE_SEURAT_OBJECT_SIZE} from '../../../modules/local/reduceseuratobjectsize.nf'
include {SAMPLE_PREPROCESSING} from '../../../modules/local/samplepreprocessing.nf'
include { SPLIT_SEURAT_OBJECT_INTO_SAMPLES } from '../../../modules/local/splitseuratobjectintosamples.nf'


workflow PREP_DATA {
    take:
        input_file

    main:


    EXTRACT_METADATA(
        input_file          = input_file
    )

    ch_input_file = !params.skip_reduction ? Channel.empty() : Channel.value(input_file)

    ch_input_file.view()

    ch_input_file.ifEmpty(input_file) | REDUCE_SEURAT_OBJECT_SIZE

    seurat_obj = REDUCE_SEURAT_OBJECT_SIZE.out.rds.collect().ifEmpty(input_file)

    seurat_obj.view()

    SPLIT_SEURAT_OBJECT_INTO_SAMPLES(
        input_file          = seurat_obj,
        sample_var          = params.sample_var
    )


    SAMPLE_PREPROCESSING(
        input_file          = SPLIT_SEURAT_OBJECT_INTO_SAMPLES.out.rds
                                .flatten()
                                .map(file -> tuple(file.simpleName, file)),
        annot               = params.annot,
        min_cells           = params.min_cells,
        is_confident        = params.is_confident
    )

    mtx = SAMPLE_PREPROCESSING.out.tsv
        .join(  SAMPLE_PREPROCESSING.out.mtx  )

    mtx.view()

    emit:
    metadata_csv            = EXTRACT_METADATA.out.csv
    metadata_rds            = EXTRACT_METADATA.out.rds
    mtx_dir                 = mtx
    seurat_obj              = SAMPLE_PREPROCESSING.out.rds
}
