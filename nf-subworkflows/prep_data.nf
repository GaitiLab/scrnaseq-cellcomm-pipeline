include { extract_metadata; reduce_seurat_object_size; sample_preprocessing; split_seurat_object_into_samples } from "../nf-modules/prep_data.nf"

workflow PREP_DATA {
    take: 
        input_file

    main:
    extract_metadata(
        input_file          = input_file
    )

    if (!params.skip_reduction) {
        reduce_seurat_object_size(
            input_file      = input_file, 
            annot           = params.annot
        ).set { input_file }
    }

    split_seurat_object_into_samples(
        input_file          = input_file, 
        sample_var          = params.sample_var
    )

    sample_preprocessing(
        input_file          = split_seurat_object_into_samples.out
                                .flatten()
                                .map(file -> tuple(file.simpleName, file)), 
        annot               = params.annot,
        min_cells           = params.min_cells,
        is_confident        = params.is_confident
    )
    emit:  
    metadata_csv            = extract_metadata.out.metadata_csv
    metadata_rds            = extract_metadata.out.metadata_rds
    mtx_dir                 = sample_preprocessing.out.mtx_dir
    seurat_obj              = sample_preprocessing.out.seurat_obj
}