include { GET_METADATA; REDUCE_SEURAT_OBJECT_SIZE; PREPROCESSING; SPLIT_SEURAT_OBJECT } from "../nf-modules/prep_data.nf"

workflow PREP_DATA {
    take: 
        input_file
        liana_db

    main:
    GET_METADATA(
        input_file          = input_file
    )

    if (!params.skip_reduction) {
        REDUCE_SEURAT_OBJECT_SIZE(
            input_file      = input_file, 
            annot           = params.annot
        ).set { input_file }
    }

    SPLIT_SEURAT_OBJECT(
        input_file          = input_file, 
        split_varname       = params.split_varname
    )
    
    PREPROCESSING(
        input_file          = SPLIT_SEURAT_OBJECT.out.flatten(), 
        interactions_db     = liana_db,
        annot               = params.annot,
        min_cells           = params.min_cells,
        min_cell_types      = params.min_cell_types, 
        is_confident        = params.is_confident
    )

    emit:  
    metadata_csv            = GET_METADATA.out.metadata_csv
    metadata_rds            = GET_METADATA.out.metadata_rds
    mtx_dir                 = PREPROCESSING.out.mtx_dir.flatten()
    seurat_obj              = PREPROCESSING.out.seurat_obj.flatten()
}