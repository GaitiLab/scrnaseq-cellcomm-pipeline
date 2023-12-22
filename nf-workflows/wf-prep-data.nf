// include { ADAPT_ANNOTATION; GET_METADATA; REDUCE_SEURAT_OBJECT_SIZE; SPLIT_SEURAT_OBJECT } from "../nf-modules/prep_data.nf"

// workflow ALL_PREP_DATA {
//     take:
//         input_file

//     main:

//     ADAPT_ANNOTATION(input_file = params.input_file, samples_oi = params.samples_oi)

//     GET_METADATA(input_file = ADAPT_ANNOTATION.out)

//     REDUCE_SEURAT_OBJECT_SIZE(input_file = ADAPT_ANNOTATION.out)

//     SPLIT_SEURAT_OBJECT(input_file = REDUCE_SEURAT_OBJECT_SIZE.out, split_varname = params.split_varname)

//     PREPROCESSING(input_file = SPLIT_SEURAT_OBJECT.out.flatten(), interactions_db = params.liana_db,
//     annot = params.annot,
//     min_cells = params.min_cells,
//     min_cell_types = params.min_cell_types)

//     emit:
//     GET_METADATA.
// }

