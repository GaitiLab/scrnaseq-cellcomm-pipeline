include { EXTRACT_METADATA                 } from '../../../modules/local/extractmetadata.nf'
include { REDUCE_SEURAT_OBJECT_SIZE        } from '../../../modules/local/reduceseuratobjectsize.nf'
include { SAMPLE_PREPROCESSING             } from '../../../modules/local/samplepreprocessing.nf'
include { SPLIT_SEURAT_OBJECT_INTO_SAMPLES } from '../../../modules/local/splitseuratobjectintosamples.nf'
include { CREATE_SAMPLE_SHEET              } from '../../../modules/local/createsamplesheet.nf'

workflow PREP_DATA {
    take:
    input_file    
    sample_var    
    annot         
    min_cells     
    is_confident  
    skip_reduction

    main:
    seurat_obj = Channel.empty()

    EXTRACT_METADATA(input_file)

    CREATE_SAMPLE_SHEET(
        EXTRACT_METADATA.out.rds,
        sample_var,
        annot,
        min_cells,
        is_confident
    )

    sample_sheet = CREATE_SAMPLE_SHEET.out.csv
        | splitCsv(header: true)
        | map { row -> row.Sample }

    if (!skip_reduction) {
        REDUCE_SEURAT_OBJECT_SIZE(input_file)
        seurat_obj = REDUCE_SEURAT_OBJECT_SIZE.out.rds
    }
    else {
        seurat_obj = input_file
    }

    SPLIT_SEURAT_OBJECT_INTO_SAMPLES(
        seurat_obj,
        sample_var
    )

    SPLIT_SEURAT_OBJECT_INTO_SAMPLES.out.rds
        .flatten()
        .map { file -> [file.simpleName, file] }
        .set {
            seurat_objects
        }

    // Only preprocess samples that have at least 2 cell types each having at least min_cells
    seurat_objects = sample_sheet.join(seurat_objects)

    SAMPLE_PREPROCESSING(
        seurat_objects,
        annot,
        min_cells,
        is_confident
    )

    mtx = SAMPLE_PREPROCESSING.out.tsv.join(SAMPLE_PREPROCESSING.out.mtx)

    emit:
    metadata_csv = EXTRACT_METADATA.out.csv
    metadata_rds = EXTRACT_METADATA.out.rds
    mtx_dir      = mtx
    seurat_obj   = SAMPLE_PREPROCESSING.out.rds
}
