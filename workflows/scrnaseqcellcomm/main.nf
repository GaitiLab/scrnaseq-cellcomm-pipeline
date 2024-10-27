/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { PREP_DATA } from '../../subworkflows/local/prep_data/main.nf'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow SCRNASEQCELLCOMM {
    take:
        input_file
        metadata_csv
        metadata_rds
        annot
        min_cells
        sample_var
        cellphone_db
        cellchat_db
        liana_db
        cell2cell_db
        ref_db
        min_pct
        n_perm
        condition_var
        patient_var
        is_confident
        min_patients
        alpha
        interactions_excel_name
        skip_reduction

    main:

    def scrnaseqcellcomm_modules = params.scrnaseqcellcomm_modules ? params.scrnaseqcellcomm_modules.split(',').collect{ it.trim().toLowerCase() } : []

    if(scrnaseqcellcomm_modules.contains("prep_data")) {
        PREP_DATA(
            input_file                      = input_file,
        )

        PREP_DATA.out.metadata_csv.set  { metadata_csv }
        PREP_DATA.out.metadata_rds.set  { metadata_rds }
        PREP_DATA.out.mtx_dir.set       { preprocessing_mtx_dir }
        PREP_DATA.out.seurat_obj.set    { preprocessing_seurat_obj }



    }

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
