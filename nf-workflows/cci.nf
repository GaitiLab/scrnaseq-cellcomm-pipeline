include { INFER_CCI } from "../nf-subworkflows/infer_cci.nf"
include { POSTPROCESSING_CCI } from "../nf-subworkflows/postprocessing_cci.nf"

workflow CCI_SIMPLE_PIPELINE {
    take:
    preprocessing_seurat_obj   
    preprocessing_mtx_dir     
    metadata_csv                
    cellchat_db
    cell2cell_db
    liana_db
    cellphone_db
    ref_db
    output_dir

    main:
    cci_liana = Channel.empty()
    cci_cell2cell = Channel.empty()
    cci_cpdb = Channel.empty()
    cci_cellchat = Channel.empty()
    
    if (params.init_step == 2) {
        INFER_CCI( 
            preprocessing_seurat_obj    = preprocessing_seurat_obj,
            preprocessing_mtx_dir       = preprocessing_mtx_dir,
            metadata_csv                = metadata_csv,
            cellchat_db                 = cellchat_db,
            cell2cell_db                = cell2cell_db,
            liana_db                    = liana_db,
            cellphone_db                = cellphone_db
            )

            INFER_CCI.out.cci_liana.set     {   cci_liana       }
            INFER_CCI.out.cci_cell2cell.set {   cci_cell2cell   }
            INFER_CCI.out.cci_cpdb.set      {   cci_cpdb        }
            INFER_CCI.out.cci_cellchat.set  {   cci_cellchat    }
    } 
    if (params.init_step == 3) {
        Channel.fromPath("${output_dir}/201_cci_liana/liana__*.rds", type: 'file') \
            | map { file -> def sample_id = input_file.simpleName.split( '__' )[0]
                return tuple(sample_id, file)
                }. set {   cci_liana    }
        Channel.fromPath("${output_dir}/200_cci_cellchat/^cellchat__*((?!raw_obj).)*.rds", type: 'file') \
            | map { file -> def sample_id = input_file.simpleName.split( '__' )[0]
                        return tuple(sample_id, file)
                        }.set { cci_cellchat }
        Channel.fromPath("${output_dir}/202_cci_cell2cell/cell2cell__*.csv", type: 'file') |
            | map { file -> def sample_id = input_file.simpleName.split( '__' )[0]
                            return tuple(sample_id, file)
                        }.set { cci_cell2cell }
        Channel.fromPath("${output_dir}/203_cci_cpdb/statistical_analysis_means__*.txt", type: 'file') \
            | map { file -> def sample_id = input_file.simpleName.split( '__' )[0]
                    return tuple(sample_id, file)
                        }.set { cci_cpdb }
    }

    POSTPROCESSING_CCI( 
        liana_obj                   = cci_liana,
        cell2cell_obj               = cci_cell2cell,
        cpdb_obj                    = cci_cpdb, 
        cellchat_obj                = cci_cellchat,
        ref_db                      = ref_db
    )
    

    emit: 
    matched_cci                     = POSTPROCESSING_CCI.out.cci_cellchat
                                        .join(
                                        POSTPROCESSING_CCI.out.cci_liana, by: 0)
                                        .join(
                                        POSTPROCESSING_CCI.out.cci_cell2cell, by: 0)
                                        .join(
                                        POSTPROCESSING_CCI.out.cci_cpdb, by: 0)
}