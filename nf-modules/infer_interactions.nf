process INFER_CELLCHAT {
    label 'mem64'
    label 'time_4h'

    publishDir "${projectDir}/output/${params.output_run_name}", mode: "copy"

    input:
    path input_file
    path interactions_db
    val annot
    val n_perm

    output:
    tuple val(input_file.simpleName), path("200_cci_cellchat/cellchat__${input_file.simpleName}.rds"), emit: cellchat_obj
    path "200_cci_cellchat/cellchat__${input_file.simpleName}__raw_obj.rds"

    script:
    def time_out_limit = (task.time).toSeconds() - 30
    """
    #!/usr/bin/env bash
    timeout ${time_out_limit} Rscript "${projectDir}/scripts/200_cci_cellchat.R" \
    --gene_expr \$PWD/${input_file} \
    --output_dir "\$PWD/200_cci_cellchat" \
    --resource \$PWD/${interactions_db} \
    --ident_col ${annot} \
    --n_perm ${n_perm} \
    --n_cores ${task.cpus}

    """
}

process INFER_LIANA {
    label 'mem4'
    label 'time_15m'

    publishDir "${projectDir}/output/${params.output_run_name}/", mode: "copy"

    input:
    path input_file
    path interactions_db
    val annot
    val n_perm

    output:

    tuple val(input_file.simpleName), path("201_cci_liana/liana__${input_file.simpleName}.rds"), emit: liana_obj

    script:
    def time_out_limit = (task.time).toSeconds() - 30
    """
    #!/usr/bin/env bash
    timeout ${time_out_limit} Rscript "${projectDir}/scripts/201_cci_liana.R" \
    --gene_expr \$PWD/${input_file} \
    --output_dir "\$PWD/201_cci_liana" \
    --resource \$PWD/${interactions_db} \
    --ident_col ${annot} \
    --n_perm ${n_perm}

    """
}

process INFER_CELL2CELL {
    label 'mem16'
    label 'time_2h'

    publishDir "${projectDir}/output/${params.output_run_name}/", mode: "copy"

    input:
    path input_dir
    path meta
    path interactions_db
    val annot
    val n_perm

    output:
    path "202_cci_cell2cell/cell2cell__${input_dir.simpleName}.pickle"
    tuple val(input_dir.simpleName), path("202_cci_cell2cell/cell2cell__${input_dir.simpleName}.csv"), emit: cell2cell_obj

    script:
    def time_out_limit = (task.time).toSeconds() - 30
    """
    #!/usr/bin/env bash

    timeout ${time_out_limit} python3 "${projectDir}/Python/202_cci_cell2cell.py" \
    --input_dir \$PWD/${input_dir} \
    --output_dir "\$PWD/202_cci_cell2cell" \
    --annot ${annot} \
    --meta \$PWD/${meta} \
    --interactions \$PWD/${interactions_db} \
    --nperm ${n_perm} \
    --sample_id ${input_dir.simpleName}

    """
}

process INFER_CPDB {
    label 'cpdb_env'
    label 'mem16'
    label 'time_1h'

    publishDir "${projectDir}/output/${params.output_run_name}/", mode: "copy"

    input:
    path input_dir
    path meta
    path interactions_db
    val annot
    val n_perm
    val min_pct
    val alpha

    output:

    tuple val(input_dir.simpleName), 
    path("203_cci_cpdb/statistical_analysis_interaction_scores__${input_dir.simpleName}.txt"), 
    path("203_cci_cpdb/statistical_analysis_pvalues__${input_dir.simpleName}.txt"), 
    path("203_cci_cpdb/statistical_analysis_significant_means__${input_dir.simpleName}.txt"), 
    path("203_cci_cpdb/statistical_analysis_deconvoluted__${input_dir.simpleName}.txt"), 
    path("203_cci_cpdb/statistical_analysis_deconvoluted_percents__${input_dir.simpleName}.txt"), 
    path("203_cci_cpdb/statistical_analysis_means__${input_dir.simpleName}.txt"), emit: cpdb_obj
    path "203_cci_cpdb/${input_dir.simpleName}_counts.h5ad"
    path "203_cci_cpdb/${input_dir.simpleName}_metadata.tsv"

    script:
    def time_out_limit = (task.time).toSeconds() - 30
    """
    #!/usr/bin/env bash
    timeout ${time_out_limit} python3 "${projectDir}/Python/203_cci_cpdb.py" \
    --output_dir \$PWD/203_cci_cpdb \
    --meta \$PWD/$meta \
    --annot $annot \
    --cpdb_file_path \$PWD/$interactions_db \
    --n_perm $n_perm \
    --input_dir \$PWD/$input_dir \
    --min_pct $min_pct \
    --threads ${task.cpus} \
    --alpha $alpha \
    --sample_id ${input_dir.simpleName}

    """

}