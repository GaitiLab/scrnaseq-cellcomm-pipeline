process infer_cellchat {
    label 'mem_32G'
    label 'time_12h'

    publishDir params.output_dir, mode: "copy"

    input:
    tuple val(sample_id), path(input_file)
    path interactions_db
    val annot
    val n_perm
    val min_cells

    output:
    tuple val(sample_id), path("200_cci_cellchat/cellchat__${sample_id}.rds"), emit: cellchat_obj
    path "200_cci_cellchat/cellchat__${sample_id}__raw_obj.rds"

    script:
    """
    #!/usr/bin/env bash
    Rscript "${projectDir}/scripts/200_cci_cellchat.R" \
        --n_perm ${n_perm} \
        --interactions_db \$PWD/${interactions_db} \
        --annot ${annot} \
        --gene_expr \$PWD/${input_file} \
        --min_cells ${min_cells} \
        --output_dir "\$PWD/200_cci_cellchat" \
        --n_cores ${task.cpus}
    """

    stub:
    """
    #!/usr/bin/env bash
    mkdir -p 200_cci_cellchat
    touch "200_cci_cellchat/cellchat__${sample_id}.rds"
    touch "200_cci_cellchat/cellchat__${sample_id}__raw_obj.rds"
    """     
}

process infer_liana {
    label 'mem_8G'
    label 'time_30m'

    publishDir params.output_dir, mode: "copy"

    input:
    tuple val(sample_id), path(input_file)
    path interactions_db
    val annot
    val n_perm
    val min_cells
    val min_pct

    output:
    tuple val(sample_id), path("201_cci_liana/liana__${sample_id}.rds"), emit: liana_obj

    script:
    """
    #!/usr/bin/env bash
    Rscript "${projectDir}/scripts/201_cci_liana.R" \
        --n_perm ${n_perm} \
        --interactions_db \$PWD/${interactions_db} \
        --annot ${annot} \
        --gene_expr \$PWD/${input_file} \
        --min_cells ${min_cells} \
        --min_pct ${min_pct} \
        --output_dir "\$PWD/201_cci_liana"
    """

    stub:
    """
    #!/usr/bin/env bash
    mkdir -p 201_cci_liana
    touch "201_cci_liana/liana__${sample_id}.rds"
    """     
}

process infer_cell2cell {
    label 'mem_64G'
    label 'time_24h'

    publishDir params.output_dir, mode: "copy"

    input:
    tuple val(sample_id), path(barcodes), 
    path(genes), path(matrix)
    path meta
    path interactions_db
    val annot
    val n_perm

    output:
    tuple val(sample_id), path("202_cci_cell2cell/cell2cell__${sample_id}.pickle")
    tuple val(sample_id), path("202_cci_cell2cell/cell2cell__${sample_id}.csv"), emit: cell2cell_obj

    script:
    """
    #!/usr/bin/env bash
    python3 "${projectDir}/Python/202_cci_cell2cell.py" \
        --input_dir \$PWD \
        --n_perm ${n_perm} \
        --interactions_db \$PWD/${interactions_db} \
        --annot ${annot} \
        --sample_id ${sample_id} \
        --meta \$PWD/${meta} \
        --output_dir "\$PWD/202_cci_cell2cell"
    """

    stub:
    """
    #!/usr/bin/env bash
    mkdir -p 202_cci_cell2cell
    touch "202_cci_cell2cell/cell2cell__${sample_id}.pickle"
    touch "202_cci_cell2cell/cell2cell__${sample_id}.csv"
    """  
}

process infer_cpdb {
    label 'cpdb_env'
    label 'mem_16G'
    label 'time_1h'

    publishDir params.output_dir, mode: "copy"

    input:
    tuple val(sample_id), path(barcodes), 
    path(genes), path(matrix)
    path meta
    path interactions_db
    val annot
    val n_perm
    val min_pct

    output:
    tuple val(sample_id), 
    path("203_cci_cpdb/statistical_analysis_interaction_scores__${sample_id}.txt"), 
    path("203_cci_cpdb/statistical_analysis_pvalues__${sample_id}.txt"), 
    path("203_cci_cpdb/statistical_analysis_significant_means__${sample_id}.txt"),
    path("203_cci_cpdb/statistical_analysis_means__${sample_id}.txt"), emit: cpdb_obj
    path "203_cci_cpdb/statistical_analysis_deconvoluted__${sample_id}.txt"
    path "203_cci_cpdb/statistical_analysis_deconvoluted_percents__${sample_id}.txt"
    path "203_cci_cpdb/${sample_id}_counts.h5ad"
    path "203_cci_cpdb/${sample_id}_metadata.tsv"

    script:
    """
    #!/usr/bin/env bash

    mkdir -p ${sample_id}

    mv ${barcodes} ${sample_id}/barcodes.tsv
    mv ${genes} ${sample_id}/genes.tsv
    mv ${matrix} ${sample_id}/matrix.mtx

    python3 "${projectDir}/Python/203_cci_cpdb.py" \
        --input_dir \$PWD/${sample_id} \
        --n_perm ${n_perm} \
        --interactions_db \$PWD/${interactions_db} \
        --annot ${annot} \
        --sample_id ${sample_id} \
        --meta \$PWD/${meta} \
        --min_pct ${min_pct} \
        --n_cores ${task.cpus} \
        --output_dir \$PWD/203_cci_cpdb
    """

    stub:
    """
    #!/usr/bin/env bash
    mkdir -p 203_cci_cpdb
    touch "203_cci_cpdb/statistical_analysis_interaction_scores__${sample_id}.txt"
    touch "203_cci_cpdb/statistical_analysis_pvalues__${sample_id}.txt"
    touch "203_cci_cpdb/statistical_analysis_significant_means__${sample_id}.txt"
    touch "203_cci_cpdb/statistical_analysis_means__${sample_id}.txt"
    touch "203_cci_cpdb/statistical_analysis_deconvoluted__${sample_id}.txt
    touch "203_cci_cpdb/statistical_analysis_deconvoluted_percents__${sample_id}.txt"
    touch "203_cci_cpdb/${sample_id}_counts.h5ad"
    touch "203_cci_cpdb/${sample_id}_metadata.tsv"
    """      
}