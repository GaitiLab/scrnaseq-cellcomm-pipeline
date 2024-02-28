process POSTPROCESSING_CELLCHAT {
    label 'time_10m'
    label 'mem2'

    publishDir "${projectDir}/output/${params.run_name}", mode: "copy"

    input:
    tuple val(sample_id), val(full_id), path(input_interactions)
    path ref_db


    output:
    tuple val(sample_id), path("300_postproc_cellchat/cellchat__${full_id}__postproc.rds")

    script:
    def time_out_limit = (task.time).toSeconds() - 30

    """
    #!/usr/bin/env bash

    timeout ${time_out_limit} Rscript "${projectDir}/scripts/300_postproc_cellchat.R" \
    --output_dir "\$PWD/300_postproc_cellchat/" \
    --input_interactions \$PWD/$input_interactions \
    --ref_db \$PWD/${ref_db} \
    --sample_id ${full_id}
    """
}

process POSTPROCESSING_LIANA {
    label 'time_10m'
    label 'mem2'

    publishDir "${projectDir}/output/${params.run_name}", mode: "copy"

    input:
    tuple val(sample_id), val(full_id), path(input_interactions)
    path ref_db 

    output:
    tuple val(sample_id), path("301_postproc_liana/liana__${full_id}__postproc.rds")

    script:
    def time_out_limit = (task.time).toSeconds() - 30

    """
    #!/usr/bin/env bash

    timeout ${time_out_limit} Rscript "${projectDir}/scripts/301_postproc_liana.R" \
    --output_dir "\$PWD/301_postproc_liana/" \
    --input_interactions \$PWD/$input_interactions \
    --sample_id ${full_id} \
    --ref_db \$PWD/${ref_db}
    """
}


process POSTPROCESSING_CELL2CELL {
    label 'time_10m'
    label 'mem2'

    publishDir "${projectDir}/output/${params.run_name}", mode: "copy"

    input:
    tuple val(sample_id), val(full_id), path(input_interactions)
    path ref_db

    output:
    tuple val(sample_id), path("302_postproc_cell2cell/cell2cell__${full_id}__postproc.rds")

    script:
    def time_out_limit = (task.time).toSeconds() - 30

    """
    #!/usr/bin/env bash

    timeout ${time_out_limit} Rscript "${projectDir}/scripts/302_postproc_cell2cell.R" \
    --output_dir "\$PWD/302_postproc_cell2cell/" \
    --input_interactions \$PWD/$input_interactions \
    --sample_id ${full_id} \
    --ref_db \$PWD/${ref_db}
    """
}

process POSTPROCESSING_CPDB {
    label 'time_10m'
    label 'mem2'

    publishDir "${projectDir}/output/${params.run_name}", mode: "copy"

    input:
    tuple val(sample_id), val(full_id), path(interaction_scores), 
    path(pvalues), path(significant_means), 
    path(deconvoluted), 
    path(deconvoluted_percents), 
    path(means)
    path ref_db

    output:
    tuple val(sample_id), path("303_postproc_cpdb/cpdb__${full_id}__postproc.rds")

    script:
    def time_out_limit = (task.time).toSeconds() - 30

    """
    #!/usr/bin/env bash

    timeout ${time_out_limit} Rscript "${projectDir}/scripts/303_postproc_cellphonedb.R" \
    --output_dir "\$PWD/303_postproc_cpdb/" \
    --sample_id ${full_id} \
    --interaction_scores \$PWD/${interaction_scores} \
    --pval \$PWD/${pvalues} \
    --sign_means \$PWD/${significant_means} \
    --means \$PWD/${means} \
    --ref_db \$PWD/${ref_db}

    """
}



process COMBINING_CELLCHAT_RUNS {
    label 'time_10m'
    label 'mem2'

    publishDir "${projectDir}/output/${params.run_name}", mode: "copy"

    input:
    tuple val(sample_id), path(postprocessed_file)
    val n_repeats

    output:
    tuple val(sample_id), path("300_postproc_cellchat/cellchat__${sample_id}__postproc.rds")

    script:
    def time_out_limit = (task.time).toSeconds() - 30

    """
    #!/usr/bin/env bash

    timeout ${time_out_limit} Rscript "${projectDir}/scripts/304_combine_cellchat.R" \
    --output_dir "\$PWD/300_postproc_cellchat" \
    --sample_id ${sample_id} \
    --input_dir "\$PWD" \
    --n_repeats ${n_repeats}

    """
}

process COMBINING_LIANA_RUNS {
    label 'time_10m'
    label 'mem2'

    publishDir "${projectDir}/output/${params.run_name}", mode: "copy"

    input:
    tuple val(sample_id), path(postprocessed_file)
    val n_repeats

    output:
    tuple val(sample_id), path("301_postproc_liana/liana__${sample_id}__postproc.rds")

    script:
    def time_out_limit = (task.time).toSeconds() - 30

    """
    #!/usr/bin/env bash

    timeout ${time_out_limit} Rscript "${projectDir}/scripts/305_combine_liana.R" \
    --output_dir "\$PWD/301_postproc_liana" \
    --sample_id ${sample_id} \
    --input_dir "\$PWD" \
    --n_repeats ${n_repeats}


    """
}

process COMBINING_CPDB_RUNS {
    label 'time_10m'
    label 'mem2'

    publishDir "${projectDir}/output/${params.run_name}", mode: "copy"

    input:
    tuple val(sample_id), path(postprocessed_file)
    val n_repeats 

    output:
    tuple val(sample_id), path("303_postproc_cpdb/cpdb__${sample_id}__postproc.rds")

    script:
    def time_out_limit = (task.time).toSeconds() - 30

    """
    #!/usr/bin/env bash

    timeout ${time_out_limit} Rscript "${projectDir}/scripts/307_combine_cpdb.R" \
    --output_dir "\$PWD/303_postproc_cpdb" \
    --sample_id ${sample_id} \
    --input_dir "\$PWD" \
    --n_repeats ${n_repeats}

    """
}


process COMBINING_CELL2CELL_RUNS {
    label 'time_10m'
    label 'mem2'

    publishDir "${projectDir}/output/${params.run_name}", mode: "copy"

    input:
    tuple val(sample_id), path(postprocessed_file)
    val n_repeats

    output:
    tuple val(sample_id), path("302_postproc_cell2cell/cell2cell__${sample_id}__postproc.rds")

    script:
    def time_out_limit = (task.time).toSeconds() - 30

    """
    #!/usr/bin/env bash

    timeout ${time_out_limit} Rscript "${projectDir}/scripts/306_combine_cell2cell.R" \
    --output_dir "\$PWD/302_postproc_cell2cell" \
    --sample_id ${sample_id} \
    --input_dir "\$PWD" \
    --n_repeats ${n_repeats}


    """
}