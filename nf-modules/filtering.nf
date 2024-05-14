process POSTPROCESSING_CELLCHAT {
    label 'time_10m'
    label 'mem_4G'

    publishDir params.output_dir, mode: "copy"

    input:
    tuple val(sample_id), path(input_interactions)
    path ref_db

    output:
    tuple val(sample_id), path("300_postproc_cellchat/cellchat__${sample_id}__postproc.rds")

    script:
    """
    #!/usr/bin/env bash

    Rscript "${projectDir}/scripts/300_postproc_cellchat.R" \
    --output_dir "\$PWD/300_postproc_cellchat/" \
    --input_interactions \$PWD/$input_interactions \
    --ref_db \$PWD/${ref_db} \
    --sample_id ${sample_id}
    """

    stub:
    """
    #!/usr/bin/env bash
    mkdir -p 300_postproc_cellchat
    touch "300_postproc_cellchat/cellchat__${sample_id}__postproc.rds"
    """    
}

process POSTPROCESSING_LIANA {
    label 'time_10m'
    label 'mem_4G'

    publishDir params.output_dir, mode: "copy"

    input:
    tuple val(sample_id), path(input_interactions)
    path ref_db 

    output:
    tuple val(sample_id), path("301_postproc_liana/liana__${sample_id}__postproc.rds")

    script:

    """
    #!/usr/bin/env bash

    Rscript "${projectDir}/scripts/301_postproc_liana.R" \
    --output_dir "\$PWD/301_postproc_liana/" \
    --input_interactions \$PWD/$input_interactions \
    --sample_id ${sample_id} \
    --ref_db \$PWD/${ref_db}
    """

    stub:
    """
    #!/usr/bin/env bash
    mkdir -p 301_postproc_liana
    touch "301_postproc_liana/liana__${sample_id}__postproc.rds"
    """ 
}


process POSTPROCESSING_CELL2CELL {
    label 'time_10m'
    label 'mem_4G'

    publishDir params.output_dir, mode: "copy"

    input:
    tuple val(sample_id), path(input_interactions)
    path ref_db

    output:
    tuple val(sample_id), path("302_postproc_cell2cell/cell2cell__${sample_id}__postproc.rds")

    script:
    """
    #!/usr/bin/env bash

    Rscript "${projectDir}/scripts/302_postproc_cell2cell.R" \
    --output_dir "\$PWD/302_postproc_cell2cell/" \
    --input_interactions \$PWD/$input_interactions \
    --sample_id ${sample_id} \
    --ref_db \$PWD/${ref_db}
    """

    stub:
    """
    #!/usr/bin/env bash
    mkdir -p 302_postproc_cell2cell
    touch "302_postproc_cell2cell/cell2cell__${sample_id}__postproc.rds"
    """ 
}

process POSTPROCESSING_CPDB {
    label 'time_10m'
    label 'mem_4G'

    publishDir params.output_dir, mode: "copy"

    input:
    tuple val(sample_id), path(interaction_scores), 
    path(pvalues), path(significant_means), path(means)
    path ref_db

    output:
    tuple val(sample_id), path("303_postproc_cpdb/cpdb__${sample_id}__postproc.rds")

    script:
    """
    #!/usr/bin/env bash

    Rscript "${projectDir}/scripts/303_postproc_cellphonedb.R" \
    --output_dir "\$PWD/303_postproc_cpdb/" \
    --sample_id ${sample_id} \
    --interaction_scores \$PWD/${interaction_scores} \
    --pval \$PWD/${pvalues} \
    --sign_means \$PWD/${significant_means} \
    --means \$PWD/${means} \
    --ref_db \$PWD/${ref_db}

    """

    stub:
    """
    #!/usr/bin/env bash
    mkdir -p 303_postproc_cellphonedb
    touch "303_postproc_cellphonedb/cpdb__${sample_id}__postproc.rds"
    """ 
}