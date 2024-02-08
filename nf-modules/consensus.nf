
process CONSENSUS {
    label "mem4"
    label "time_15m"

    publishDir "${projectDir}/output/${params.run_name}", mode: "copy"

    input:
    tuple val(sample_id), path(cellchat_obj), path(liana_obj), path(cell2cell_obj), path(cpdb_obj)
    val alpha

    output:
    path("400_consensus/${sample_id}__interactions_mvoted.rds"
    ), emit: mvoted_interactions
    path("400_consensus/${sample_id}__signif_interactions.rds"
    ), emit: signif_interactions

    script:
    def time_out_limit = (task.time).toSeconds() - 30

    """
    #!/usr/bin/env bash

    timeout ${time_out_limit} Rscript "${projectDir}/scripts/400_consensus.R" \
    --output_dir "\$PWD/400_consensus" \
    --sample_id ${sample_id} \
    --alpha ${alpha} \
    --cellchat_obj \$PWD/${cellchat_obj} \
    --liana_obj \$PWD/${liana_obj} \
    --cell2cell_obj \$PWD/${cell2cell_obj} \
    --cpdb_obj \$PWD/${cpdb_obj}
    """
}


process COMBINE_SAMPLES {
    label "mem2"
    label "time_10m"

    publishDir "${projectDir}/output/${params.run_name}", mode: "copy"

    input:
    path "*__interactions_mvoted.rds"
    path "*__signif_interactions.rds"
    path metadata
    path meta_vars_oi
    val condition_varname
    val sample_varname
    val patient_varname


    output:
    path "401_combine_samples/401_samples_interactions_mvoted.rds", emit: mvoted_interactions
    path "401_combine_samples/401_samples_sign_interactions.rds"

    script:
    def time_out_limit = (task.time).toSeconds() - 30
    """
    #!/usr/bin/env bash

    timeout ${time_out_limit} Rscript "${projectDir}/scripts/401_combine_samples.R" \
    --output_dir \$PWD/401_combine_samples \
    --input_dir \$PWD \
    --metadata \$PWD/${metadata} \
    --meta_vars_oi \$PWD/${meta_vars_oi} \
    --condition_varname ${condition_varname} \
    --sample_varname ${sample_varname} \
    --patient_varname ${patient_varname}
    """
}


process AGGREGATION_SAMPLE {
    label "mem2"
    label "time_10m"

    publishDir "${projectDir}/output/${params.run_name}", mode: "copy"

    input:
    path input_file
    path metadata
    val min_cells
    val min_frac_samples
    val annot
    val condition_varname
    val sample_varname

    output:
    path "402_aggregation/402_number_of_interactions.xlsx"
    path "402_aggregation/402_sample_interactions_mvoted_w_filters.rds"

    script:
    def time_out_limit = (task.time).toSeconds() - 30
    """
    #!/usr/bin/env bash

    timeout ${time_out_limit} Rscript "${projectDir}/scripts/402a_aggregation_sample.R" \
    --output_dir \$PWD/402_aggregation \
    --input_file \$PWD/${input_file} \
    --metadata \$PWD/${metadata} \
    --min_cells ${min_cells} \
    --min_frac_samples ${min_frac_samples} \
    --annot ${annot} \
    --condition_varname ${condition_varname} \
    --sample_varname ${sample_varname}
    """
}

process AGGREGATION_PATIENT {
    label "mem2"
    label "time_10m"

    publishDir "${projectDir}/output/${params.run_name}", mode: "copy"

    input:
    path input_file
    val annot
    val condition_varname
    val min_patients

    output:
    path "402_aggregation/402_patient_interactions_mvoted_w_filters.rds"

    script:
    def time_out_limit = (task.time).toSeconds() - 30
    """
    #!/usr/bin/env bash

    timeout ${time_out_limit} Rscript "${projectDir}/scripts/402b_aggregation_patient.R" \
    --output_dir \$PWD/402_aggregation \
    --input_file \$PWD/${input_file} \
    --annot ${annot} \
    --condition_varname ${condition_varname} \
    --min_patients ${min_patients}
    """
}
