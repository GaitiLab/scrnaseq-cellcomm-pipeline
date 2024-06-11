# Unload all previously loaded packages + remove previous environment
rm(list = ls(all = TRUE))
pacman::p_unload()

require(GaitiLabUtils)
# Set working directory
set_wd()

# Load libraries
pacman::p_load(glue, data.table, tidyverse, stringr)
if (!interactive()) {
    # Define input arguments when running from bash
    parser <- setup_default_argparser(
        description = "Post-filtering/formatting",
        default_output = "output/402_aggregation_and_filtering"
    )
    parser$add_argument("--input_file",
        type = "character", default = NULL,
        help = "Path to '401_samples_interactions_mvoted.rds' file"
    )
    parser$add_argument("--min_patients", type = "integer", default = 2, help = "Minimum number of patients for an interaction to be kept")
    parser$add_argument("--condition_var", type = "character", help = "Name of condition variable", default = "Condition_dummy")
    args <- parser$parse_args()
} else {
    # Provide arguments here for local runs
    args <- list()
    args$log_level <- 5
    args$output_dir <- "output/test_individual_scripts/402_aggregation_and_filtering"
    args$input_file <- "output/test_individual_scripts/401_combine_samples/401_samples_interactions_mvoted.rds"
    args$condition_var <- "Condition"
    args$min_patients <- 1
}

# Set up logging
logr <- init_logging(log_level = args$log_level)
log_info(ifelse(interactive(),
    "Running interactively...",
    "Running from command line/terminal..."
))

log_info("Create output directory...")
create_dir(args$output_dir)

# Load additional libraries
scrnaseq.cellcomm::filter_by_detection_in_multi_samples(
    input_file = args$input_file,
    min_patients = args$min_patients,
    condition_var = args$condition_var,
    output_dir = args$output_dir
)

# TODO might need to change this to too, to not get a skewed view of the data
# glob_min_samples <- 1

# log_info("Loading input file...")
# input_file <- readRDS(args$input_file)

# # # r$> head(input_file)
# # # A tibble: 6 × 11
# #   Sample         source_target        complex_interaction n_methods in_liana in_cellchat in_cell2cell in_cpdb lenient_voting stringent_voting Region_Grouped
# #   <chr>          <chr>                <chr>                   <int>    <dbl>       <dbl>        <dbl>   <dbl> <lgl>          <lgl>            <chr>
# # 1 6234_2895153_A Microglia__Microglia A2M__LRP1                   4        1           1            1       1 TRUE           TRUE             TE
# # 2 6234_2895153_A Microglia__Microglia ACTR2__ADRB2                1        0           0            1       0 FALSE          FALSE            TE
# # 3 6234_2895153_A Microglia__Microglia ADAM10__AXL                 4        1           1            1       1 TRUE           TRUE             TE
# # 4 6234_2895153_A Microglia__Microglia ADAM10__CADM1               1        0           0            1       0 FALSE          FALSE            TE
# # 5 6234_2895153_A Microglia__Microglia ADAM10__CD44                2        0           0            1       1 FALSE          FALSE            TE
# # 6 6234_2895153_A Microglia__Microglia ADAM10__GPNMB               4        1           1            1       1 TRUE           TRUE             TE

# cols_oi <- c("Patient", "Sample", args$condition_var, "source_target", "complex_interaction")

# if (sum(str_detect(colnames(input_file), "Patient")) == 0) {
#     log_info("Patient column missing, assuming Patient = Sample...")
#     input_file <- input_file %>% mutate(Patient = Sample)
# }

# # ---- Lenient voting ---- #
# # Check detection of the same interaction in the same patient across samples for the same region
# lenient_voting <- input_file %>%
#     filter(lenient_voting) %>%
#     select(all_of(cols_oi)) %>%
#     group_by_at(vars((c(args$condition_var, "Patient", "source_target", "complex_interaction")))) %>%
#     reframe(lenient_detected_same_patient = paste0(Sample, collapse = ", "), lenient_N_samples_same_patient = n())
# head(lenient_voting %>% arrange(desc(lenient_N_samples_same_patient)))
# # # A tibble: 6 × 6
# #   Region_Grouped Patient source_target        complex_interaction lenient_detected_same_patient      lenient_N_samples_same_patient
# #   <chr>          <chr>   <chr>                <chr>               <chr>                                                       <int>
# # 1 PT             6419    Astrocyte__Astrocyte ALDH1A1__RORB       6419_cortex, 6419_enhancing_border                              2
# # 2 PT             6419    Astrocyte__Astrocyte ANGPTL4__SDC2       6419_cortex, 6419_enhancing_border                              2
# # 3 PT             6419    Astrocyte__Astrocyte ANGPTL4__SDC4       6419_cortex, 6419_enhancing_border                              2
# # 4 PT             6419    Astrocyte__Astrocyte CADM1__CADM1        6419_cortex, 6419_enhancing_border                              2
# # 5 PT             6419    Astrocyte__Astrocyte CNTN1__NOTCH2       6419_cortex, 6419_enhancing_border                              2
# # 6 PT             6419    Astrocyte__Astrocyte CNTN1__PTPRZ1       6419_cortex, 6419_enhancing_border                              2


# lenient_voting_by_region <- lenient_voting %>%
#     group_by_at(vars((c(args$condition_var, "source_target", "complex_interaction")))) %>%
#     reframe(
#         lenient_condition_patients = paste0(Patient, collapse = ", "), lenient_condition_n_patients = n(),
#         lenient_condition_n_samples = sum(lenient_N_samples_same_patient),
#         lenient_condition_samples = paste0(lenient_detected_same_patient, collapse = ", ")
#     )
# # r$> head(lenient_voting_by_region)
# # # A tibble: 6 × 7
# #   Region_Grouped source_target        complex_interaction lenient_condition lenient_condition_n_patients lenient_condition_n_samples lenient_condition_samples
# #   <chr>          <chr>                <chr>               <chr>                                    <int>                       <int> <chr>
# # 1 PT             Astrocyte__Astrocyte ALDH1A1__RORB       6237, 6419, 6509                             3                           4 6237_2222190_A, 6419_cortex, 6419_enhancing_border, 6509_cortex
# # 2 PT             Astrocyte__Astrocyte ANGPTL4__SDC2       6419                                         1                           2 6419_cortex, 6419_enhancing_border
# # 3 PT             Astrocyte__Astrocyte ANGPTL4__SDC4       6419                                         1                           2 6419_cortex, 6419_enhancing_border
# # 4 PT             Astrocyte__Astrocyte ANOS1__SDC2         6237                                         1                           1 6237_2222190_A
# # 5 PT             Astrocyte__Astrocyte APOE__ABCA1         6419, 6509                                   2                           2 6419_cortex, 6509_cortex
# # 6 PT             Astrocyte__Astrocyte APOE__LRP1          6419, 6509                                   2                           2 6419_cortex, 6509_cortex
# n_before <- nrow(lenient_voting_by_region)

# log_info(glue("Only keep interactions that are found in at least {args$min_patients} patients..."))
# lenient_voting_by_region <- lenient_voting_by_region %>% filter(lenient_condition_n_patients >= args$min_patients)
# n_after <- nrow(lenient_voting_by_region)
# # r$> head(lenient_voting_by_region)
# # # A tibble: 6 × 7
# #   Region_Grouped source_target        complex_interaction lenient_condition_patients lenient_condition_n_patients lenient_condition_n_samples lenient_condition_samples
# #   <chr>          <chr>                <chr>               <chr>                                    <int>                       <int> <chr>
# # 1 PT             Astrocyte__Astrocyte ALDH1A1__RORB       6237, 6419, 6509                             3                           4 6237_2222190_A, 6419_cortex, 6419_enhancing_border, 6509_cortex
# # 2 PT             Astrocyte__Astrocyte ANGPTL4__SDC2       6419                                         1                           2 6419_cortex, 6419_enhancing_border
# # 3 PT             Astrocyte__Astrocyte ANGPTL4__SDC4       6419                                         1                           2 6419_cortex, 6419_enhancing_border
# # 4 PT             Astrocyte__Astrocyte ANOS1__SDC2         6237                                         1                           1 6237_2222190_A
# # 5 PT             Astrocyte__Astrocyte APOE__ABCA1         6419, 6509                                   2                           2 6419_cortex, 6509_cortex
# # 6 PT             Astrocyte__Astrocyte APOE__LRP1          6419, 6509                                   2                           2 6419_cortex, 6509_cortex

# log_info(glue("Before filtering: {n_before}"))
# log_info(glue("After filtering: {n_after}"))

# # ---- Stringent voting ---- #
# stringent_voting <- input_file %>%
#     filter(stringent_voting) %>%
#     select(all_of(cols_oi)) %>%
#     group_by_at(vars(
#         (c(args$condition_var, "Patient", "source_target", "complex_interaction"))
#     )) %>%
#     reframe(stringent_detected_same_patient = paste0(Sample, collapse = ", "), stringent_N_samples_same_patient = n())
# # # A tibble: 6 × 6
# #   Region_Grouped Patient source_target        complex_interaction stringent_detected_same_patient stringent_N_samples_same_patient
# #   <chr>          <chr>   <chr>                <chr>               <chr>                                                      <int>
# # 1 PT             6234    Microglia__Microglia A2M__LRP1           6234_2895153_B                                                 1
# # 2 PT             6234    Microglia__Microglia ADAM10__AXL         6234_2895153_B                                                 1
# # 3 PT             6234    Microglia__Microglia ADAM10__CD44        6234_2895153_B                                                 1
# # 4 PT             6234    Microglia__Microglia ADAM10__GPNMB       6234_2895153_B                                                 1
# # 5 PT             6234    Microglia__Microglia ADAM10__IL6R        6234_2895153_B                                                 1
# # 6 PT             6234    Microglia__Microglia ADAM10__NOTCH1      6234_2895153_B                                                 1

# stringent_voting_by_region <- stringent_voting %>%
#     group_by_at(vars((c(args$condition_var, "source_target", "complex_interaction")))) %>%
#     reframe(stringent_condition_patients = paste0(Patient, collapse = ", "), stringent_condition_n_patients = n(), stringent_condition_n_samples = sum(stringent_N_samples_same_patient), stringent_condition_samples = paste0(stringent_detected_same_patient, collapse = ", "))
# # r$> head(stringent_voting_by_region)
# # # A tibble: 6 × 7
# #   Region_Grouped source_target        complex_interaction stringent_condition_patients stringent_condition_n_patients stringent_condition_n_samples stringent_condition_samples
# #   <chr>          <chr>                <chr>               <chr>                                        <int>                         <int> <chr>
# # 1 PT             Astrocyte__Astrocyte ALDH1A1__RORB       6237, 6419, 6509                                 3                             4 6237_2222190_A, 6419_cortex, 6419_enhancing_border, 6509_cortex
# # 2 PT             Astrocyte__Astrocyte ANGPTL4__SDC2       6419                                             1                             1 6419_enhancing_border
# # 3 PT             Astrocyte__Astrocyte ANGPTL4__SDC4       6419                                             1                             1 6419_enhancing_border
# # 4 PT             Astrocyte__Astrocyte APOE__ABCA1         6419, 6509                                       2                             2 6419_cortex, 6509_cortex
# # 5 PT             Astrocyte__Astrocyte APOE__LRP1          6419, 6509                                       2                             2 6419_cortex, 6509_cortex
# # 6 PT             Astrocyte__Astrocyte APOE__LRP4          6419, 6509                                       2                             2 6419_cortex, 6509_cortex

# n_before <- nrow(stringent_voting_by_region)

# log_info(glue("Only keep interactions that are found in at least {args$min_patients} patients..."))
# stringent_voting_by_region <- stringent_voting_by_region %>% filter(stringent_condition_n_patients >= args$min_patients)
# # r$> head(stringent_voting_by_region)
# # # A tibble: 6 × 7
# #   Region_Grouped source_target        complex_interaction stringent_condition_patients stringent_condition_n_patients stringent_condition_n_samples stringent_samples
# #   <chr>          <chr>                <chr>               <chr>                                        <int>                         <int> <chr>
# # 1 PT             Astrocyte__Astrocyte ALDH1A1__RORB       6237, 6419, 6509                                 3                             4 6237_2222190_A, 6419_cortex, 6419_enhancing_border, 6509_cortex
# # 2 PT             Astrocyte__Astrocyte APOE__ABCA1         6419, 6509                                       2                             2 6419_cortex, 6509_cortex
# # 3 PT             Astrocyte__Astrocyte APOE__LRP1          6419, 6509                                       2                             2 6419_cortex, 6509_cortex
# # 4 PT             Astrocyte__Astrocyte APOE__LRP4          6419, 6509                                       2                             2 6419_cortex, 6509_cortex
# # 5 PT             Astrocyte__Astrocyte BMP7__BMPR1A        6237, 6419                                       2                             2 6237_2222190_A, 6419_cortex
# # 6 PT             Astrocyte__Astrocyte BMP7__BMPR1B        6237, 6419, 6509                                 3                             3 6237_2222190_A, 6419_cortex, 6509_cortex

# n_after <- nrow(stringent_voting_by_region)

# log_info(glue("Before filtering: {n_before}"))
# log_info(glue("After filtering: {n_after}"))

# # ---- Combine stringent and lenient voting ---- #
# log_info("Combine stringent and lenient voting...")
# # Add column for visualization scripts later on
# lenient_voting_by_region <- lenient_voting_by_region %>% mutate(lenient_condition = TRUE)
# stringent_voting_by_region <- stringent_voting_by_region %>% mutate(stringent_condition = TRUE)

# combined_voting <- merge(lenient_voting_by_region, stringent_voting_by_region, by = c(args$condition_var, "source_target", "complex_interaction"), all = TRUE)
# # r$> head(combined_voting)
# #   Region_Grouped        source_target complex_interaction lenient_condition_patients lenient_condition_n_patients lenient_condition_n_samples                                       lenient_condition_samples stringent_condition_patients stringent_condition_n_patients
# # 1             PT Astrocyte__Astrocyte       ALDH1A1__RORB  6237, 6419, 6509                            3                           4 6237_2222190_A, 6419_cortex, 6419_enhancing_border, 6509_cortex    6237, 6419, 6509                              3
# # 2             PT Astrocyte__Astrocyte         APOE__ABCA1        6419, 6509                            2                           2                                        6419_cortex, 6509_cortex          6419, 6509                              2
# # 3             PT Astrocyte__Astrocyte          APOE__LRP1        6419, 6509                            2                           2                                        6419_cortex, 6509_cortex          6419, 6509                              2
# # 4             PT Astrocyte__Astrocyte          APOE__LRP4        6419, 6509                            2                           2                                        6419_cortex, 6509_cortex          6419, 6509                              2
# # 5             PT Astrocyte__Astrocyte          APOE__SDC2  6237, 6419, 6509                            3                           3                        6237_2222190_A, 6419_cortex, 6509_cortex                <NA>                             NA
# # 6             PT Astrocyte__Astrocyte        BMP7__BMPR1A        6237, 6419                            2                           2                                     6237_2222190_A, 6419_cortex          6237, 6419                              2
# #   stringent_condition_n_samples                                     stringent_condition_samples
# # 1                             4 6237_2222190_A, 6419_cortex, 6419_enhancing_border, 6509_cortex
# # 2                             2                                        6419_cortex, 6509_cortex
# # 3                             2                                        6419_cortex, 6509_cortex
# # 4                             2                                        6419_cortex, 6509_cortex
# # 5                            NA                                                            <NA>
# # 6                             2                                     6237_2222190_A, 6419_cortex
# log_info("Save results...")
# saveRDS(combined_voting, glue("{args$output_dir}/402a_filtering_detect_in_multi_samples.rds"))
# log_info("COMPLETED!")
