#' @title Compute correlation
#' @description Compute correlation between two vectors in a dataframe
#' @param df dataframe
#' @param i vector of indices
#' @param v1 name of group1
#' @param v2 name of group2
#' @param method method to use for computing correlation
#' @export
#' @importFrom dplyr pull sym %>%
compute_corr <- function(df, i, v1, v2, method = "pearson") {
    # Select data for computation
    df_sub <- df[i, ]
    return(cor(df_sub %>% pull(!!sym(v1)), df_sub %>% pull(!!sym(v2)), method = method))
}

#' @title Run bootstrap correlation analysis
#' @param pair named vector of length 2 (names: Var1 and Var2) with columns to select for analysis
#' @param df dataframe
#' @param stat_func function to use for computation of statistic, i.e. correlation
#' @param n_iter number of iterations to use for bootstrapping
#' @param method method to use for computing correlation
#' @export
#' @importFrom boot boot
run_bootstrap_corr <- function(pair, df, stat_func = compute_corr, n_iter = 500, method = "pearson") {
    v1 <- pair[["Var1"]]
    v2 <- pair[["Var2"]]

    bootcorr <- boot(df, statistic = stat_func, R = n_iter, v1 = v1, v2 = v2, method = method)

    # https://stackoverflow.com/questions/19963512/where-does-the-bootstrap-standard-error-live-in-the-boot-class
    # Compute mean and SE
    corr_mean <- mean(bootcorr$t)
    corr_obs <- bootcorr$t0
    corr_se <- sd(bootcorr$t)

    return(data.frame(
        pair = paste0(v1, "__", v2),
        state = v1,
        gene = v2,
        corr_mean = corr_mean,
        corr_obs = corr_obs,
        corr_se = corr_se
    ))
}
