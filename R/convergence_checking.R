#' Converts a single tibble to posterior format
#'
#' @param res_df one of the three tibbles returned by running mcmc
#' @param parameter_cols the column(s) containing the parameter draws
#' @param names_cols another column specifying an index for parameter draws
#' @param names_prefix the prefix for naming the columns of the output tibble
#'
#' @return a tibble in the format required by posterior package
convert_df_to_posterior_format <- function(res_df, parameter_cols, names_cols,
                                           names_prefix) {
  res_df %>%
    dplyr::rename(
    .chain=chain,
    .iteration=iteration
  ) %>%
    dplyr::mutate(.draw=.iteration) %>%
    tidyr::pivot_wider(id_cols = c(.chain, .iteration, .draw),
                values_from = all_of(parameter_cols),
                names_from = all_of(names_cols),
                names_prefix = names_prefix)
}

#' Converts results to format required for posterior package
#'
#' @param results an object returned by running mcmc
#'
#' @return a named list of three tibbles: "cases", "Rt" and "reporting" which contain estimates of the model parameters
#' @export
convert_results_to_posterior_format <- function(results) {

  cases_df <- results$cases
  Rt_df <- results$Rt
  rep_df <- results$reporting

  cases_df <- convert_df_to_posterior_format(cases_df, "cases_true", "time_onset", "cases_true_")
  Rt_df <- convert_df_to_posterior_format(Rt_df, "Rt", "Rt_index", "Rt_")
  rep_df <- rep_df %>%
    dplyr::rename(
      .chain=chain,
      .iteration=iteration
    ) %>%
    dplyr::mutate(.draw=.iteration) %>%
    dplyr::relocate(c(.chain, .iteration, .draw)) %>%
    tidyr::pivot_wider(
      id_cols = c(.chain, .iteration, .draw),
      names_from=reporting_piece_index,
      values_from = c(mean, sd))

  is_negative_binomial <- FALSE
  if("overdispersion" %in% names(results)) {
    overdispersion_df <- results$overdispersion
    overdispersion_df <- overdispersion_df %>%
      dplyr::rename(
        .chain=chain,
        .iteration=iteration
      ) %>%
      dplyr::mutate(.draw=.iteration) %>%
      dplyr::relocate(c(.chain, .iteration, .draw))
    is_negative_binomial <- TRUE
  }

  list_results <- list(
    cases=cases_df,
    Rt=Rt_df,
    reporting=rep_df
  )

  if(is_negative_binomial)
    list_results$overdispersion <- overdispersion_df
  list_results
}
