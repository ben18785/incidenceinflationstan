test_that("check that convert_results_to_posterior_format converts an MCMC result
          into posterior format", {

    days_total <- 100
    r_params <- list(mean=10, sd=3)
    s_params <- list(mean=5, sd=3)

    v_Rt <- c(rep(1.5, 40), rep(0.4, 20), rep(1.5, 40))
    Rt_function <- stats::approxfun(1:days_total, v_Rt)
    kappa <- 10
    df <- generate_snapshots(days_total, Rt_function, s_params, r_params,
                             kappa=kappa)

    snapshot_with_Rt_index_df <- df
    initial_cases_true <- df %>%
      dplyr::select(time_onset, cases_true) %>%
      unique()
    snapshot_with_Rt_index_df <- snapshot_with_Rt_index_df %>%
      dplyr::select(-cases_true)
    initial_Rt <- tidyr::tribble(~Rt_index, ~Rt,
                                 1, 1.5,
                                 2, 1.5,
                                 3, 0.4,
                                 4, 1.5,
                                 5, 1.5)
    Rt_indices <- unlist(purrr::map(seq(1, 5, 1), ~rep(., 20)))

    Rt_index_lookup <- tidyr::tibble(
      time_onset=seq_along(Rt_indices),
      Rt_index=Rt_indices)
    snapshot_with_Rt_index_df <- snapshot_with_Rt_index_df %>%
      dplyr::left_join(Rt_index_lookup)

    initial_reporting_parameters <- list(mean=5, sd=3)
    serial_parameters <- list(mean=5, sd=3)
    Rt_prior <- list(shape=1, rate=1)
    priors <- list(Rt=Rt_prior,
                   reporting=list(mean_mu=5,
                                  mean_sigma=10,
                                  sd_mu=3,
                                  sd_sigma=5),
                   max_cases=5000)

    # multiple chains in serial
    niter <- 2
    nchains <- 2
    res <- mcmc(niterations=niter,
                snapshot_with_Rt_index_df,
                priors,
                serial_parameters,
                initial_cases_true,
                initial_reporting_parameters,
                initial_Rt,
                reporting_metropolis_parameters=list(mean_step=0.25, sd_step=0.1),
                serial_max=40, p_gamma_cutoff=0.99, maximise=FALSE,
                nchains=nchains)

    results <- convert_results_to_posterior_format(res)
    cases_df <- results$cases
    Rt_df <- results$Rt
    rep_df <- results$reporting

    case_names <- purrr::map_chr(seq(1, 100, 1), ~paste0("cases_true_", .))
    expect_true(all.equal(colnames(cases_df), c(".chain", ".iteration", ".draw", case_names)))
    expect_true(all.equal(colnames(Rt_df), c(".chain", ".iteration", ".draw", "Rt_1", "Rt_2", "Rt_3", "Rt_4", "Rt_5")))
    expect_true(all.equal(colnames(rep_df), c(".chain", ".iteration", ".draw", "reporting_piece_index", "mean", "sd")))
})
