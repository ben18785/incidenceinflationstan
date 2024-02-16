test_that("sample_true_cases_single_onset produces reasonable case draws", {
  observation_matrix <- dplyr::tibble(time_reported=c(1, 3, 5),
                                      cases_reported=c(1, 1, 1))
  reporting_parameters <- list(location=5, scale=3)
  max_cases <- 100
  s_params <- list(mean=10, sd=1)
  Rt <- 2
  t_max <- 30
  cases_history <- rep(4, t_max)
  day_onset <- 0
  w <- weights_series(t_max, s_params)
  mean_cases <- expected_cases(Rt, w, cases_history)
  case <- sample_true_cases_single_onset(observation_df=observation_matrix,
                                 cases_history=cases_history,
                                 max_cases=max_cases,
                                 Rt=Rt,
                                 day_onset=day_onset,
                                 serial_parameters=s_params,
                                 reporting_parameters=reporting_parameters)
  expect_true(case >= min(observation_matrix$cases_reported))
  expect_true(case <= max_cases)

  # when reporting mean is low, expect to have pretty much seen all cases
  reporting_parameters <- list(location=1, scale=1)
  cases <- sample_true_cases_single_onset(observation_df=observation_matrix,
                            cases_history=cases_history,
                            max_cases=max_cases,
                            Rt=Rt,
                            day_onset=day_onset,
                            serial_parameters=s_params,
                            reporting_parameters=reporting_parameters,
                            ndraws=20)
  expect_true(any(cases <= 5))

  # check that cases jump if serial interval looks further back
  cases_history <- c(rep(4, t_max / 2), rep(1000, t_max / 2))
  s_params <- list(mean=1, sd=1)
  cases <- sample_true_cases_single_onset(observation_df=observation_matrix,
                             cases_history=cases_history,
                             max_cases=max_cases,
                             Rt=Rt,
                             day_onset=day_onset,
                             serial_parameters=s_params,
                             reporting_parameters=reporting_parameters,
                             ndraws=20)
  s_params <- list(mean=25, sd=1)
  cases1 <- sample_true_cases_single_onset(observation_df=observation_matrix,
                             cases_history=cases_history,
                             max_cases=max_cases,
                             Rt=Rt,
                             day_onset=day_onset,
                             serial_parameters=s_params,
                             reporting_parameters=reporting_parameters,
                             ndraws=20)
  expect_true(mean(cases1) > mean(cases))

  max_cases <- 2
  expect_warning(sample_true_cases_single_onset(observation_df=observation_matrix,
                                                cases_history=cases_history,
                                                max_cases=max_cases,
                                                Rt=Rt,
                                                day_onset=day_onset,
                                                serial_parameters=s_params,
                                                reporting_parameters=reporting_parameters,
                                                ndraws=20))
  max_cases <- 0
  expect_error(sample_true_cases_single_onset(
    observation_df=observation_matrix,
    cases_history=cases_history,
    max_cases=max_cases,
    Rt=Rt,
    day_onset=day_onset,
    serial_parameters=s_params,
    reporting_parameters=reporting_parameters,
    ndraws=20))
})

test_that("sample_true_cases_single_onset produces consistent maximum", {
  observation_matrix <- dplyr::tibble(time_reported=c(1, 3, 5),
                                      cases_reported=c(1, 1, 1))
  reporting_parameters <- list(location=5, scale=3)
  max_cases <- 100
  s_params <- list(mean=10, sd=1)
  Rt <- 2
  t_max <- 30
  cases_history <- rep(4, t_max)
  day_onset <- 0
  max_observed_cases <- max(observation_matrix$cases_reported)
  possible_cases <- max_observed_cases:max_cases
  logps <- conditional_cases_logp(possible_cases, observation_matrix, cases_history,
                                  Rt, day_onset, s_params, reporting_parameters,
                                  is_gamma_delay = TRUE)

  map_val <- possible_cases[which.max(logps)]
  case <- sample_true_cases_single_onset(
    observation_df=observation_matrix,
    cases_history=cases_history,
    max_cases=max_cases,
    Rt=Rt,
    day_onset=day_onset,
    serial_parameters=s_params,
    reporting_parameters=reporting_parameters,
    maximise=T)
  expect_equal(case, map_val)
})

test_that("max_uncertain_days selects reasonable gamma distribution quanties", {
  r_params <- list(location=7, scale=0.01)
  expect_equal(round(max_uncertain_days(0.5, r_params, is_gamma_delay = TRUE)), 7)
  expect_true(max_uncertain_days(0.95, r_params, is_gamma_delay = TRUE) > 7)

  r_params <- list(location=7, scale=5)
  expect_true(max_uncertain_days(0.95, r_params, is_gamma_delay = TRUE) > 10)
  expect_equal(max_uncertain_days(0, r_params, is_gamma_delay = TRUE), 0)
})

test_that("sample_cases_history adds cases_estimated that look reasonable", {
  days_total <- 100
  kappa <- 1000
  r_params <- list(location=10, scale=3)
  v_Rt <- c(rep(1.5, 40), rep(0.4, 20), rep(1.5, 40))
  Rt_function <- stats::approxfun(1:days_total, v_Rt)
  s_params <- list(mean=5, sd=3)
  df <- generate_snapshots(days_total, Rt_function,
                           s_params, r_params,
                           kappa=kappa, thinned=T) %>%
    dplyr::mutate(reporting_piece_index=1)
  max_cases <- 5000
  r_params <- dplyr::tibble(location=10, scale=3, reporting_piece_index=1)
  df_est <- sample_cases_history(df, max_cases, Rt_function, s_params, r_params)
  expect_true(all.equal(df_est %>% dplyr::select(-cases_estimated),
                        df))
  expect_true(max(df_est$cases_estimated) < max_cases)
  expect_true(min(df_est$cases_estimated) > 0)
  expect_equal(sum(is.na(df_est$cases_estimated)), 0)

  df_group <- df_est %>%
    dplyr::group_by(time_onset) %>%
    dplyr::summarise(
      cases_reported=max(cases_reported),
      cases_estimated=min(cases_estimated),
      cases_true=mean(cases_true)) %>%
    dplyr::mutate(diff=cases_estimated-cases_reported) %>%
    dplyr::mutate(diff_true=cases_estimated-cases_true)
  expect_equal(sum(df_group$diff < 0), 0)
  expect_true(max(abs(df_group$diff_true)) < 300)

  # throws error when reporting_piece_index not in observation_onset_df
  df_tmp <- df %>%
    dplyr::select(-reporting_piece_index)
  expect_error(sample_cases_history(df_tmp, max_cases, Rt_function, s_params, r_params))
})

test_that("sample_cases_history yields a single case history when
maximising", {
  days_total <- 30
  kappa <- 1000
  r_params <- list(location=10, scale=3)
  v_Rt <- c(rep(1.5, 10), rep(0.4, 10), rep(1.5, 10))
  Rt_function <- stats::approxfun(1:days_total, v_Rt)
  s_params <- list(mean=5, sd=3)
  df <- generate_snapshots(days_total, Rt_function,
                           s_params, r_params,
                           kappa=kappa, thinned=T) %>%
    dplyr::mutate(reporting_piece_index=1)
  max_cases <- 5000
  r_params <- dplyr::tibble(
    reporting_piece_index=1,
    location=r_params$location,
    scale=r_params$scale
  )
  f_est <- function(i) {
    df_est <- sample_cases_history(df, max_cases, Rt_function,
                         s_params, r_params,
                         maximise = T)
    df_est$cases_estimated
  }
  cases <- purrr::map(seq(1, 4, 1), f_est)
  expect_true(all.equal(cases[[1]], cases[[2]]))
  expect_true(all.equal(cases[[2]], cases[[3]]))
  expect_true(all.equal(cases[[3]], cases[[4]]))
})

test_that("sample_or_maximise_from_gamma returns either draws or maximum of
          gamma", {
  shape <- 5
  rate <- 5
  ndraws <- 1
  val <- sample_or_maximise_gamma(shape, rate, ndraws)
  expect_equal(length(val), ndraws)
  ndraws <- 20
  vals <- sample_or_maximise_gamma(shape, rate, ndraws)
  expect_equal(length(vals), ndraws)

  # maximise
  val <- sample_or_maximise_gamma(shape, rate, ndraws,
                                  maximise=T)
  expect_equal(val, (shape - 1) / rate)
})

# tests for sampling Rt
days_total <- 100
Rt_1 <- 1.5
Rt_2 <- 1.0
Rt_3 <- 1.3
v_Rt <- c(rep(Rt_1, 40), rep(Rt_2, 20), rep(Rt_3, 40))
Rt_function <- stats::approxfun(1:days_total, v_Rt)
s_params <- list(mean=5, sd=3)
r_params <- list(location=10, scale=3)
kappa <- 1000
df <- generate_snapshots(days_total, Rt_function, s_params, r_params,
                         kappa=kappa) %>%
  dplyr::select(time_onset, cases_true)
Rt_indices <- unlist(purrr::map(seq(1, 5, 1), ~rep(., 20)))
Rt_index_lookup <- dplyr::tibble(
  time_onset=seq_along(Rt_indices),
  Rt_index=Rt_indices)
df <- df %>%
  dplyr::left_join(Rt_index_lookup, by = "time_onset") %>%
  dplyr::select(time_onset, cases_true, Rt_index) %>%
  unique()
Rt_prior <- list(shape=1, rate=1)

test_that("sample_Rt_single_piece returns reasonable values: these are
          basically functional tests", {

  ndraws <- 42
  Rt_vals <- sample_Rt_single_piece(
    2, df,
    Rt_prior, s_params,
    ndraws = ndraws,
    serial_max = 20)
  expect_equal(ndraws, length(Rt_vals))

  f_Rt <- function(piece) {
    sample_Rt_single_piece(
      piece, df,
      Rt_prior, s_params,
      ndraws = 1000,
      serial_max = 20)
  }

  Rt_vals <- sample_Rt_single_piece(
    2, df,
    Rt_prior, s_params,
    ndraws = 1000,
    serial_max = 40)

  # test for maximisation
  f_Rt <- function(i) sample_Rt_single_piece(
    2, df,
    Rt_prior, s_params,
    ndraws = 1000,
    serial_max = 20,
    maximise = T)
  expect_equal(length(f_Rt(1)), 1)
  Rt_vals1 <- purrr::map_dbl(seq(1, 3, 1), f_Rt)
  # no variation in maximum so should be no sd
  expect_equal(sd(Rt_vals1), 0)

  # mode of sampling distribution should be close to max
  mode <- function(x) {
    ux <- unique(x)
    ux[which.max(tabulate(match(x, ux)))]
  }
  Rt_vals <- sample_Rt_single_piece(
    2, df,
    Rt_prior, s_params,
    ndraws = 1000,
    serial_max = 20)
  expect_true(abs(mode(Rt_vals) - Rt_vals[1]) < 0.2)
})

test_that("sample_Rt returns sensible values", {
  ndraws <- 300
  Rt_df <- sample_Rt(df,
                     Rt_prior, s_params,
                     ndraws = ndraws,
                     serial_max = 20)

  # check shapes of outputs
  cnames <- colnames(Rt_df)
  expect_true(all.equal(cnames,
                        c("Rt_index", "draw_index", "Rt")))
  npieces <- max(df$Rt_index)
  expect_equal(npieces * ndraws, nrow(Rt_df))

  # check substance of outputs
  f_Rt <- function(piece) {
    val <- sample_Rt_single_piece(
      piece, df,
      Rt_prior, s_params,
      ndraws = 1000,
      serial_max = 20)
    mean(val)
  }
  indices <- seq(1, 5, 1)
  Rt_vals <- purrr::map_dbl(indices, f_Rt)
  single_df <- dplyr::tibble(Rt_index=indices,
                             Rt_single=Rt_vals)
  Rt_df <- Rt_df %>%
    dplyr::group_by(.data$Rt_index) %>%
    dplyr::summarise(Rt=mean(Rt)) %>%
    dplyr::left_join(single_df, by = "Rt_index") %>%
    dplyr::mutate(diff=Rt-Rt_single)
  expect_true(abs(sum(Rt_df$diff)) < 0.2)

  # check maximisation
  f_Rt <- function(piece) {
    val <- sample_Rt_single_piece(
      piece, df,
      Rt_prior, s_params,
      serial_max = 20,
      maximise = T)
    mean(val)
  }
  indices <- seq(1, 5, 1)
  Rt_vals <- purrr::map_dbl(indices, f_Rt)
  single_df <- dplyr::tibble(Rt_index=indices,
                             Rt_single=Rt_vals)
  Rt_df <- sample_Rt(df,
                     Rt_prior, s_params,
                     ndraws = ndraws,
                     serial_max = 20,
                     maximise = T)
  expect_equal(nrow(Rt_df), max(df$Rt_index))
  Rt_df <- Rt_df %>%
    dplyr::left_join(single_df, by = "Rt_index") %>%
    dplyr::mutate(diff=Rt-Rt_single)
  expect_true(abs(sum(Rt_df$diff)) < 0.0001)
})

days_total <- 30
df <- generate_snapshots(days_total, Rt_function,
                         s_params, r_params,
                         kappa=kappa, thinned=T) %>%
  dplyr::mutate(reporting_piece_index=1)
stan_model <- cmdstanr::cmdstan_model("../../src/stan/conditional_renewal.stan",
                                      compile_model_methods=TRUE,
                                      force_recompile=TRUE)


test_that("mcmc produces outputs of correct shape", {

  niter <- 2
  days_total <- 100
  r_params <- list(location=10, scale=3)
  s_params <- list(mean=5, sd=3)
  v_Rt <- c(rep(1.5, 40), rep(0.4, 20), rep(1.5, 40))
  Rt_function <- stats::approxfun(1:days_total, v_Rt)
  Rt_prior <- list(location=1, scale=1, is_gamma=TRUE)
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
    dplyr::left_join(Rt_index_lookup, by="time_onset") %>%
    dplyr::mutate(reporting_piece_index=1)

  initial_reporting_parameters <- dplyr::tibble(location=5, scale=3, reporting_piece_index=1)
  serial_parameters <- list(mean=5, sd=3)
  priors <- list(Rt=Rt_prior,
                 reporting=list(mean_mu=5,
                                mean_sigma=10,
                                sd_mu=3,
                                sd_sigma=5),
                 max_cases=5000)

  # test throws errors
  ## wrongly named cols
  wrong_df <- snapshot_with_Rt_index_df %>%
    dplyr::rename(time_onset_wrong=time_onset)
  expect_error(mcmc(niterations=niter,
                    wrong_df,
                    priors,
                    serial_parameters,
                    initial_cases_true,
                    initial_reporting_parameters,
                    initial_Rt,
                    stan_model),
               "Incorrect column names in snapshot_with_Rt_index_df")

  ## too many cols
  wrong_df <- snapshot_with_Rt_index_df %>%
    dplyr::mutate(unnecessary_col="hi")
  expect_error(mcmc(niterations=niter,
                    wrong_df,
                    priors,
                    serial_parameters,
                    initial_cases_true,
                    initial_reporting_parameters,
                    initial_Rt,
                    stan_model),
               "There should be either four or five columns in snapshot_with_Rt_index_df")

  # five columns but not reporting_piece_index
  wrong_df <- snapshot_with_Rt_index_df %>%
    dplyr::rename(reporting_piece_index_wrong=reporting_piece_index)
  expect_error(mcmc(niterations=niter,
                    wrong_df,
                    priors,
                    serial_parameters,
                    initial_cases_true,
                    initial_reporting_parameters,
                    initial_Rt,
                    stan_model),
               "The reporting delay indices must be provided through a column named 'reporting_piece_index'")

  # overdispersion parameter <= 0
  expect_error(mcmc(niter,
                    snapshot_with_Rt_index_df,
                    priors,
                    serial_parameters,
                    initial_cases_true,
                    initial_reporting_parameters,
                    initial_Rt,
                    stan_model,
                    is_negative_binomial = TRUE,
                    initial_overdispersion = -1))

  # test MCMC sampling
  niter <- 5
  res <- mcmc(niter,
              snapshot_with_Rt_index_df,
              priors,
              serial_parameters,
              initial_cases_true,
              initial_reporting_parameters,
              initial_Rt,
              stan_model)

  ## check outputs
  ### overall
  expect_equal(length(res), 3)

  ### cases
  cases_df <- res$cases
  expect_true(all.equal(c("time_onset", "cases_true", "iteration", "chain"),
                        colnames(cases_df)))
  expect_equal(min(cases_df$iteration), 1)
  expect_equal(max(cases_df$iteration), niter)
  expect_equal(min(cases_df$time_onset), min(snapshot_with_Rt_index_df$time_onset))
  expect_equal(max(cases_df$time_onset), max(snapshot_with_Rt_index_df$time_onset))
  expect_equal(min(cases_df$chain), 1)
  expect_equal(max(cases_df$chain), 1)

  ## Rt
  rt_df <- res$Rt
  expect_true(all.equal(c("iteration", "Rt_index", "Rt", "chain"),
                        colnames(rt_df)))
  expect_equal(min(rt_df$iteration), 1)
  expect_equal(max(rt_df$iteration), niter)
  expect_equal(min(rt_df$Rt_index), min(initial_Rt$Rt_index))
  expect_equal(max(rt_df$Rt_index), max(initial_Rt$Rt_index))
  expect_equal(min(rt_df$chain), 1)
  expect_equal(max(rt_df$chain), 1)

  # reporting delays
  reporting_df <- res$reporting
  expect_true(all.equal(c("reporting_piece_index", "location", "scale", "iteration", "chain"),
                        colnames(reporting_df)))
  expect_equal(min(reporting_df$iteration), 1)
  expect_equal(max(reporting_df$iteration), niter)
  expect_equal(min(reporting_df$chain), 1)
  expect_equal(max(reporting_df$chain), 1)

  # test optimisation
  niter <- 2 # needed since maximisation is iterative
  res <- mcmc(niterations=niter,
              snapshot_with_Rt_index_df,
              priors,
              serial_parameters,
              initial_cases_true,
              initial_reporting_parameters,
              initial_Rt,
              stan_model,
              maximise=TRUE)
  ### overall
  expect_equal(length(res), 3)

  ### cases
  cases_df <- res$cases
  expect_true(all.equal(c("time_onset", "cases_true", "iteration", "chain"),
                        colnames(cases_df)))
  expect_equal(min(cases_df$iteration), 1)
  expect_equal(max(cases_df$iteration), niter)
  expect_equal(min(cases_df$time_onset), min(snapshot_with_Rt_index_df$time_onset))
  expect_equal(max(cases_df$time_onset), max(snapshot_with_Rt_index_df$time_onset))
  expect_equal(min(cases_df$chain), 1)
  expect_equal(max(cases_df$chain), 1)

  ## Rt
  rt_df <- res$Rt
  expect_true(all.equal(c("iteration", "Rt_index", "Rt", "chain"),
                        colnames(rt_df)))
  expect_equal(min(rt_df$iteration), 1)
  expect_equal(max(rt_df$iteration), niter)
  expect_equal(min(rt_df$Rt_index), min(initial_Rt$Rt_index))
  expect_equal(max(rt_df$Rt_index), max(initial_Rt$Rt_index))
  expect_equal(min(rt_df$chain), 1)
  expect_equal(max(rt_df$chain), 1)

  # reporting delays
  reporting_df <- res$reporting
  expect_true(all.equal(c("reporting_piece_index", "location", "scale", "iteration", "chain"),
                        colnames(reporting_df)))
  expect_equal(min(reporting_df$iteration), 1)
  expect_equal(max(reporting_df$iteration), niter)
  expect_equal(min(reporting_df$chain), 1)
  expect_equal(max(reporting_df$chain), 1)


  # test MCMC sampling with NB model
  niter <- 5


  # sampling works if given overdispersion prior
  priors$overdispersion <- list(location=0, scale=2)
  res <- mcmc(niter,
              snapshot_with_Rt_index_df,
              priors,
              serial_parameters,
              initial_cases_true,
              initial_reporting_parameters,
              initial_Rt,
              stan_model,
              is_negative_binomial=TRUE,
              initial_overdispersion = 5)
  expect_equal(length(res), 4)
  expect_true("other" %in% names(res))
  other <- res$other
  expect_equal(nrow(other), 5)
  expect_true("overdispersion" %in% names(other))
  expect_true("iteration" %in% names(other))
})

test_that("multiple chains works", {

  days_total <- 100
  r_params <- list(location=10, scale=3)
  s_params <- list(mean=5, sd=3)
  v_Rt <- c(rep(1.5, 40), rep(0.4, 20), rep(1.5, 40))
  Rt_function <- stats::approxfun(1:days_total, v_Rt)
  Rt_prior <- list(location=1, scale=1, is_gamma=TRUE)
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

  initial_reporting_parameters <- list(location=5, scale=3)
  serial_parameters <- list(mean=5, sd=3)
  priors <- list(Rt=Rt_prior,
                 reporting=list(mean_mu=5,
                                mean_sigma=10,
                                sd_mu=3,
                                sd_sigma=5),
                 max_cases=5000)

  # multiple chains in serial
  niter <- 2
  nchains <- 2
  res <- mcmc(niter,
              snapshot_with_Rt_index_df,
              priors,
              serial_parameters,
              initial_cases_true,
              initial_reporting_parameters,
              initial_Rt,
              stan_model,
              nchains=nchains)
  cases_df <- res$cases
  Rt_df <- res$Rt
  rep_df <- res$reporting

  expect_equal(max(cases_df$chain), nchains)
  expect_equal(max(Rt_df$chain), nchains)
  expect_equal(max(rep_df$chain), nchains)

  nchains <- 3
  res <- mcmc(niter,
              snapshot_with_Rt_index_df,
              priors,
              serial_parameters,
              initial_cases_true,
              initial_reporting_parameters,
              initial_Rt,
              stan_model,
              nchains=nchains)
  cases_df <- res$cases
  Rt_df <- res$Rt
  rep_df <- res$reporting

  expect_equal(max(cases_df$chain), nchains)
  expect_equal(max(Rt_df$chain), nchains)
  expect_equal(max(rep_df$chain), nchains)


  # negative binomial model
  res <- mcmc(niter,
              snapshot_with_Rt_index_df,
              priors,
              serial_parameters,
              initial_cases_true,
              initial_reporting_parameters,
              initial_Rt,
              stan_model,
              is_negative_binomial = TRUE,
              initial_overdispersion = 10,
              nchains=nchains)
  cases_df <- res$cases
  Rt_df <- res$Rt
  rep_df <- res$reporting
  other_df <- res$other

  expect_equal(max(cases_df$chain), nchains)
  expect_equal(max(Rt_df$chain), nchains)
  expect_equal(max(rep_df$chain), nchains)
  expect_equal(max(other_df$chain), nchains)

  # negative binomial model and RW prior on Rt
  priors$Rt$is_gamma <- FALSE
  res <- mcmc(niter,
              snapshot_with_Rt_index_df,
              priors,
              serial_parameters,
              initial_cases_true,
              initial_reporting_parameters,
              initial_Rt,
              stan_model,
              is_negative_binomial = TRUE,
              initial_overdispersion = 10,
              initial_sigma = 10,
              nchains=nchains)
  cases_df <- res$cases
  Rt_df <- res$Rt
  rep_df <- res$reporting
  other_df <- res$other

  expect_equal(max(cases_df$chain), nchains)
  expect_equal(max(Rt_df$chain), nchains)
  expect_equal(max(rep_df$chain), nchains)
  expect_equal(max(other_df$chain), nchains)

  # multiple chains in parallel: TODO but waiting on cmdstanr changes
  # library(doParallel)
  # cl <- makeCluster(2)
  # registerDoParallel(cl)
  # res <- mcmc(niterations=niter,
  #             snapshot_with_Rt_index_df,
  #             priors,
  #             serial_parameters,
  #             initial_cases_true,
  #             initial_reporting_parameters,
  #             initial_Rt,
  #             reporting_metropolis_parameters=list(mean_step=0.25, sd_step=0.1),
  #             serial_max=40, p_gamma_cutoff=0.99, maximise=FALSE,
  #             nchains=nchains, is_parallel=TRUE, is_negative_binomial=FALSE)
  # stopCluster(cl)
  #
  # expect_equal(max(cases_df$chain), nchains)
  # expect_equal(max(Rt_df$chain), nchains)
  # expect_equal(max(rep_df$chain), nchains)

})

test_that("prepare_stan_data_and_init works", {

  days_total <- 100
  r_params <- list(location=10, scale=3)
  s_params <- list(mean=5, sd=3)
  v_Rt <- c(rep(1.5, 40), rep(0.4, 20), rep(1.5, 40))
  Rt_function <- stats::approxfun(1:days_total, v_Rt)
  Rt_prior <- list(location=1, scale=5, is_gamma=TRUE)
  kappa <- 10
  df <- generate_snapshots(days_total, Rt_function, s_params, r_params,
                           kappa=kappa)
  Rt_indices <- unlist(purrr::map(seq(1, 5, 1), ~rep(., 20)))
  Rt_index_lookup <- dplyr::tibble(
    time_onset=seq_along(Rt_indices),
    Rt_index=Rt_indices)
  df <- df %>%
    dplyr::left_join(Rt_index_lookup, by = "time_onset") %>%
    mutate(reporting_piece_index=1)

  current_values <- list(
    R=c(1.5, 1.5, 0.5, 1.5, 1.5),
    theta=matrix(c(10, 3), ncol = 2)
  )
  priors <- list(Rt=Rt_prior,
                 reporting=list(mean_mu=5,
                                mean_sigma=10,
                                sd_mu=3,
                                sd_sigma=5),
                 max_cases=5000)

  # without kappa or sigma
  res <- prepare_stan_data_and_init(df,
                      current_values,
                      priors,
                      is_negative_binomial=FALSE,
                      is_rw_prior=FALSE,
                      serial_parameters=s_params,
                      serial_max=40,
                      is_gamma_delay=TRUE)
  data_stan <- res$data
  expect_equal(data_stan$N, days_total)
  expect_equal(length(data_stan$w), 40)
  expect_equal(length(data_stan), 24)

  init_values <- res$init()
  expect_equal(init_values$R, current_values$R)
  expect_equal(init_values$theta, current_values$theta)

  # with kappa and/or sigma
  current_values$kappa <- 1
  res <- prepare_stan_data_and_init(df,
                                    current_values,
                                    priors,
                                    is_negative_binomial=TRUE,
                                    is_rw_prior=FALSE,
                                    serial_parameters=s_params,
                                    serial_max=40,
                                    is_gamma_delay=TRUE)
  data_stan <- res$data
  expect_equal(data_stan$N, days_total)
  expect_equal(length(data_stan$w), 40)
  expect_equal(length(data_stan), 24)

  init_values <- res$init()
  expect_equal(init_values$R, current_values$R)
  expect_equal(init_values$theta, current_values$theta)
  expect_equal(init_values$kappa[1], current_values$kappa)


  current_values$sigma <- 2
  res <- prepare_stan_data_and_init(df,
                                    current_values,
                                    priors,
                                    is_negative_binomial=TRUE,
                                    is_rw_prior=TRUE,
                                    serial_parameters=s_params,
                                    serial_max=40,
                                    is_gamma_delay=TRUE)
  data_stan <- res$data
  expect_equal(data_stan$N, days_total)
  expect_equal(length(data_stan$w), 40)
  expect_equal(length(data_stan), 24)

  init_values <- res$init()
  expect_equal(init_values$R, current_values$R)
  expect_equal(init_values$theta, current_values$theta)
  expect_equal(init_values$kappa[1], current_values$kappa)
  expect_equal(init_values$sigma[1], current_values$sigma)
})

test_that("check that get_step_size works ok", {

  days_total <- 100
  r_params <- list(location=10, scale=3)
  s_params <- list(mean=5, sd=3)
  v_Rt <- c(rep(1.5, 40), rep(0.4, 20), rep(1.5, 40))
  Rt_function <- stats::approxfun(1:days_total, v_Rt)
  Rt_prior <- list(location=1, scale=5, is_gamma=TRUE)
  kappa <- 10
  df <- generate_snapshots(days_total, Rt_function, s_params, r_params,
                           kappa=kappa)
  Rt_indices <- unlist(purrr::map(seq(1, 5, 1), ~rep(., 20)))
  Rt_index_lookup <- dplyr::tibble(
    time_onset=seq_along(Rt_indices),
    Rt_index=Rt_indices)
  df <- df %>%
    dplyr::left_join(Rt_index_lookup, by = "time_onset") %>%
    mutate(reporting_piece_index=1)

  current_values <- list(
    R=c(1.5, 1.5, 0.5, 1.5, 1.5),
    theta=matrix(c(10, 3), ncol = 2)
  )
  priors <- list(Rt=Rt_prior,
                 reporting=list(mean_mu=5,
                                mean_sigma=10,
                                sd_mu=3,
                                sd_sigma=5),
                 max_cases=5000)

  step_size <- get_step_size(df=df,
                             current_values=current_values,
                             priors=priors,
                             is_negative_binomial=FALSE,
                             is_rw_prior=FALSE,
                             serial_parameters=s_params,
                             serial_max=40,
                             is_gamma_delay=TRUE,
                             stan_model=stan_model,
                             n_warmup = 2)
  # really I am just testing that this runs
  expect_true(step_size > 0)
})

test_that("stan_metropolis runs and extracting works" , {

  days_total <- 100
  r_params <- list(location=10, scale=3)
  s_params <- list(mean=5, sd=3)
  v_Rt <- c(rep(1.5, 40), rep(0.4, 20), rep(1.5, 40))
  Rt_function <- stats::approxfun(1:days_total, v_Rt)
  Rt_prior <- list(location=1, scale=5, is_gamma=TRUE)
  kappa <- 10
  df <- generate_snapshots(days_total, Rt_function, s_params, r_params,
                           kappa=kappa)
  Rt_indices <- unlist(purrr::map(seq(1, 5, 1), ~rep(., 20)))
  Rt_index_lookup <- dplyr::tibble(
    time_onset=seq_along(Rt_indices),
    Rt_index=Rt_indices)
  df <- df %>%
    dplyr::left_join(Rt_index_lookup, by = "time_onset") %>%
    mutate(reporting_piece_index=1)

  current_values <- list(
    R=c(1.5, 1.5, 0.5, 1.5, 1.5),
    theta=matrix(c(10, 3), ncol = 2)
  )
  priors <- list(Rt=Rt_prior,
                 reporting=list(mean_mu=5,
                                mean_sigma=10,
                                sd_mu=3,
                                sd_sigma=5),
                 max_cases=5000)

  res <- stan_metropolis(0.25, df=df,
                             current_values=current_values,
                             priors=priors,
                             is_negative_binomial=FALSE,
                             is_rw_prior=FALSE,
                             serial_parameters=s_params,
                             serial_max=40,
                             is_gamma_delay=TRUE,
                             stan_model=stan_model,
                             n_iterations = 2)
  expect_equal(length(res$current$R), 5)
  expect_equal(length(res$current$theta), 2)

  # with kappa
  current_values$kappa <- 1
  res <- stan_metropolis(0.25, df=df,
                         current_values=current_values,
                         priors=priors,
                         is_negative_binomial=TRUE,
                         is_rw_prior=FALSE,
                         serial_parameters=s_params,
                         serial_max=40,
                         is_gamma_delay=TRUE,
                         stan_model=stan_model,
                         n_iterations = 3)
  expect_equal(length(res$current$R), 5)
  expect_equal(length(res$current$theta), 2)
  expect_equal(length(res$current$kappa), 1)

  # with sigma
  current_values$sigma <- 2
  priors$Rt$is_gamma <- FALSE
  res <- stan_metropolis(0.25, df=df,
                         current_values=current_values,
                         priors=priors,
                         is_negative_binomial=TRUE,
                         is_rw_prior=TRUE,
                         serial_parameters=s_params,
                         serial_max=40,
                         is_gamma_delay=TRUE,
                         stan_model=stan_model,
                         n_iterations = 3)
  expect_equal(length(res$current$R), 5)
  expect_equal(length(res$current$theta), 2)
  expect_equal(length(res$current$kappa), 1)
  expect_equal(length(res$current$sigma), 1)
})


