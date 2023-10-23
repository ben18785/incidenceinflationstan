test_that("sample_true_cases_single_onset produces reasonable case draws", {
  observation_matrix <- dplyr::tibble(time_reported=c(1, 3, 5),
                                      cases_reported=c(1, 1, 1))
  reporting_parameters <- list(mean=5, sd=3)
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
  reporting_parameters <- list(mean=1, sd=1)
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
  reporting_parameters <- list(mean=5, sd=3)
  max_cases <- 100
  s_params <- list(mean=10, sd=1)
  Rt <- 2
  t_max <- 30
  cases_history <- rep(4, t_max)
  day_onset <- 0
  max_observed_cases <- max(observation_matrix$cases_reported)
  possible_cases <- max_observed_cases:max_cases
  logps <- conditional_cases_logp(possible_cases, observation_matrix, cases_history,
                                  Rt, day_onset, s_params, reporting_parameters)

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
  r_params <- list(mean=7, sd=0.01)
  expect_equal(round(max_uncertain_days(0.5, r_params)), 7)
  expect_true(max_uncertain_days(0.95, r_params) > 7)

  r_params <- list(mean=7, sd=5)
  expect_true(max_uncertain_days(0.95, r_params) > 10)
  expect_equal(max_uncertain_days(0, r_params), 0)
})

test_that("sample_cases_history adds cases_estimated that look reasonable", {
  days_total <- 100
  kappa <- 1000
  r_params <- list(mean=10, sd=3)
  v_Rt <- c(rep(1.5, 40), rep(0.4, 20), rep(1.5, 40))
  Rt_function <- stats::approxfun(1:days_total, v_Rt)
  s_params <- list(mean=5, sd=3)
  df <- generate_snapshots(days_total, Rt_function,
                           s_params, r_params,
                           kappa=kappa, thinned=T)
  max_cases <- 5000
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
})

test_that("sample_cases_history yields a single case history when
maximising", {
  days_total <- 30
  kappa <- 1000
  r_params <- list(mean=10, sd=3)
  v_Rt <- c(rep(1.5, 10), rep(0.4, 10), rep(1.5, 10))
  Rt_function <- stats::approxfun(1:days_total, v_Rt)
  s_params <- list(mean=5, sd=3)
  df <- generate_snapshots(days_total, Rt_function,
                           s_params, r_params,
                           kappa=kappa, thinned=T)
  max_cases <- 5000
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
r_params <- list(mean=10, sd=3)
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
  # too few data so returns prior draws
  Rt_vals <- f_Rt(1)
  prior_mean <- Rt_prior$shape / Rt_prior$rate
  expect_true(abs(mean(Rt_vals) - prior_mean) < 0.2)

  Rt_vals <- f_Rt(2)
  expect_true(abs(mean(Rt_vals) - Rt_1) < 0.2)
  Rt_vals <- f_Rt(3)
  expect_true(abs(mean(Rt_vals) - Rt_2) < 0.2)
  Rt_vals <- f_Rt(4)
  expect_true(abs(mean(Rt_vals) - Rt_3) < 0.2)
  Rt_vals <- f_Rt(5)
  expect_true(abs(mean(Rt_vals) - Rt_3) < 0.2)

  # increase serial_max so that more data points are
  # discarded: resulting in prior draws
  Rt_vals <- sample_Rt_single_piece(
    2, df,
    Rt_prior, s_params,
    ndraws = 1000,
    serial_max = 40)
  expect_true(abs(mean(Rt_vals) - prior_mean) < 0.2)

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

test_that("prior_reporting_parameters works ok", {
  current_reporting_parameters <- list(mean=3, sd=2)
  prior_params <- list(mean_mu=3, mean_sigma=2,
                       sd_mu=5, sd_sigma=1)
  val <- prior_reporting_parameters(current_reporting_parameters,
                                    prior_params)
  val1 <- dgamma_mean_sd(current_reporting_parameters$mean,
                         prior_params$mean_mu,
                         prior_params$mean_sigma,
                         log=TRUE)
  val2 <- dgamma_mean_sd(current_reporting_parameters$sd,
                         prior_params$sd_mu,
                         prior_params$sd_sigma,
                         log=TRUE)
  expect_equal(val, val1 + val2)
})

test_that("propose_reporting_parameters works ok", {
  r_params <- list(mean=4, sd=4)
  met_params <- list(mean_step=0.1, sd_step=0.2)
  r_params1 <- propose_reporting_parameters(
    r_params, met_params)
  expect_true(abs(r_params$mean - r_params1$mean) < 3)
  expect_true(abs(r_params$sd - r_params1$sd) < 3)

  met_params <- list(mean_step=0.000001, sd_step=0.2)
  r_params1 <- propose_reporting_parameters(
    r_params, met_params)
  expect_true(abs(r_params$mean - r_params1$mean) < 0.1)
})

days_total <- 30
df <- generate_snapshots(days_total, Rt_function,
                         s_params, r_params,
                         kappa=kappa, thinned=T)

test_that("metropolis_step works as expected", {
  met_params <- list(mean_step=0.01, sd_step=0.01)
  prior_params <- list(mean_mu=2, mean_sigma=100,
                       sd_mu=2, sd_sigma=100)
  r_params1 <- metropolis_step(df, r_params,
                               prior_params,
                               met_params)
  expect_true(abs(r_params1$mean - r_params$mean) < 0.2)
  expect_true(abs(r_params1$sd - r_params$sd) < 0.2)
})

test_that("metropolis_steps returns multiple steps", {
  met_params <- list(mean_step=0.01, sd_step=0.01)
  ndraws <- 10
  rep_prior_params <- list(mean_mu=5, sd_mu=3,
                           mean_sigma=5, sd_sigma=3)
  output <- metropolis_steps(
    snapshot_with_true_cases_df=df,
    current_reporting_parameters=r_params,
    prior_parameters=rep_prior_params,
    metropolis_parameters=met_params,
    ndraws=ndraws)
  expect_equal(nrow(output), ndraws)
  expect_true(all.equal(colnames(output),
                        c("draw_index", "mean", "sd")))
})

test_that("maximise_reporting_logp maximises prob", {
  prior_params <- list(mean_mu=2, mean_sigma=100,
                       sd_mu=2, sd_sigma=100)
  output <- maximise_reporting_logp(df, r_params, prior_params)
  expect_true(abs(output$mean - r_params$mean) < 0.4)
  expect_true(abs(output$sd - r_params$sd) < 0.4)
  expect_equal(nrow(output), 1)
})

test_that("sample_reporting produces output of correct shape", {
  prior_params <- list(mean_mu=2, mean_sigma=100,
                       sd_mu=2, sd_sigma=100)
  met_params <- list(mean_step=0.01, sd_step=0.01)
  output <- sample_reporting(df, r_params, prior_params, met_params)
  expect_equal(nrow(output), 1)
  expect_equal(max(output$draw_index), 1)

  ndraws <- 2
  output <- sample_reporting(df, r_params, prior_params, met_params,
                             ndraws=ndraws)
  expect_equal(nrow(output), ndraws)
  expect_equal(max(output$draw_index), ndraws)

  output <- sample_reporting(df, r_params, prior_params, met_params,
                             maximise=T)
  expect_equal(nrow(output), 1)
  expect_equal(max(output$draw_index), 1)
})


test_that("mcmc produces outputs of correct shape", {

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
                    reporting_metropolis_parameters=list(mean_step=0.25, sd_step=0.1),
                    serial_max=40, p_gamma_cutoff=0.99, maximise=FALSE))

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
                    reporting_metropolis_parameters=list(mean_step=0.25, sd_step=0.1),
                    serial_max=40, p_gamma_cutoff=0.99, maximise=FALSE))

  # test MCMC sampling
  niter <- 5
  res <- mcmc(niterations=niter,
              snapshot_with_Rt_index_df,
              priors,
              serial_parameters,
              initial_cases_true,
              initial_reporting_parameters,
              initial_Rt,
              reporting_metropolis_parameters=list(mean_step=0.25, sd_step=0.1),
              serial_max=40, p_gamma_cutoff=0.99, maximise=FALSE)

  ## check outputs
  ### overall
  expect_equal(length(res), 3)

  ### cases
  cases_df <- res$cases
  expect_true(all.equal(c("time_onset", "cases_true", "iteration"),
                        colnames(cases_df)))
  expect_equal(min(cases_df$iteration), 1)
  expect_equal(max(cases_df$iteration), niter)
  expect_equal(min(cases_df$time_onset), min(df$time_onset))
  expect_equal(max(cases_df$time_onset), max(df$time_onset))

  ## Rt
  rt_df <- res$Rt
  expect_true(all.equal(c("iteration", "Rt_index", "Rt"),
                        colnames(rt_df)))
  expect_equal(min(rt_df$iteration), 1)
  expect_equal(max(rt_df$iteration), niter)
  expect_equal(min(rt_df$Rt_index), min(initial_Rt$Rt_index))
  expect_equal(max(rt_df$Rt_index), max(initial_Rt$Rt_index))

  # reporting delays
  reporting_df <- res$reporting
  expect_true(all.equal(c("iteration", "mean", "sd"),
                        colnames(reporting_df)))
  expect_equal(min(reporting_df$iteration), 1)
  expect_equal(max(reporting_df$iteration), niter)

  # test optimisation
  niter <- 2 # needed since maximisation is iterative
  res <- mcmc(niterations=niter,
              snapshot_with_Rt_index_df,
              priors,
              serial_parameters,
              initial_cases_true,
              initial_reporting_parameters,
              initial_Rt,
              reporting_metropolis_parameters=list(mean_step=0.25, sd_step=0.1),
              serial_max=40, p_gamma_cutoff=0.99, maximise=TRUE)
  ### overall
  expect_equal(length(res), 3)

  ### cases
  cases_df <- res$cases
  expect_true(all.equal(c("time_onset", "cases_true", "iteration"),
                        colnames(cases_df)))
  expect_equal(min(cases_df$iteration), 1)
  expect_equal(max(cases_df$iteration), niter)
  expect_equal(min(cases_df$time_onset), min(df$time_onset))
  expect_equal(max(cases_df$time_onset), max(df$time_onset))

  ## Rt
  rt_df <- res$Rt
  expect_true(all.equal(c("iteration", "Rt_index", "Rt"),
                        colnames(rt_df)))
  expect_equal(min(rt_df$iteration), 1)
  expect_equal(max(rt_df$iteration), niter)
  expect_equal(min(rt_df$Rt_index), min(initial_Rt$Rt_index))
  expect_equal(max(rt_df$Rt_index), max(initial_Rt$Rt_index))

  # reporting delays
  reporting_df <- res$reporting
  expect_true(all.equal(c("iteration", "mean", "sd"),
                        colnames(reporting_df)))
  expect_equal(min(reporting_df$iteration), 1)
  expect_equal(max(reporting_df$iteration), niter)
})
