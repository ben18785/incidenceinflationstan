#' Samples a case count arising on a given onset day
#'
#' The distribution being drawn from is given by:
#' \deqn{p(cases_true_t|data, Rt, reporting_params, serial_params) \propto p(data|cases_true, reporting_params)
#'  p(cases_true_t|cases_true_t_1, cases_true_t_2, ..., Rt, serial_params)}
#'
#' @param max_cases maximum possible cases thought to arise on a given day
#' @inheritParams conditional_cases_logp
#' @param ndraws number of draws of cases
#' @param maximise rather than sample a case count give the case count with the
#' maximum probability (by default is FALSE)
#'
#' @return a sampled case count arising on a given onset day
sample_true_cases_single_onset <- function(
  observation_df, cases_history, max_cases,
  Rt, day_onset, serial_parameters, reporting_parameters,
  ndraws=1, maximise=FALSE) {
  max_observed_cases <- max(observation_df$cases_reported)
  if(max_observed_cases > max_cases)
    stop("Max possible cases should be (much) greater than max observed cases.")
  possible_cases <- max_observed_cases:max_cases
  logps <- conditional_cases_logp(possible_cases, observation_df, cases_history,
                                  Rt, day_onset, serial_parameters, reporting_parameters)
  probs <- exp(logps - matrixStats::logSumExp(logps))
  if(dplyr::last(probs) > 0.01)
    warning(paste0("Cases too few for onset day: ", day_onset,
                   ". Increase max_cases."))
  if(maximise)
    possible_cases[which.max(probs)]
  else
    sample(possible_cases, ndraws, prob=probs, replace=TRUE)
}

#' Calculates max number of days we are uncertain about reporting
#'
#' @param p_gamma_cutoff a p value (0 <= p <= 1) indicating the threshold above which
#' we deem certainty
#' @inheritParams conditional_cases_logp
#'
#' @return a number of days
max_uncertain_days <- function(p_gamma_cutoff, reporting_parameters) {
  r_mean <- reporting_parameters$mean
  r_sd <- reporting_parameters$sd
  days_from_end <- qgamma_mean_sd(p_gamma_cutoff, r_mean, r_sd)
  days_from_end
}

#' Draws a possible history (or histories) of cases
#'
#' The distribution being drawn from at each time t is given by:
#' \deqn{p(cases_true_t|data, Rt, reporting_params, serial_params) \propto p(data|cases_true, reporting_params)
#'  p(cases_true_t|cases_true_t_1, cases_true_t_2, ..., Rt, serial_params)}
#'
#' @param observation_onset_df a tibble with three columns: time_onset, time_reported, cases_reported
#' @inheritParams sample_true_cases_single_onset
#' @inheritParams true_cases
#' @inheritParams max_uncertain_days
#'
#' @return a tibble with an extra cases_estimated column
#' @export
#' @importFrom rlang .data
sample_cases_history <- function(
  observation_onset_df, max_cases,
  Rt_function, serial_parameters, reporting_parameters,
  p_gamma_cutoff=0.99,
  maximise=FALSE) {

  uncertain_period <- max_uncertain_days(p_gamma_cutoff, reporting_parameters)
  start_uncertain_period <- max(observation_onset_df$time_onset) - uncertain_period
  observation_history_df <- observation_onset_df %>%
    dplyr::group_by(.data$time_onset) %>%
    dplyr::mutate(cases_estimated=ifelse(.data$time_onset < start_uncertain_period,
                                    max(.data$cases_reported), NA)) %>%
    dplyr::ungroup()
  onset_times <- unique(observation_history_df$time_onset)
  onset_times_uncertain_period <- onset_times[onset_times >= start_uncertain_period]

  for(i in seq_along(onset_times_uncertain_period)) {
    onset_time <- onset_times_uncertain_period[i]
    snapshots_at_onset_time_df <- observation_history_df %>%
      dplyr::filter(.data$time_onset==onset_time) %>%
      dplyr::select("time_reported", "cases_reported")
    pre_observation_df <- observation_history_df %>%
      dplyr::filter(.data$time_onset < onset_time) %>%
      dplyr::select("time_onset", "cases_estimated") %>%
      unique() %>%
      dplyr::arrange(dplyr::desc(.data$time_onset))
    cases_history <- pre_observation_df$cases_estimated
    Rt <- Rt_function(onset_time)
    case <- sample_true_cases_single_onset(
      observation_df=snapshots_at_onset_time_df,
      cases_history=cases_history,
      max_cases=max_cases,
      Rt=Rt,
      day_onset=onset_time,
      serial_parameters=serial_parameters,
      reporting_parameters=reporting_parameters,
      ndraws=1,
      maximise=maximise)
    index_onset_time <- which(observation_history_df$time_onset==onset_time)
    observation_history_df$cases_estimated[index_onset_time] <- case
  }
  observation_history_df
}

#' Draws from the gamma distribution or returns the value which maximises
#' it
#'
#' @param shape the shape parameter of a gamma distribution
#' @param rate the rate parameter of a gamma distribution
#' @param ndraws number of draws if maximise=FALSE
#' @param maximise whether to return the mode of the gamma distribution
#'
#' @return a value or (if ndraws > 1) a vector of values
sample_or_maximise_gamma <- function(shape, rate, ndraws, maximise=FALSE) {
  if(maximise)
    (shape - 1) / rate
  else
    stats::rgamma(ndraws, shape, rate)
}

#' Sample a single Rt value corresponding to a single piecewise-
#' constant element of an Rt vector
#'
#' @param Rt_piece_index the index of the Rt piece being sampled
#' @param cases_history_df a tibble with three columns: time_onset, cases_true
#' and Rt_index
#' @param Rt_prior_parameters a list with elements 'shape' and 'rate' describing
#' the gamma prior for Rt
#' @inheritParams sample_cases_history
#' @param ndraws number of draws of Rt
#' @inheritParams true_cases
#'
#' @return a draw (or draws) for Rt
#' @importFrom rlang .data
sample_Rt_single_piece <- function(
    Rt_piece_index, cases_history_df,
    Rt_prior_parameters, serial_parameters,
    serial_max=40, ndraws=1,
    maximise=FALSE) {
  short_df <- cases_history_df %>%
    dplyr::filter(.data$Rt_index <= Rt_piece_index)
  time_max_post_initial_period <- max(short_df$time_onset) - serial_max

  posterior_shape <- Rt_prior_parameters$shape
  posterior_rate <- Rt_prior_parameters$rate

  short_df <- short_df %>%
    dplyr::mutate(time_after_start = .data$time_onset - serial_max) %>%
    dplyr::mutate(is_observed_data=dplyr::if_else(
      .data$Rt_index == Rt_piece_index, 1, 0))
  onset_times <- short_df %>%
    dplyr::filter(.data$is_observed_data == 1) %>%
    dplyr::pull(.data$time_onset)

  w <- weights_series(serial_max, serial_parameters)
  for(i in seq_along(onset_times)) {
    onset_time <- onset_times[i]
    if(onset_time > 1) {
      true_cases <- short_df %>%
        dplyr::filter(.data$time_onset == onset_time) %>%
        dplyr::pull(.data$cases_true)
      posterior_shape <- posterior_shape + true_cases
      cases_history <- short_df %>%
        dplyr::filter(.data$time_onset < onset_time) %>%
        dplyr::arrange(dplyr::desc(.data$time_onset)) %>%
        dplyr::pull(.data$cases_true)
      diff_time <- serial_max - length(cases_history)
      if(diff_time <= 0)
        cases_history <- cases_history[1:serial_max]
      else
        cases_history <- c(cases_history, rep(0, diff_time))
      posterior_rate <- posterior_rate + sum(w * cases_history)
    }
  }
  sample_or_maximise_gamma(
    posterior_shape, posterior_rate, ndraws, maximise)
}


#' Sample piecewise-constant Rt values
#'
#' Models the renewal process as from a Poisson:
#' \deqn{cases_true_t ~ Poisson(Rt * \sum_tau=1^t_max w_t cases_true_t-tau))}
#' If an Rt value is given a gamma prior, this results in a posterior
#' distribution:
#' \deqn{Rt ~ gamma(alpha + cases_true_t, beta + \sum_tau=1^t_max w_t cases_true_t-tau))}
#' where alpha and beta are the shape and rate parameters of the gamma
#' prior distribution. Here, we assume that Rt is constant over a set of onset
#' times 'onset_time_set'. This means that the posterior for a single Rt value is given
#' by:
#' \deqn{Rt ~ gamma(alpha + \sum_{t in onset_time_set} cases_true_t,
#'       beta + \sum_{t in onset_time_set}\sum_tau=1^t_max w_t cases_true_t-tau))}
#' This function either returns a draw (or draws if ndraws>1) from
#' this posterior, or it returns the Rt set that maximises it
#' (if maximise=TRUE).
#'
#' @inheritParams sample_Rt_single_piece
#' @return a tibble with three columns: "Rt_piece_index", "draw_index", "Rt"
#' @export
sample_Rt <- function(cases_history_df,
                      Rt_prior_parameters,
                      serial_parameters,
                      serial_max=40,
                      ndraws=1,
                      maximise=FALSE) {
  Rt_piece_indices <- unique(cases_history_df$Rt_index)
  num_Rt_pieces <- length(Rt_piece_indices)
  if(maximise)
    ndraws <- 1
  draw_indices <- seq(1, ndraws, 1)
  m_draws <- matrix(nrow = num_Rt_pieces * ndraws,
                    ncol = 3)
  k <- 1
  for(i in seq_along(Rt_piece_indices)) {
    Rt_piece_index <- Rt_piece_indices[i]
    Rt_vals <- sample_Rt_single_piece(
      Rt_piece_index, cases_history_df,
      Rt_prior_parameters, serial_parameters,
      serial_max, ndraws, maximise=maximise)
    for(j in 1:ndraws) {
      m_draws[k, ] <- c(Rt_piece_index, j, Rt_vals[j])
      k <- k + 1
    }
  }
  colnames(m_draws) <- c("Rt_index", "draw_index", "Rt")
  m_draws <- m_draws %>%
    dplyr::as_tibble()
  m_draws
}

#' Propose new reporting parameters using normal kernel
#' centered at current values
#'
#' @param current_reporting_parameters named list of 'mean' and 'sd' of gamma distribution
#' characterising the reporting delay distribution
#' @param metropolis_parameters named list of 'mean_step', 'sd_step' containing
#' step sizes for Metropolis step
#'
#' @return list of reporting parameters
propose_reporting_parameters <- function(
  current_reporting_parameters,
  metropolis_parameters) {
  mean_now <- current_reporting_parameters$mean
  sd_now <- current_reporting_parameters$sd
  mean_stepsize <- metropolis_parameters$mean_step
  sd_stepsize <- metropolis_parameters$sd_step
  mean_proposed <- stats::rnorm(1, mean_now, mean_stepsize)
  sd_proposed <- stats::rnorm(1, sd_now, sd_stepsize)
  list(mean=mean_proposed, sd=sd_proposed)
}

#' Gamma prior for reporting parameters
#'
#' @inheritParams propose_reporting_parameters
#' @param prior_parameters named list with elements 'mean_mu', 'mean_sigma', 'sd_mu',
#' 'sd_sigma' representing the gamma prior parameters for the mean and sd
#' parameters of the reporting parameters (itself described by a gamma
#' distribution)
#'
#' @return a log-probability density
prior_reporting_parameters <- function(
  current_reporting_parameters,
  prior_parameters) {
  mean <- current_reporting_parameters$mean
  sd <- current_reporting_parameters$sd
  logp_mean <- dgamma_mean_sd(mean,
                              prior_parameters$mean_mu,
                              prior_parameters$mean_sigma,
                              log=TRUE)
  logp_sigma <- dgamma_mean_sd(sd,
                               prior_parameters$sd_mu,
                               prior_parameters$sd_sigma,
                               log=TRUE)
  logp_mean + logp_sigma
}


#' Sample reporting parameters using a single Metropolis step
#'
#' @inheritParams observation_process_all_times_logp
#' @inheritParams propose_reporting_parameters
#' @inheritParams prior_reporting_parameters
#'
#' @return list of reporting parameters
metropolis_step <- function(snapshot_with_true_cases_df,
                            current_reporting_parameters,
                            prior_parameters,
                            metropolis_parameters) {
  proposed_reporting_parameters <- propose_reporting_parameters(
    current_reporting_parameters,
    metropolis_parameters)
  logp_current <- observation_process_all_times_logp(
    snapshot_with_true_cases_df=snapshot_with_true_cases_df,
    reporting_parameters=current_reporting_parameters
  ) + prior_reporting_parameters(current_reporting_parameters,
                                 prior_parameters)
  logp_proposed <- observation_process_all_times_logp(
    snapshot_with_true_cases_df=snapshot_with_true_cases_df,
    reporting_parameters=proposed_reporting_parameters
  ) + prior_reporting_parameters(proposed_reporting_parameters,
                                 prior_parameters)

  log_r <- logp_proposed - logp_current
  log_u <- log(stats::runif(1))
  # nocov start
  if(log_r > log_u)
    proposed_reporting_parameters
  else
    current_reporting_parameters
  # nocov end
}

#' Sample reporting parameters using Metropolis MCMC
#'
#' @inheritParams metropolis_step
#' @param ndraws number of iterates of the Markov chain to simulate
#'
#' @return a tibble with three columns: "draw_index", "mean", "sd"
metropolis_steps <- function(
  snapshot_with_true_cases_df,
  current_reporting_parameters,
  prior_parameters,
  metropolis_parameters,
  ndraws) {

  m_reporting <- matrix(ncol = 3, nrow = ndraws)
  reporting_parameters <- current_reporting_parameters
  for(i in 1:ndraws) {
    reporting_parameters <- metropolis_step(
      snapshot_with_true_cases_df,
      reporting_parameters,
      prior_parameters,
      metropolis_parameters
    )
    m_reporting[i, ] <- c(i,
                          reporting_parameters$mean,
                          reporting_parameters$sd)
  }
  colnames(m_reporting) <- c("draw_index", "mean", "sd")
  m_reporting <- m_reporting %>%
    dplyr::as_tibble()
  m_reporting
}

#' Select reporting parameters by maximising log-probability
#'
#' @inheritParams metropolis_step
#'
#' @return a tibble with three columns: "draw_index", "mean, "sd"
maximise_reporting_logp <- function(
  snapshot_with_true_cases_df,
  current_reporting_parameters,
  prior_parameters) {

  objective_function <- function(theta) {
    -observation_process_all_times_logp(
      snapshot_with_true_cases_df,
      list(mean=theta[1], sd=theta[2])) +
      prior_reporting_parameters(
        list(mean=theta[1], sd=theta[2]),
        prior_parameters)
  }

  start_point <- c(current_reporting_parameters$mean,
                   current_reporting_parameters$sd)
  theta <- stats::optim(start_point, objective_function)$par
  reporting_parameters <- list(mean=theta[1],
                               sd=theta[2])
  dplyr::tibble(draw_index=1,
         mean=reporting_parameters$mean,
         sd=reporting_parameters$sd)
}

#' Draw reporting parameter values either by sampling or by
#' maximising
#'
#' @inheritParams metropolis_steps
#' @param maximise if true choose reporting parameters by maximising
#' log-probability; else (default) use Metropolis MCMC
#' to draw parameters
#'
#' @return a tibble with three columns: "draw_index", "mean, "sd"
#' @export
sample_reporting <- function(
  snapshot_with_true_cases_df,
  current_reporting_parameters,
  prior_parameters,
  metropolis_parameters,
  maximise=FALSE,
  ndraws=1) {
  if(maximise)
    reporting_parameters <- maximise_reporting_logp(
      snapshot_with_true_cases_df,
      current_reporting_parameters,
      prior_parameters)
  else
    reporting_parameters <- metropolis_steps(
      snapshot_with_true_cases_df,
      current_reporting_parameters,
      prior_parameters,
      metropolis_parameters,
      ndraws=ndraws)

  reporting_parameters
}

#' Runs MCMC or optimisation to estimate Rt, cases and reporting parameters
#'
#' @param niterations number of MCMC iterations to run or number of iterative maximisations to run
#' @param snapshot_with_Rt_index_df a tibble with
#' four columns: time_onset, time_reported, cases_reported, Rt_index
#' @param priors a named list with: 'Rt', 'reporting', 'max_cases'. These take
#' the form: 'Rt' is a named list with elements 'shape' and 'rate' describing
#' the gamma prior for each Rt; 'reporting' is a named list with elements
#' 'mean_mu', 'mean_sigma', 'sd_mu', 'sd_sigma' representing the gamma
#' prior parameters for the mean and sd parameters of the reporting parameters
#' (itself described by a gamma distribution); max cases controls the upper
#' limit of the discrete uniform distribution representing the prior on true
#' cases
#' @inheritParams sample_Rt_single_piece
#' @inheritParams max_uncertain_days
#' @param reporting_metropolis_parameters named list of 'mean_step', 'sd_step' containing
#' step sizes for Metropolis step
#' @param maximise whether to estimate MAP values of parameters (if true) or
#' sample parameter values using MCMC (if false). By default this is false.
#' @param initial_cases_true a tibble with two columns: "time_onset" and "cases_true", which represents initial
#' estimates of the true number of cases with each onset time.
#' @param initial_reporting_parameters a list with two named elements: 'mean', 'sd'
#' indicating an initial guess of the mean and sd of the reporting delay distribution
#' @param initial_Rt initial guess of the Rt values in each of the piecewise segments.
#' Provided in the form of a tibble with columns: 'Rt_index' and 'Rt'
#' @param print_to_screen prints progress of MCMC sampling to screen. Defaults to true.
#' @return a named list of three tibbles: "cases", "Rt" and "reporting" which contain estimates of the model parameters
#' @importFrom rlang .data
#' @export
mcmc_single <- function(
  niterations,
  snapshot_with_Rt_index_df,
  priors,
  serial_parameters,
  initial_cases_true,
  initial_reporting_parameters,
  initial_Rt,
  reporting_metropolis_parameters=list(mean_step=0.25, sd_step=0.1),
  serial_max=40, p_gamma_cutoff=0.99, maximise=FALSE, print_to_screen=TRUE) {

  cnames <- colnames(snapshot_with_Rt_index_df)
  expected_names <- c("time_onset", "time_reported",
                      "cases_reported", "Rt_index")

  if(sum(cnames %in% expected_names) != 4)
    stop("Incorrect column names in snapshot_with_Rt_index_df")

  if(length(cnames) != 4)
    stop("There should be only four columns in snapshot_with_Rt_index_df:
          time_onset, time_reported, cases_reported, Rt_index")

  df_running <- snapshot_with_Rt_index_df %>%
    dplyr::left_join(initial_cases_true, by = "time_onset") %>%
    dplyr::left_join(initial_Rt, by = "Rt_index")
  reporting_current <- initial_reporting_parameters

  reporting_samples <- matrix(ncol = 3,
                              nrow = niterations)
  num_Rts <- nrow(initial_Rt)
  Rt_samples <- matrix(nrow = num_Rts * niterations,
                       ncol = 3)
  num_cases_true <- unique(df_running$time_onset)

  if(print_to_screen) {
    progress_bar <- utils::txtProgressBar(
      min = 0,
      max = niterations,
      style = 3,
      width = niterations, # Needed to avoid multiple printings
      char = "=")

    init <- numeric(niterations)
    end <- numeric(niterations)
  }

  max_cases <- priors$max_cases

  k <- 1
  for(i in 1:niterations) {

    if(print_to_screen)
      init[i] <- Sys.time()

    # sample incidence
    Rt_current <- df_running %>%
      dplyr::select(time_onset, Rt) %>%
      unique()
    Rt_function <- stats::approxfun(Rt_current$time_onset,
                                    Rt_current$Rt)
    df_current <- df_running %>%
      dplyr::select(time_onset, time_reported, cases_reported)

    df_temp <- sample_cases_history(
      df_current, max_cases,
      Rt_function, serial_parameters, reporting_current,
      p_gamma_cutoff=p_gamma_cutoff,
      maximise=maximise)

    # sample Rt
    cases_history_df <- df_running %>%
      dplyr::select(time_onset, Rt_index) %>%
      dplyr::left_join(df_temp, by = "time_onset", relationship = "many-to-many") %>%
      dplyr::select(-c("cases_reported", "time_reported")) %>%
      dplyr::rename(cases_true=cases_estimated) %>%
      unique()

    df_Rt <- sample_Rt(cases_history_df,
                       priors$Rt,
                       serial_parameters,
                       serial_max,
                       ndraws=1,
                       maximise=maximise)
    # nocov start
    if(nrow(df_Rt) != num_Rts)
      stop("Number of Rts outputted not equal to initial Rt dims.")
    # nocov end
    # store Rts
    for(j in 1:num_Rts) {
      Rt_index <- df_Rt$Rt_index[j]
      Rt_temp <- df_Rt$Rt[j]
      Rt_samples[k, ] <- c(i, Rt_index, Rt_temp)
      k <- k + 1
    }

    # sample reporting parameters
    df_temp <- df_temp %>%
      dplyr::rename(cases_true=cases_estimated)
    reporting_temp <- sample_reporting(
      snapshot_with_true_cases_df=df_temp,
      current_reporting_parameters=reporting_current,
      prior_parameters=priors$reporting,
      metropolis_parameters=reporting_metropolis_parameters,
      maximise=maximise,
      ndraws=1)
    reporting_temp$draw_index <- NULL
    reporting_current <- reporting_temp
    reporting_samples[i, ] <- c(i,
                                reporting_current$mean,
                                reporting_current$sd)

    # update main df used in sampling
    df_running <- df_temp %>%
      dplyr::select(-cases_true) %>%
      dplyr::left_join(cases_history_df, by = "time_onset") %>%
      dplyr::left_join(df_Rt, by = "Rt_index") %>%
      dplyr::select(-draw_index)

    # store cases
    cases_history_df <- cases_history_df %>%
      dplyr::select(-Rt_index) %>%
      dplyr::mutate(iteration=i)
    if(i == 1) {
      cases_history_samples <- cases_history_df
    } else {
      cases_history_samples <- cases_history_samples %>%
        dplyr::bind_rows(cases_history_df)
    }

    if(print_to_screen) {
      end[i] <- Sys.time()
      utils::setTxtProgressBar(progress_bar, i)
      time <- round(lubridate::seconds_to_period(sum(end - init)), 0)

      # Estimated remaining time based on the
      # mean time that took to run the previous iterations
      est <- niterations * (mean(end[end != 0] - init[init != 0])) - time
      remaining <- round(lubridate::seconds_to_period(est), 0)

      cat(paste(" // Execution time:", tolower(as.character(time)),
                " // Estimated time remaining:", tolower(as.character(remaining))), "")
    }
  }

  Rt_samples <- Rt_samples %>%
    as.data.frame()
  colnames(Rt_samples) <- c("iteration", "Rt_index", "Rt")
  reporting_samples <- reporting_samples %>%
    as.data.frame()
  colnames(reporting_samples) <- c("iteration", "mean", "sd")

  list(cases=cases_history_samples,
       Rt=Rt_samples,
       reporting=reporting_samples)
}

#' Combines Markov chains across multiple runs of mcmc_single
#'
#' @param list_of_results a list of results, where each element is a result of running mcmc_single
#' @return a named list of three tibbles: "cases", "Rt" and "reporting" which
#' contain estimates of the model parameters with chain index  included
combine_chains <- function(list_of_results) {

  for(i in seq_along(list_of_results)) {
    res <- list_of_results[[i]]
    cases_df <- res$cases
    Rt_df <- res$Rt
    reporting_df <- res$reporting
    cases_df$chain <- i
    Rt_df$chain <- i
    reporting_df$chain <- i

    if(i == 1) {
      cases_overall <- cases_df
      Rt_overall <- Rt_df
      reporting_overall <- reporting_df
    } else {
      cases_overall <- cases_overall %>% dplyr::bind_rows(cases_df)
      Rt_overall <- Rt_overall %>% dplyr::bind_rows(Rt_df)
      reporting_overall <- reporting_overall %>% dplyr::bind_rows(reporting_df)
    }
  }
  list(
    cases=cases_overall,
    Rt=Rt_overall,
    reporting=reporting_overall
  )
}

#' Runs MCMC or optimisation to estimate Rt, cases and reporting parameters
#'
#' @param niterations number of MCMC iterations to run or number of iterative maximisations to run
#' @param data a tibble with five columns: time_onset, time_reported, cases_reported, Rt_index, reporting_index
#' @param priors a named list with: 'Rt', 'reporting', 'max_cases'. These take
#' the form: 'Rt' is a named list with elements 'shape' and 'rate' describing
#' the gamma prior for each Rt; 'reporting' is a named list with elements
#' 'mean_mu', 'mean_sigma', 'sd_mu', 'sd_sigma' representing the gamma
#' prior parameters for the mean and sd parameters of the reporting parameters
#' (itself described by a gamma distribution); max cases controls the upper
#' limit of the discrete uniform distribution representing the prior on true
#' cases
#' @inheritParams sample_Rt_single_piece
#' @inheritParams max_uncertain_days
#' @param reporting_metropolis_parameters named list of 'mean_step', 'sd_step' containing
#' step sizes for Metropolis step
#' @param maximise whether to estimate MAP values of parameters (if true) or
#' sample parameter values using MCMC (if false). By default this is false.
#' @param initial_cases_true a tibble with two columns: "time_onset" and "cases_true", which represents initial
#' estimates of the true number of cases with each onset time.
#' @param initial_reporting_parameters a list with two named elements: 'mean', 'sd'
#' indicating an initial guess of the mean and sd of the reporting delay distribution
#' @param initial_Rt initial guess of the Rt values in each of the piecewise segments.
#' Provided in the form of a tibble with columns: 'Rt_index' and 'Rt'
#' @param print_to_screen prints progress of MCMC sampling to screen. Defaults to true. Disabled when is_parallel is TRUE.
#' @param nchains number of Markov chains to run. Defaults to 1
#' @param is_parallel Boolean to indicate whether or not to run chains in parallel. Defaults to FALSE.
#' @return a named list of three tibbles: "cases", "Rt" and "reporting" which contain estimates of the model parameters
#' @export
#' @importFrom rlang .data
mcmc <- function(
    niterations,
    data,
    priors,
    serial_parameters,
    initial_cases_true,
    initial_reporting_parameters,
    initial_Rt,
    reporting_metropolis_parameters=list(mean_step=0.25, sd_step=0.1),
    serial_max=40, p_gamma_cutoff=0.99, maximise=FALSE, print_to_screen=TRUE,
    nchains=1, is_parallel=FALSE) {

  if(nchains==1) {
    res <- mcmc_single(niterations,
                       snapshot_with_Rt_index_df,
                       priors,
                       serial_parameters,
                       initial_cases_true,
                       initial_reporting_parameters,
                       initial_Rt,
                       reporting_metropolis_parameters,
                       serial_max,
                       p_gamma_cutoff,
                       maximise,
                       print_to_screen)
    res$cases$chain <- 1
    res$Rt$chain <- 1
    res$reporting$chain <- 1
  } else {

    list_of_results <- vector(mode = "list", length = nchains)
    if(is_parallel) {

      if (requireNamespace("foreach", quietly = TRUE)) {

        f_run_single <- function() {
          mcmc_single(niterations,
                      snapshot_with_Rt_index_df,
                      priors,
                      serial_parameters,
                      initial_cases_true,
                      initial_reporting_parameters,
                      initial_Rt,
                      reporting_metropolis_parameters,
                      serial_max,
                      p_gamma_cutoff,
                      maximise,
                      print_to_screen)
        }

        list_of_results <- foreach::foreach(i=1:nchains, .export = "mcmc_single") %dopar% {
          res <- f_run_single()
        }

      } else {
        # nocov start
        warning("The foreach package must be installed to use this functionality")
        #Either exit or do something without rgl
        return(NULL)
        # nocov end
      }

    } else {
      for(i in seq_along(1:nchains)) {
        res <- mcmc_single(niterations,
                           snapshot_with_Rt_index_df,
                           priors,
                           serial_parameters,
                           initial_cases_true,
                           initial_reporting_parameters,
                           initial_Rt,
                           reporting_metropolis_parameters,
                           serial_max,
                           p_gamma_cutoff,
                           maximise,
                           print_to_screen)
        list_of_results[[i]] <- res
      }
    }

    res <- combine_chains(list_of_results)
  }

  res
}
