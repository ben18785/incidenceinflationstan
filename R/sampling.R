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
  ndraws=1, maximise=FALSE, is_gamma_delay=TRUE, kappa=NULL, is_negative_binomial=FALSE) {

  max_observed_cases <- max(observation_df$cases_reported)
  if(max_observed_cases > max_cases)
    stop("Max possible cases should be (much) greater than max observed cases.")
  possible_cases <- max_observed_cases:max_cases
  logps <- conditional_cases_logp(possible_cases, observation_df, cases_history,
                                  Rt, day_onset, serial_parameters, reporting_parameters,
                                  is_gamma_delay,
                                  kappa, is_negative_binomial)
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
max_uncertain_days <- function(p_gamma_cutoff, reporting_parameters, is_gamma_delay) {
  r_location <- reporting_parameters$location
  r_scale <- reporting_parameters$scale
  if(is_gamma_delay)
    days_from_end <- qgamma_mean_sd(p_gamma_cutoff, r_location, r_scale)
  else
    days_from_end <- qlnorm(p_gamma_cutoff, r_location, r_scale)
  days_from_end
}

#' Draws a possible history (or histories) of cases
#'
#' The distribution being drawn from at each onset time t is given by:
#' \deqn{p(cases_true_t|data, Rt, reporting_params, serial_params) \propto p(data|cases_true, reporting_params)
#'  p(cases_true_t|cases_true_t_1, cases_true_t_2, ..., Rt, serial_params)}
#'
#' @param observation_onset_df a tibble with four columns: time_onset, time_reported, cases_reported, reporting_piece_index
#' @param reporting_parameters a tibble with three columns: reporting_piece_index,
#' mean, sd.
#' @inheritParams observed_cases
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
  maximise=FALSE,
  is_gamma_delay=TRUE,
  kappa=NULL,
  is_negative_binomial=FALSE) {

  if(!"reporting_piece_index" %in% colnames(observation_onset_df))
    stop("observation_onset_df must contain a column: 'reporting_piece_index'.")

  # TODO make this uncertain period determination better; put into separate function
  latest_onset_time <- max(observation_onset_df$time_onset)
  df_latest_onset_time <- observation_onset_df %>%
    dplyr::filter(time_onset==latest_onset_time)
  reporting_index_latest_onset_time <- df_latest_onset_time$reporting_piece_index[1]
  reporting_latest_onset_time <- reporting_parameters %>%
    dplyr::filter(reporting_piece_index == reporting_index_latest_onset_time)
  reporting_latest_onset_time <- list(location=reporting_latest_onset_time$location[1],
                                      scale=reporting_latest_onset_time$scale[1])
  uncertain_period <- max_uncertain_days(p_gamma_cutoff, reporting_latest_onset_time, is_gamma_delay)

  start_uncertain_period <- latest_onset_time - uncertain_period

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
      dplyr::filter(.data$time_onset==onset_time)
    reporting_index <- snapshots_at_onset_time_df$reporting_piece_index[1]
    reporting_parameters_tmp <- reporting_parameters %>%
      dplyr::filter(reporting_piece_index == reporting_index)
    reporting_parameters_tmp <- list(
      location=reporting_parameters_tmp$location,
      scale=reporting_parameters_tmp$scale
      )
    snapshots_at_onset_time_df <- snapshots_at_onset_time_df %>%
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
      reporting_parameters=reporting_parameters_tmp,
      ndraws=1,
      maximise=maximise,
      is_gamma_delay=is_gamma_delay,
      kappa=kappa,
      is_negative_binomial=is_negative_binomial)
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
#' @param nresamples number of importance resamples of Rt to perform if assuming a
#' negative binomial model
#'
#' @return a draw (or draws) for Rt
#' @importFrom rlang .data
sample_Rt_single_piece <- function(
    Rt_piece_index, cases_history_df,
    Rt_prior_parameters, serial_parameters,
    kappa=NULL,
    serial_max=40, ndraws=1,
    maximise=FALSE) {

  short_df <- cases_history_df %>%
    dplyr::filter(.data$Rt_index <= Rt_piece_index)
  time_max_post_initial_period <- max(short_df$time_onset) - serial_max

  posterior_shape <- Rt_prior_parameters$shape
  posterior_rate <- Rt_prior_parameters$rate

  short_df <- short_df %>%
    dplyr::mutate(time_after_start = .data$time_onset - serial_max)
  onset_times <- short_df %>%
    dplyr::filter(.data$Rt_index == Rt_piece_index) %>%
    dplyr::pull(.data$time_onset)

  w <- weights_series(serial_max, serial_parameters)
  for(i in seq_along(onset_times)) {

    onset_time <- onset_times[i]
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

  sample_or_maximise_gamma(
    posterior_shape, posterior_rate, ndraws, maximise)
}


#' Sample piecewise-constant Rt values
#'
#' If the renewal model is specified by a Poisson (the default):
#' \deqn{cases_true_t ~ Poisson(Rt * \sum_tau=1^t_max w_t cases_true_t-tau)}
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
#' Alternatively the renewal equation may be specified by a negative binomial distribution:
#' \deqn{cases_true_t ~ NB(Rt * \sum_tau=1^t_max w_t cases_true_t-tau, kappa)}
#' where kappa is the overdispersion parameter. In this case, importance sampling
#' using the Poisson posterior as the importance distribution is used to
#' estimate a negative binomial posterior.
#'
#' @inheritParams sample_Rt_single_piece
#' @return a tibble with three columns: "Rt_piece_index", "draw_index", "Rt"
#' @export
sample_Rt <- function(cases_history_df,
                      Rt_prior_parameters,
                      serial_parameters,
                      kappa=NULL,
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
      kappa,
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

#' Prepares data into list for Stan file and also outputs initialisation function
#'
#' @param df
#' @param current_values
#' @param priors
#' @param is_negative_binomial
#' @param is_rw_prior
#' @param serial_parameters
#' @param serial_max
#' @param n_iterations
#' @param is_gamma_delay
#'
#' @return a named list of 'data' and 'init'
prepare_stan_data_and_init <- function(
    df,
    current_values,
    priors,
    is_negative_binomial,
    is_rw_prior,
    serial_parameters,
    serial_max,
    n_iterations,
    is_gamma_delay) {

  df_short <- df %>%
    dplyr::select(time_onset, Rt_index, cases_true) %>%
    unique()
  w <- weights_series(serial_max, serial_parameters)

  prior_Rt <- priors$Rt
  is_gamma_Rt_prior <- prior_Rt$is_gamma
  priors_reporting <- priors$reporting

  current_values_R <- current_values$R
  theta <- current_values$theta

  if(is_negative_binomial) {

    if("overdispersion" %in% names(priors)) {
      prior_overdispersion <- priors$overdispersion
      prior_kappa_a <- prior_overdispersion$location
      prior_kappa_b <- prior_overdispersion$scale
    } else{
      prior_kappa_a <- 0
      prior_kappa_b <- 10
    }

    current_value_overdispersion <- current_values$kappa

    init_fn <- function() {
      list(
        R=current_values_R,
        theta=theta,
        kappa=as.array(c(current_value_overdispersion))
      )
    }

  } else {
    prior_kappa_a <- 0
    prior_kappa_b <- 10

    init_fn <- function() {
      list(
        R=current_values_R,
        theta=theta
      )
    }

  }

  data_stan <- list(

    # data for renewal model
    N=nrow(df_short),
    K=max(df_short$Rt_index),
    window=df_short$Rt_index,
    C=df_short$cases_true,
    wmax=serial_max,
    w=w,

    # data for reporting delay model
    N_delay=nrow(df),
    time_reported=df$time_reported,
    time_onset=df$time_onset,
    cases_reported=df$cases_reported,
    cases_true=df$cases_true,
    n_reporting_window=max(df$reporting_piece_index),
    reporting_window=df$reporting_piece_index,

    # options
    is_poisson=ifelse(is_negative_binomial==1, 0, 1),
    is_gamma_reporting_delay=ifelse(is_gamma_delay, 1, 0),
    prior_R_choice=ifelse(is_gamma_Rt_prior, 2, 1),
    prior_R_a=prior_Rt$location,
    prior_R_b=prior_Rt$scale,
    prior_kappa_a=prior_kappa_a,
    prior_kappa_b=prior_kappa_b,
    prior_theta_1_a=priors_reporting$mean_mu,
    prior_theta_1_b=priors_reporting$mean_sigma,
    prior_theta_2_a=priors_reporting$sd_mu,
    prior_theta_2_b=priors_reporting$sd_sigma
  )

  list(data=data_stan, init=init_fn)
}

#' Initialises Stan model to allow access to log_prob of model
#'
#' Needs to be rerun each time any data being passed to the Stan file changes
#'
#' @inheritParams prepare_stan_data_and_init
#' @param stan_model an object generated by running cmdstanr::cmdstan_model
#'
#' @return an initialised model
stan_initialisation <- function(
    df,
    current_values,
    priors,
    is_negative_binomial,
    is_rw_prior,
    serial_parameters,
    serial_max,
    is_gamma_delay,
    stan_model) {

  print("buttt")
  tmp <- prepare_stan_data_and_init(
    df,
    current_values,
    priors,
    is_negative_binomial,
    is_rw_prior,
    serial_parameters,
    serial_max,
    n_iterations,
    is_gamma_delay)

  data_stan <- tmp$data
  init_fn <- tmp$init

  print(init_fn())

  # model must be run for one iteration to give access to log-prob
  model <- stan_model$sample(
    data=data_stan,
    init = init_fn,
    iter_warmup = 1,
    iter_sampling = 0,
    adapt_delta=0.9,
    refresh = 0,
    show_messages = FALSE,
    show_exceptions = FALSE,
    chains=1)
  model$init_model_methods()
  model
}


#' Runs single Adaptive-Covariance Metropolis step
#'
#' @param current a named with with current parameter values
#' @param model a stan object generated by stan_initialisation
#' @param mu a tuned hyperparameter vector
#' @param omega a tuned hyperparameter covariance matrix
#' @param log_lambda a tuned hyperparameter value
#' @param eta a hyperparameter that when larger leads to slower decay in adaptations
#' @param iteration iteration count
#'
#' @return a named list with elements: 'current', 'log_p_current', 'mu', 'omega', 'log_lambda'
metropolis_step_Rt_reporting_overdispersion <- function(current, model, mu, omega, log_lambda, eta, iteration) {

  unconstrained_current <- model$unconstrain_variables(current)
  unconstrained_proposed <- mvtnorm::rmvnorm(1, unconstrained_current, exp(log_lambda) * omega)
  constrained_proposed <- model$constrain_variables(unconstrained_proposed)

  proposed <- list(R=constrained_proposed$R,
                   theta=constrained_proposed$theta,
                   kappa=constrained_proposed$kappa,
                   sigma=constrained_proposed$sigma)

  # needs recalculation since cases likely changed
  log_p_current <- model$log_prob(unconstrained_current)
  log_p_proposed <- model$log_prob(unconstrained_proposed)

  log_r <- log_p_proposed - log_p_current

  log_u <- log(runif(1))
  if(log_r > log_u) {
    current <- proposed
    unconstrained_current <- unconstrained_proposed
    alpha <- 1
    log_p_current <- log_p_proposed
  } else {
    alpha <- 0
  }

  # adaptive covariance updates, from https://github.com/pints-team/pints/blob/main/pints/_mcmc/_haario_bardenet_ac.py
  gamma <- (iteration + 1)^(-eta)
  mu <- (1 - gamma) * mu + gamma * unconstrained_current
  m1 <- matrix(unconstrained_current - mu, nrow=length(unconstrained_current))
  m_extra <- gamma * m1 %*% t(m1)
  omega <- (1 - gamma) * omega + m_extra
  log_lambda <- log_lambda + gamma * (alpha - 0.234)

  list(current=current,
       log_p_current=log_p_current,
       mu=mu,
       omega=omega,
       log_lambda=log_lambda)
}


#' Maximises log-probability using L-BFGS-B
#'
#' @param current a named list containing the current parameter values
#' @param model an initialised Stan model
#'
#' @return a named list with a single element: 'current'
maximise_Rt_reporting_overdispersion <- function(current, model) {

  uncons_current <- model$unconstrain_variables(current)
  neg_log_like <- function(theta) {
    -model$log_prob(theta)
  }

  # gradient something led to errors being thrown so not using in optimisation
  # grad_neg_log_like <- function(theta) {
  #   -model$grad_log_prob(theta)
  # }
  # optim minimises by default (note, tried Stan's opt here but it fails very often)
  opt <- optim(uncons_current, neg_log_like, method="L-BFGS-B")
  unconstrained_proposed <- opt$par
  proposed <- model$constrain_variables(unconstrained_proposed)
  list(current=proposed)
}

#' Runs MCMC or optimisation to estimate Rt, cases and reporting parameters
#'
#' @param niterations number of MCMC iterations to run or number of iterative maximisations to run
#' @param snapshot_with_Rt_index_df a tibble with four columns: time_onset, time_reported,
#' cases_reported, Rt_index; an optional fifth column can be provided named reporting_piece_index,
#' which specifies which reporting distribubtion each onset time corresponds to
#' @param priors a named list with: 'Rt', 'reporting', 'max_cases' (and optionally
#' 'overdispersion' if using a negative binomial model). These take
#' the form: 'Rt' is a named list with elements 'location', 'scale' and 'is_gamma'. If
#' is_gamma is TRUE, then an independent gamma prior is placed on each Rt value; if
#' is_gamma is FALSE, a random walk prior is specified on the Rt values;
#' reporting' is a named list with elements
#' 'mean_mu', 'mean_sigma', 'sd_mu', 'sd_sigma' representing the
#' prior parameters for the mean and sd parameters of the reporting parameter distribution
#' (itself described by a gamma distribution by default and a lognormal if is_gamma_delay is FALSE);
#' max cases controls the upper limit of the discrete uniform distribution representing the prior on true
#' cases; 'overdispersion' is a named list specifying the mean and sd of a gamma
#' prior on this parameter
#' @inheritParams sample_Rt_single_piece
#' @inheritParams max_uncertain_days
#' @param reporting_metropolis_parameters named list of 'mean_step', 'sd_step' containing
#' step sizes for Metropolis step
#' @param maximise whether to estimate MAP values of parameters (if true) or
#' sample parameter values using MCMC (if false). By default this is false.
#' @param initial_cases_true a tibble with two columns: "time_onset" and "cases_true", which represents initial
#' estimates of the true number of cases with each onset time.
#' @param initial_reporting_parameters provides initial guesses of the mean and
#' sd of the reporting delay distribution(s). These can be either a named with two named elements
#' ('mean', 'sd') for a time-invariant reporting delay or a tibble with three columns:
#' 'reporting_piece_index', 'mean', 'sd' (where the number of indices corresponds to the number
#' provided in the data frame).
#' @param initial_Rt initial guess of the Rt values in each of the piecewise segments.
#' Provided in the form of a tibble with columns: 'Rt_index' and 'Rt'
#' @param print_to_screen prints progress of MCMC sampling to screen. Defaults to true.
#' @param initial_overdispersion the initial value of the overdispersion parameter if
#' assuming a negative binomial sampling model (default to 5)
#' @param initial_sigma the initial value of the sigma parameter that characterises
#' the random walk prior if selected for Rt
#' @param is_negative_binomial Boolean indicating whether to use a negative binomial renewal model
#' or a Poisson model (if FALSE). By default a Poisson model is assumed
#' @param is_gamma_delay Boolean indicating whether a gamma reporting delay distribution
#' is assumed (if TRUE) or a lognormal (if FALSE). By default, a gamma distribution is assumed
#' @return a named list of three tibbles: "cases", "Rt" and "reporting" which contain estimates of the model parameters
#' @importFrom rlang .data
#' @importFrom methods is
#' @export
mcmc_single <- function(
  niterations,
  snapshot_with_Rt_index_df,
  priors,
  serial_parameters,
  initial_cases_true,
  initial_reporting_parameters,
  initial_Rt,
  stan_model,
  initial_overdispersion=NULL,
  initial_sigma=NULL,
  is_negative_binomial=FALSE,
  is_gamma_delay=TRUE,
  serial_max=40, p_gamma_cutoff=0.99, maximise=FALSE, print_to_screen=TRUE) {

  cnames <- colnames(snapshot_with_Rt_index_df)
  expected_names <- c("time_onset", "time_reported",
                      "cases_reported", "Rt_index")

  if(sum(cnames %in% expected_names) != 4)
    stop("Incorrect column names in snapshot_with_Rt_index_df")

  if(length(cnames) != 4 & length(cnames) != 5)
    stop("There should be either four or five columns in snapshot_with_Rt_index_df")

  if(length(cnames) == 5 & !("reporting_piece_index" %in% cnames))
    stop("The reporting delay indices must be provided through a column named 'reporting_piece_index'")

  # assume if column not supplied then reporting delay same throughout
  if(!"reporting_piece_index" %in% cnames)
    snapshot_with_Rt_index_df$reporting_piece_index <- 1

  if(is_negative_binomial) {

    if(is.null(initial_overdispersion))
      stop("Initial overdispersion value (for a negative binomial renewal model) must not be NULL.")

    if(initial_overdispersion <= 0)
      stop("Initial overdispersion value (for a negative binomial renewal model) must be positive.")

    overdispersion_current <- initial_overdispersion
  }

  is_rw_Rt_prior <- FALSE
  if(!priors$Rt$is_gamma) { # denotes RW prior on Rt

    is_rw_Rt_prior <- TRUE

    if(is.null(initial_sigma))
      stop("Initial sigma value must not be NULL.")

    if(initial_sigma <= 0)
      stop("Initial sigma value must be positive.")

    sigma_current <- initial_sigma
  }

  reporting_current <- initial_reporting_parameters
  if(methods::is(reporting_current, "list")) { # if only a single list provided
    reporting_current <- dplyr::tibble(
      reporting_piece_index=1,
      location=reporting_current$location,
      scale=reporting_current$scale
    )
  }
  num_reporting_parameters <- max(reporting_current$reporting_piece_index)

  df_running <- snapshot_with_Rt_index_df %>%
    dplyr::left_join(initial_cases_true, by = "time_onset") %>%
    dplyr::left_join(initial_Rt, by = "Rt_index")

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
      dplyr::select(time_onset, time_reported, cases_reported, reporting_piece_index)

    df_temp <- sample_cases_history(
      df_current, max_cases,
      Rt_function, serial_parameters,
      reporting_current,
      p_gamma_cutoff=p_gamma_cutoff,
      maximise=maximise,
      is_gamma_delay=is_gamma_delay,
      kappa=overdispersion_current,
      is_negative_binomial=is_negative_binomial)

    cases_history_df <- df_running %>%
      dplyr::select(time_onset, Rt_index) %>%
      dplyr::left_join(df_temp, by = "time_onset", relationship = "many-to-many") %>%
      dplyr::select(-c("cases_reported", "time_reported")) %>%
      dplyr::rename(cases_true=cases_estimated) %>%
      unique()

    df_Rt_index <- df_running %>%
      dplyr::select(time_onset, Rt_index) %>%
      unique()

    df_temp <- df_temp %>%
      dplyr::left_join(df_Rt_index, by="time_onset")

    current_Rt <- df_running %>%
      dplyr::select(Rt_index, Rt) %>%
      unique()

    theta <- matrix(c(reporting_current$location,
                      reporting_current$scale),
                    ncol = 2)
    current <- list(R=current_Rt$Rt,
                    theta=theta)
    if(is_negative_binomial)
      current$kappa <- overdispersion_current
    if(is_rw_Rt_prior)
      current$sigma <- sigma_current

    print(current)

    # just runs model for one iteration to give access to log_p
    model <- stan_initialisation(
      df_temp %>%
        dplyr::rename(cases_true=cases_estimated),
      current,
      priors,
      is_negative_binomial,
      is_rw_prior,
      serial_parameters,
      serial_max,
      is_gamma_delay,
      stan_model)

    print("biya")

    if(i == 1 & !maximise) { # setup parameters for adaptive covariance MCMC

      mu <- model$unconstrain_variables(current)
      omega <- diag(1.0, nrow = length(mu), ncol = length(mu))
      log_lambda <- 0
      eta <- 0.6 # standard value used in adaptive covariance
    }

    if(!maximise) { # MCMC

      res <- metropolis_step_Rt_reporting_overdispersion(
        current, model,
        mu, omega, log_lambda,
        eta, i)
      mu <- res$mu
      omega <- res$omega
      log_lambda <- res$log_lambda

    } else {

      res <- maximise_Rt_reporting_overdispersion(
        current, model)

    }

    current <- res$current
    overdispersion_current <- current$kappa
    sigma_current <- current$sigma

    df_Rt <- dplyr::tibble(Rt=current$R, Rt_index=seq_along(Rt))
    df_running <- df_temp %>%
      dplyr::left_join(df_Rt, by="Rt_index")

    reporting_tmp <- dplyr::tibble(
      location=current$theta[, 1],
      scale=current$theta[, 2],
      reporting_piece_index=seq_along(location),
      iteration=i
    )
    reporting_current <- reporting_tmp %>%
      dplyr::filter(iteration == i)

    if(i == 1)
      reporting_samples <- reporting_tmp
    else
      reporting_samples <- reporting_samples %>% bind_rows(reporting_tmp)

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

    if(is_negative_binomial) {

      overdipersion_tmp <- dplyr::tibble(
        overdispersion=current$kappa,
        iteration=i)

      if(i == 1)
        overdispersion_samples <- overdipersion_tmp
      else
        overdispersion_samples <- overdispersion_samples %>%
        dplyr::bind_rows(overdipersion_tmp)
    }

    if(is_rw_Rt_prior) {

      sigma_tmp <- dplyr::tibble(
        sigma=current$sigma,
        iteration=i)

      if(i == 1)
        sigma_samples <- sigma_tmp
      else
        sigma_samples <- sigma_samples %>%
          dplyr::bind_rows(sigma_tmp)
    }

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

  cases_samples <- cases_history_samples %>%
    dplyr::select(-reporting_piece_index)

  reporting_samples <- reporting_samples %>%
    dplyr::relocate("reporting_piece_index", "location", "scale", "iteration")

  list_results <-list(
    cases=cases_samples,
    Rt=Rt_samples,
    reporting=reporting_samples)

  if(is_negative_binomial) {
    if(!is_rw_Rt_prior)
      list_results$other <- overdispersion_samples
    else
      list_results$other <- overdispersion_samples %>% left_join(sigma_samples, by="iteration")
  } else{
    if(is_rw_Rt_prior)
      list_results$other <- sigma_samples
  }

  list_results
}

#' Combines Markov chains across multiple runs of mcmc_single
#'
#' @param list_of_results a list of results, where each element is a result of running mcmc_single
#' @param is_either_nb_or_sigma Boolean indicating whether model is negative binomial or includes sigma parameter
#' @return a named list of three tibbles: "cases", "Rt" and "reporting" (and if running
#' a negative binomial model an additional 'overdispersion' element) which
#' contain estimates of the model parameters with chain index  included
combine_chains <- function(list_of_results, is_either_nb_or_sigma=FALSE) {

  for(i in seq_along(list_of_results)) {
    res <- list_of_results[[i]]
    cases_df <- res$cases
    Rt_df <- res$Rt
    reporting_df <- res$reporting
    cases_df$chain <- i
    Rt_df$chain <- i
    reporting_df$chain <- i

    if(is_either_nb_or_sigma) {
      other_df <- res$other
      other_df$chain <- i
    }

    if(i == 1) {
      cases_overall <- cases_df
      Rt_overall <- Rt_df
      reporting_overall <- reporting_df
      if(is_either_nb_or_sigma) {
        other_overall <- other_df
      }
    } else {
      cases_overall <- cases_overall %>% dplyr::bind_rows(cases_df)
      Rt_overall <- Rt_overall %>% dplyr::bind_rows(Rt_df)
      reporting_overall <- reporting_overall %>% dplyr::bind_rows(reporting_df)

      if(is_either_nb_or_sigma) {
        other_overall <- other_overall %>% dplyr::bind_rows(other_df)
      }

    }
  }
  list_combined <- list(
    cases=cases_overall,
    Rt=Rt_overall,
    reporting=reporting_overall
  )

  if(is_either_nb_or_sigma)
    list_combined$other <- other_overall

  list_combined
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
#' @inheritParams mcmc_single
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
    stan_model,
    nchains=1,
    initial_overdispersion=NULL,
    initial_sigma=NULL,
    is_negative_binomial=FALSE,
    is_gamma_delay=TRUE,
    serial_max=40, p_gamma_cutoff=0.99, maximise=FALSE, print_to_screen=TRUE,
    is_parallel=FALSE) {

  is_either_nb_or_sigma <- is_negative_binomial | !is.null(initial_sigma)

  if(nchains==1) {
    res <- mcmc_single(niterations,
                       data,
                       priors,
                       serial_parameters,
                       initial_cases_true,
                       initial_reporting_parameters,
                       initial_Rt,
                       stan_model,
                       initial_overdispersion,
                       initial_sigma,
                       is_negative_binomial,
                       is_gamma_delay,
                       serial_max,
                       p_gamma_cutoff,
                       maximise,
                       print_to_screen)

    res$cases$chain <- 1
    res$Rt$chain <- 1
    res$reporting$chain <- 1
    if(is_either_nb_or_sigma)
      res$other$chain <- 1
  } else {

    list_of_results <- vector(mode = "list", length = nchains)
    if(is_parallel) {

      if (requireNamespace("foreach", quietly = TRUE, .packages = "Rcpp", .noexport = c(stan_model))) {

        f_run_single <- function() {
          mcmc_single(niterations,
                      data,
                      priors,
                      serial_parameters,
                      initial_cases_true,
                      initial_reporting_parameters,
                      initial_Rt,
                      stan_model,
                      initial_overdispersion,
                      initial_sigma,
                      is_negative_binomial,
                      is_gamma_delay,
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
                           data,
                           priors,
                           serial_parameters,
                           initial_cases_true,
                           initial_reporting_parameters,
                           initial_Rt,
                           stan_model,
                           initial_overdispersion,
                           initial_sigma,
                           is_negative_binomial,
                           is_gamma_delay,
                           serial_max,
                           p_gamma_cutoff,
                           maximise,
                           print_to_screen)
        list_of_results[[i]] <- res
      }
    }

    res <- combine_chains(list_of_results, is_either_nb_or_sigma)
  }

  res
}
