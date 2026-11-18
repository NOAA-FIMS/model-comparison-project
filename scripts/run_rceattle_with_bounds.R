#' Fit one Rceattle recruitment scenario
#'
#' @param data_list Rceattle data_list.
#' @param inits Initialized parameter list from seed_rceattle_inits().
#' @param scenario One of RCEATTLE_SCENARIOS.
#' @param srr_pars Beverton-Holt alpha / beta from rceattle_bh_pars().
#' @param initMode Initial-age-structure mode. 0 ("FreeParams") estimates the
#'   year-1 age structure as free parameters, which is what WHAM
#'   (`N1_model = 0`, 12 free age-specific N1) and FIMS (`log_init_naa` as fixed
#'   effects) both do.
#' @param time_limit Wall-clock seconds after which a fit is abandoned (returns
#'   NULL, recorded as non-converged). Guards against the pathological
#'   random-effects fit on a deterministic OM case (true sigmaR = 0), where the
#'   Laplace inner problem is ill-conditioned and the optimiser can stall. The
#'   budget covers the estimateMode = 3 map build as well as the fit itself.
#' @return A fitted object of class "Rceattle" (with `$run_time` populated by a
#'   wrapping system.time()), or NULL on error / time-out.
#' @keywords internal
#' @noRd
fit_rceattle_scenario <- function(data_list, inits, scenario, srr_pars,
                                  initMode = 0, time_limit = 300) {
  # Set at function scope, not inside the tryCatch expression below: an on.exit
  # registered from a lazily-evaluated tryCatch argument attaches to tryCatch's
  # frame rather than this one, so the limit could outlive the call.
  setTimeLimit(elapsed = time_limit, transient = TRUE)
  on.exit(setTimeLimit(), add = TRUE)
  
  random_rec <- scenario %in%
    c("random_effects", "random_effects_sigmaR_constant")
  
  # Optimiser controls. getJointPrecision is off because nothing downstream
  # reads it. A loopnum restart begins where the last pass stopped, so it does
  # not rescue a failed fit; the starting values in seed_rceattle_inits() do.
  ctl <- Rceattle::fit_control(
    getsd = TRUE, 
    getJointPrecision = FALSE,
    loopnum = 1, 
    newtonsteps = 3, 
    phase = FALSE, 
    bias_adjust_obs = FALSE,   # the OM draws observations with no bias correction
    nlminb_control = list(iter.max = 10000, eval.max = 10000, trace = 0),
    verbose = 1
  )
  
  # Beverton-Holt with alpha fixed at the value implied by the OM's steepness;
  # srr_est_mode = 0 is Rceattle's "fix alpha to the prior mean". WHAM does the
  # same thing explicitly, fixing log_SR_a from the true h and Phi.0 and
  # estimating only beta (run_wham.R). FIMS holds steepness constant too, but at
  # its own package default of 0.75 rather than at the OM's value -- the two
  # coincide for these cases and would not for an OM with a different steepness.
  #
  # Estimating alpha and beta jointly is not an option here: it is singular on
  # this data (reciprocal condition number ~2e-39).
  common <- list(
    data_list  = data_list,
    inits      = inits,
    msmMode    = 0,                                  # single species
    random_rec = random_rec,
    recFun     = Rceattle::build_srr(
      srr_fun        = 2,                            # Beverton-Holt
      srr_est_mode   = 0,                            # fix alpha at the prior mean
      srr_prior      = srr_pars$alpha,
      srr_alpha_init = srr_pars$alpha,
      srr_beta_init  = srr_pars$beta),
    initMode   = initMode
  )
  
  # Build parameter bounds (and map) directly from Rceattle's dry-run mode
  # (estimateMode = 3) so no external reference fit is required.
  build <- tryCatch(
    do.call(Rceattle::fit_mod, c(common, list(
      estimateMode = 3,
      fit_control  = Rceattle::fit_control(getsd = FALSE, verbose = 0)))),
    error = function(e) {
      message(sprintf("  [%s] bounds/map build error: %s", scenario, conditionMessage(e)))
      NULL
    }
  )
  if (is.null(build)) return(NULL)
  
  # Floor active init_dev slots (1:11, which parameterize Ages 2:12) at
  # log(0.01) = -4.605 so N_{a,1} >= 0.01 (safely above Rceattle's 0.001 floor),
  # while preserving the -999 sentinel in the unused Age12 slot.
  bounds <- build$bounds
  if (!is.null(bounds$lower$init_dev)) {
    active_init <- (inits$init_dev != -999)
    bounds$lower$init_dev[active_init] <- log(0.01)
  }
  
  # Parameter map. fit_mod() builds the map internally when `map = NULL`, which
  # is what the two default scenarios want (random_effects: R_log_sd estimated;
  # fixed_effects: random_rec = FALSE already maps R_log_sd off and treats
  # rec_dev as penalised fixed effects). For the sigmaR-constant scenario we
  # keep rec_dev as random effects but hold R_log_sd fixed at its initial value
  # log(sigma_rec_prior).
  map <- NULL
  if (scenario == "random_effects_sigmaR_constant") {
    map <- build$map
    map$mapFactor$R_log_sd <- factor(rep(NA, length(map$mapFactor$R_log_sd)))
    if (!is.null(map$mapList)) map$mapList$R_log_sd[] <- NA
  }
  
  fit <- NULL
  fit_without_sdreport <- NULL
  # create ctl_without_sdreport by setting getsd = FALSE in ctl
  ctl_without_sdreport <- ctl
  ctl_without_sdreport$getsd <- FALSE
  fit_without_sdreport <- tryCatch(
    do.call(Rceattle::fit_mod, c(common, list(
      map          = map,
      bounds = bounds,
      estimateMode = 1,                        # hindcast (no HCR projection)
      fit_control  = ctl_without_sdreport))),
    error = function(e) {
      message(sprintf("  [%s] fit_mod error: %s", scenario, conditionMessage(e)))
      NULL
    }
  )
  
  elapsed <- system.time({
    fit <- tryCatch(
      do.call(Rceattle::fit_mod, c(common, list(
        map          = map,
        bounds = bounds,
        estimateMode = 1,                        # hindcast (no HCR projection)
        fit_control  = ctl))),
      error = function(e) {
        message(sprintf("  [%s] fit_mod error: %s", scenario, conditionMessage(e)))
        NULL
      }
    )
  })
  
  if (is.null(fit)) return(NULL)
  
  rceattle_obj <- fit$obj
  tmb_sdreport_time <- system.time(
    TMB::sdreport(rceattle_obj, 
                  skip.delta.method = TRUE,
                  getReportCovariance = FALSE)
  )
  
  if (is.null(fit_without_sdreport)) {
    warning(sprintf(
      "[%s] fit without sdreport failed; optimization runtime is unavailable.",
      scenario
    ), call. = FALSE)
    attr(fit, "fit_optimization") <- NA_real_
  } else {
    # run_time is a difftime, which R reports in minutes once past 60 s.
    attr(fit, "fit_optimization") <- as.numeric(fit_without_sdreport$run_time, units = "secs")
  }
  
  attr(fit, "tmb_sdreport") <- tmb_sdreport_time[["elapsed"]]
  fit
}