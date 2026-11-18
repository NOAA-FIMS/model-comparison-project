maindir <- maindir
subdir <- "Rceattle"
oms_sim_num <- sim_input[["om_sim_num"]]
case_name <- sim_input[["case_name"]]
casedir <- file.path(maindir, case_name)
em_bias_cor <- sim_input[["em_bias_cor"]]
survey_units <- "numbers"
save_full_fit <- TRUE
time_limit <-300
# Scenario identifiers, reused for the file names.
RCEATTLE_SCENARIOS <- c(
  "random_effects",
  "random_effects_sigmaR_constant",
  "fixed_effects"
)

om_sim_list <- c(81, 82, 88, 218)

log_file <- file.path(casedir, "output", subdir, "rceattle_random_effects_console_log.txt")
con <- file(log_file, open = "wt")

# Redirect both standard output (cat/print) and messages/warnings to the file
# split = TRUE lets cat()/print() still appear on your console while logging
sink(con, split = TRUE)
sink(con, type = "message")

tryCatch({
  for (om_sim in om_sim_list) {
    cat("\n#### sim_id = ", om_sim, "\n")
    setup <- tryCatch({
      load(file = file.path(casedir, "output", "OM", paste0("OM", om_sim, ".RData")))
      data_list <- ASSAMC:::om_to_rceattle(om_input, om_output, em_input,
                                           survey_units = survey_units)
      list(
        outdir    = file.path(casedir, "output", subdir, paste0("s", om_sim)),
        data_list = data_list,
        inits     = ASSAMC:::seed_rceattle_inits(data_list, om_input),
        srr_pars  = ASSAMC:::rceattle_bh_pars(om_input),
        nyrs      = length(om_input[["year"]])
      )
    }, error = function(e) {
      message(sprintf("  [OM %d] skipped: %s", om_sim, conditionMessage(e)))
      conditionMessage(e)
    })
    
    if (is.character(setup)) next
    
    msgs <- character(0)
    for (scenario in "random_effects") {
      res <- tryCatch({
        fit <- :fit_rceattle_scenario(
        #fit <- ASSAMC:::fit_rceattle_scenario(
          data_list  = setup$data_list,
          inits      = setup$inits,
          scenario   = scenario,
          srr_pars   = setup$srr_pars,
          time_limit = time_limit
        )
        if (!is.null(fit) && scenario == "fixed_effects") {
          write.csv(
            fit$opt$par,
            row.names = names(fit$opt$par),
            file = file.path(casedir, "output", subdir, paste0("s", om_sim),
                             "parameters_fixed_effects_rceattle.csv")
          )
        }
        ASSAMC:::save_rceattle_outputs(fit, scenario, setup$outdir, nyrs = setup$nyrs,
                                       save_full_fit = save_full_fit)
        NA_character_
      }, error = function(e) {
        message(sprintf("  [OM %d / %s] failed: %s",
                        om_sim, scenario, conditionMessage(e)))
        sprintf("%s (%s)", conditionMessage(e), scenario)
      })
      if (!is.na(res)) msgs <- c(msgs, res)
    }
    if (length(msgs)) {
      message(sprintf("OM %d: %s", om_sim, paste(msgs, collapse = "; ")))
    }
  }
}, finally = {
  # Always restore console output/messages and close the file connection
  sink(type = "message")
  sink()
  close(con)
})
