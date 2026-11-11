# ==============================================================================
# Model Comparison Project: Estimation Models (EMs) and Configuration Settings
# For inclusion in manuscript.qmd (rendered via knitr::kable or tinytable)
# ==============================================================================

library(tibble)
library(dplyr)
library(knitr)

em_comparison_table <- tibble::tribble(
  ~Category, ~Feature_or_Setting, ~ASAP, ~BAM, ~SS3, ~FIMS, ~Rceattle, ~WHAM,

  # --- Framework & Software ---
  "Platform", "Underlying software engine", 
  "ADMB", "ADMB", "ADMB", "TMB & C++", "TMB / RTMB", "TMB",

  "Platform", "Model execution mode in study", 
  "Executable (.exe via Wine/Linux)", "Executable (.exe via Wine/Linux)", "Executable (.exe via Wine/Linux)", "R package (FIMS::fit_fims)", "R package (Rceattle::fit_mod)", "R package (wham::fit_wham)",

  # --- Recruitment Deviations & State-Space Scenarios ---
  "Recruitment Configuration", "Penalized Maximum Likelihood (PML) scenario", 
  "Fixed effects rec_devs with quadratic penalty (sigmaR = 0.2)", 
  "Fixed effects rec_devs with quadratic penalty (sigmaR = 0.2)", 
  "Fixed effects rec_devs with quadratic penalty (sigmaR = 0.2)", 
  "fixed_effects: rec_devs estimated as fixed effects, log_sd fixed constant", 
  "fixed_effects: random_rec = FALSE, rec_devs as penalized fixed effects", 
  "fixed_effects: log_NAA removed from random effects, sigma fixed",

  "Recruitment Configuration", "State-space: unknown rec SD (estimated sigmaR)", 
  "Not supported / not evaluated", 
  "Not supported / not evaluated", 
  "Not supported / not evaluated", 
  "random_effects: rec_devs as TMB random effects, log_sd estimated", 
  "random_effects: random_rec = TRUE, R_log_sd estimated as random walk/RE SD", 
  "random_effects: NAA_re = list(sigma='rec'), log_NAA_sigma estimated",

  "Recruitment Configuration", "State-space: known rec SD (fixed sigmaR)", 
  "Not supported / not evaluated", 
  "Not supported / not evaluated", 
  "Not supported / not evaluated", 
  "random_effects_sigmaR_constant: rec_devs as random effects, log_sd fixed at OM truth", 
  "random_effects_sigmaR_constant: random_rec = TRUE, R_log_sd mapped to NA (fixed at OM truth)", 
  "random_effects_sigmaR_constant: NAA_re = list(sigma='rec'), log_NAA_sigma mapped to NA (fixed at OM truth)",

  # --- Population Dynamics & Initialization ---
  "Population Dynamics", "Age classes modeled", 
  "Ages 1-12+ (plus group at 12)", "Ages 1-12+ (plus group at 12)", "Ages 1-12+ (configured from age 1)", "Ages 1-12+ (plus group at 12)", "Ages 1-12+ (plus group at 12)", "Ages 1-12+ (plus group at 12)",

  "Population Dynamics", "Stock-recruitment curve", 
  "Beverton-Holt (steepness h = 0.75 fixed, R0 estimated)", 
  "Beverton-Holt (median-unbiased; h = 0.75 fixed, R0 estimated)", 
  "Beverton-Holt (mean-unbiased; h = 0.75 fixed, R0 estimated)", 
  "Beverton-Holt (steepness h = 0.75 constant, R0 estimated)", 
  "Beverton-Holt (R = alpha*S / (1 + beta*S); alpha fixed from OM h, beta estimated)", 
  "Beverton-Holt (R = alpha*S / (1 + beta*S); log_SR_a fixed from OM h, beta estimated)",

  "Population Dynamics", "Initial population (Year 1 N-at-age)", 
  "Deviations from equilibrium N1 based on initial F1", 
  "Equilibrium age structure adjusted for initial F1", 
  "Initial equilibrium F1 with recruitment offset", 
  "Free parameters: log_init_naa initialized at OM truth", 
  "Free parameters: initMode = 0 (12 free age-specific N1)", 
  "Free parameters: N1_model = 0 (12 free age-specific N1)",

  "Population Dynamics", "Recruitment bias adjustment / integration", 
  "None (internal machinery absent)", 
  "Median-to-mean adjustment on equilibrium recruitment", 
  "Internal mean-unbiased bias adjustment ramp", 
  "Laplace approximation marginal likelihood (random effects)", 
  "Laplace approximation marginal likelihood (random effects)", 
  "Laplace approximation marginal likelihood (random effects)",

  # --- Fishery & Survey Settings ---
  "Fishery & Survey", "Fishing mortality (F) formulation", 
  "Apical F (separable F with annual devs)", 
  "Apical F (separable F with annual devs)", 
  "Continuous hybrid F / apical F", 
  "Annual log_Fmort estimated per year", 
  "Annual log_F estimated per year", 
  "Annual log_F_tot estimated per year",

  "Fishery & Survey", "Selectivity form (Fishery & Survey)", 
  "Logistic (2 parameters: mu1 = inflection, mu2 = slope)", 
  "Logistic (2 parameters: x1 = slope, x2 = inflection)", 
  "Logistic (2 parameters: nu1 = inflection, nu2 = ln(19)/slope)", 
  "Logistic (inflection point and slope)", 
  "Logistic (inflection point and log slope)", 
  "Logistic (inflection point and slope)",

  "Fishery & Survey", "Catchability (q) for Survey", 
  "Estimated freely (log scale)", 
  "Estimated freely (analytical / log scale)", 
  "Estimated freely (log scale)", 
  "Estimated freely (log_q)", 
  "Estimated freely (closed-form MLE / log_q)", 
  "Estimated freely (log scale)",

  # --- Observation Likelihoods ---
  "Likelihoods", "Landings / Catch", 
  "Lognormal (CV = 0.05)", "Lognormal (CV = 0.05)", "Lognormal (CV = 0.05)", "Lognormal (CV = 0.05)", "Lognormal (CV = 0.05)", "Lognormal (CV = 0.05)",

  "Likelihoods", "Survey index", 
  "Lognormal (CV = 0.20, numbers)", "Lognormal (CV = 0.20, numbers)", "Lognormal (CV = 0.20, numbers)", "Lognormal (CV = 0.20, biomass)", "Lognormal (CV = 0.20, numbers)", "Lognormal (CV = 0.20, numbers)",

  "Likelihoods", "Age composition", 
  "Multinomial (effN = 200)", "Multinomial (effN = 200)", "Multinomial (effN = 200)", "Multinomial (effN = 200)", "Multinomial (effN = 200)", "Multinomial (effN = 200)"
)

# Save table as CSV for manuscript.qmd
write.csv(
  em_comparison_table,
  file.path("C0", "figure", "em_comparison_table.csv"),
  row.names = FALSE
)
