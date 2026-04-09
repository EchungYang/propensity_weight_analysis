# =============================================================================
# Propensity Score Weighting and Consent Weighting for Multiply Imputed Data
# =============================================================================
# This script performs propensity score weighting and consent weighting
# on multiply imputed datasets in separate, verifiable steps.
#
# Required packages: data.table, mice
# =============================================================================

library(data.table)
library(mice)

# =============================================================================
# STEP 1: Extract imputed datasets from mice object
# =============================================================================
# Convert mice object to a list of completed datasets
# Input: imp1 (mice mids object)
# Output: implist (list of data.tables)

extract_imputed_datasets <- function(mids_object) {
  m <- mids_object$m
  lapply(1:m, function(i) as.data.table(complete(mids_object, i)))
}

# Usage:
# implist <- extract_imputed_datasets(imp1)

# =============================================================================
# STEP 2: Add consent indicator column
# =============================================================================
# Mutate a column based on whether participants consented for data linkage
# Input: dt (data.table), consent_ids (vector of pidp values who consented)
# Output: data.table with consent column (1 = consented, 0 = not)

add_consent_indicator <- function(dt, consent_ids) {
  dt <- copy(dt)
  dt[, consent := fifelse(pidp %in% consent_ids, 1L, 0L)]
  return(dt)
}

# Usage:
# implist <- lapply(implist, add_consent_indicator, consent_ids = consent)

# =============================================================================
# STEP 3: Filter to pre-COVID data
# =============================================================================
# Keep only pre-COVID observations
# Input: dt (data.table)
# Output: filtered data.table

filter_precovid <- function(dt) {
  dt <- copy(dt)
  dt[precovid == TRUE]
}

# Usage:
# implist <- lapply(implist, filter_precovid)

# =============================================================================
# STEP 4: Filter to mentally distressed or MH service users
# =============================================================================
# Select mentally distressed people (GHQ > 10) or mental health service users
# Input: dt (data.table)
# Output: filtered data.table

filter_mental_health_sample <- function(dt, ghq_threshold = 10) {
  dt <- copy(dt)
  dt[ghq_scale8 > ghq_threshold | t0mhs == 1]
}

# Usage:
# implist <- lapply(implist, filter_mental_health_sample)

# =============================================================================
# STEP 5: Process HSU variable and filter
# =============================================================================
# Recode HSU to binary (1 if > 0, else 0) and remove HSU cases
# Input: dt (data.table)
# Output: filtered data.table with hsu == 0

process_and_filter_hsu <- function(dt) {
  dt <- copy(dt)
  dt[, hsu := fifelse(hsu > 0, 1L, 0L, na = 0L)]
  dt[hsu == 0]
}

# Usage:
# implist <- lapply(implist, process_and_filter_hsu)

# =============================================================================
# STEP 6: Create interaction term
# =============================================================================
# Create age * ghq_scale8 interaction term
# Input: dt (data.table)
# Output: data.table with int1 column

create_interaction_term <- function(dt) {
  dt <- copy(dt)
  dt[, int1 := age * ghq_scale8]
  return(dt)
}

# Usage:
# implist <- lapply(implist, create_interaction_term)

# =============================================================================
# STEP 7: Fit propensity score model for T0 treatment
# =============================================================================
# Fit logistic regression for treatment propensity scores
# Input: dt (data.table), formula (model formula)
# Output: data.table with pst0 (propensity score) column

fit_propensity_model <- function(dt, formula = NULL) {
  dt <- copy(dt)
  
  # Default formula for T0 MHS propensity model

  # Note: ghq_scale8 is intentionally excluded from the propensity model 
  # because it is used as a selection criterion (GHQ > 10) in Step 4.
  # Including it would cause near-collinearity issues. The model uses
  # other health and socioeconomic predictors instead.
  if (is.null(formula)) {
    formula <- t0mhs ~ age + sf12_pcs9 + sf12_mcs9 +
      lt_sick9 + hh_size9 + totincome9 + scisolate9 + finnow9
  }
  
  psfit <- glm(formula, data = dt, family = binomial("logit"))
  dt[, pst0 := psfit$fitted.values]
  
  return(dt)
}

# Usage:
# implist <- lapply(implist, fit_propensity_model)

# =============================================================================
# STEP 8: Fit consent weighting model
# =============================================================================
# Fit logistic regression for consent propensity and calculate IPCW
# Input: dt (data.table), formula (model formula)
# Output: data.table with p_consent, pconsent, and ipcw columns

fit_consent_model <- function(dt, formula = NULL) {
  dt <- copy(dt)
  
  # Default consent formula
  if (is.null(formula)) {
    formula <- consent ~ hiqual9 + totincome9 + imd_mean + age + 
      sf12_pcs9 + health_satisf9
  }
  
  consent_fit <- glm(formula, data = dt, family = binomial())
  
  dt[, p_consent := predict(consent_fit, type = "response")]
  dt[, pconsent := mean(consent, na.rm = TRUE)]
  
  # IPCW for Sample 2 (consented only)
  dt[, ipcw := fifelse(consent == 1, pconsent / p_consent, NA_real_)]
  
  return(dt)
}

# Usage:
# implist <- lapply(implist, fit_consent_model)

# =============================================================================
# STEP 9: Calculate stabilized IPT weights for T0
# =============================================================================
# Calculate stabilized inverse probability of treatment weights
# Input: dt (data.table)
# Output: data.table with pt0 and iptwt0 columns

calculate_stabilized_iptw <- function(dt) {
  dt <- copy(dt)
  
  # Marginal probability of treatment
  dt[, pt0 := mean(t0mhs, na.rm = TRUE)]
  
  # Stabilized IPTW formula:
  # numerator: t0mhs * pt0 + (1 - t0mhs) * (1 - pt0)
  # denominator: t0mhs * pst0 + (1 - t0mhs) * (1 - pst0)
  dt[, iptwt0 := (t0mhs * pt0 + (1 - t0mhs) * (1 - pt0)) /
       (t0mhs * pst0 + (1 - t0mhs) * (1 - pst0))]
  
  return(dt)
}

# Usage:
# implist <- lapply(implist, calculate_stabilized_iptw)

# =============================================================================
# STEP 10: Combine weights (IPTW * IPCW)
# =============================================================================
# Multiply propensity weight by consent weight
# Input: dt (data.table)
# Output: data.table with iptwt0t1 (combined weight) column

combine_weights <- function(dt) {
  dt <- copy(dt)
  dt[, iptwt0t1 := iptwt0 * ipcw]
  return(dt)
}

# Usage:
# implist <- lapply(implist, combine_weights)

# =============================================================================
# FULL PIPELINE (for reference)
# =============================================================================
# Run all steps sequentially with ability to verify at each step

run_full_pipeline <- function(mids_object, consent_ids) {
  
  message("Step 1: Extracting imputed datasets...")
  implist <- extract_imputed_datasets(mids_object)
  
  message("Step 2: Adding consent indicator...")
  implist <- lapply(implist, add_consent_indicator, consent_ids = consent_ids)
  
  message("Step 3: Filtering to pre-COVID data...")
  implist <- lapply(implist, filter_precovid)
  
  message("Step 4: Filtering to mental health sample...")
  implist <- lapply(implist, filter_mental_health_sample)
  
  message("Step 5: Processing and filtering HSU...")
  implist <- lapply(implist, process_and_filter_hsu)
  
  message("Step 6: Creating interaction term...")
  implist <- lapply(implist, create_interaction_term)
  
  message("Step 7: Fitting propensity score model...")
  implist <- lapply(implist, fit_propensity_model)
  
  message("Step 8: Fitting consent model...")
  implist <- lapply(implist, fit_consent_model)
  
  message("Step 9: Calculating stabilized IPTW...")
  implist <- lapply(implist, calculate_stabilized_iptw)
  
  message("Step 10: Combining weights...")
  implist <- lapply(implist, combine_weights)
  
  message("Pipeline complete!")
  return(implist)
}

# =============================================================================
# EXAMPLE USAGE (Step-by-step verification)
# =============================================================================
# 
# # Step 1: Extract datasets
# implist <- extract_imputed_datasets(imp1)
# print(paste("Number of imputed datasets:", length(implist)))
# print(paste("Observations per dataset:", nrow(implist[[1]])))
# 
# # Step 2: Add consent
# implist <- lapply(implist, add_consent_indicator, consent_ids = consent)
# print(table(implist[[1]]$consent))
# 
# # Step 3: Filter pre-COVID
# implist <- lapply(implist, filter_precovid)
# print(paste("After pre-COVID filter:", nrow(implist[[1]])))
# 
# # Step 4: Filter mental health sample
# implist <- lapply(implist, filter_mental_health_sample)
# print(paste("After MH filter:", nrow(implist[[1]])))
# 
# # Step 5: Process HSU
# implist <- lapply(implist, process_and_filter_hsu)
# print(paste("After HSU filter:", nrow(implist[[1]])))
# 
# # Step 6: Create interaction
# implist <- lapply(implist, create_interaction_term)
# print(summary(implist[[1]]$int1))
# 
# # Step 7: Propensity model
# implist <- lapply(implist, fit_propensity_model)
# print(summary(implist[[1]]$pst0))
# 
# # Step 8: Consent model
# implist <- lapply(implist, fit_consent_model)
# print(summary(implist[[1]]$p_consent))
# 
# # Step 9: Stabilized IPTW
# implist <- lapply(implist, calculate_stabilized_iptw)
# print(summary(implist[[1]]$iptwt0))
# 
# # Step 10: Combined weights
# implist <- lapply(implist, combine_weights)
# print(summary(implist[[1]]$iptwt0t1))
# 
# # Or run the full pipeline at once:
# # implist <- run_full_pipeline(imp1, consent)
