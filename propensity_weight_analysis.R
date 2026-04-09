# =============================================================================
# Propensity Score Weighting and Consent Weighting for Multiply Imputed Data
# =============================================================================
# This script performs propensity score weighting and consent weighting
# on multiply imputed datasets in separate, verifiable steps.
#
# Required packages: data.table, mice, lme4, pbapply, optimx
# =============================================================================

library(data.table)
library(mice)
library(lme4)
library(pbapply)
library(optimx)

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
# STEP 11: Reshape wide format to long format
# =============================================================================
# Transpose wide format data into long format for longitudinal analysis
# Input: dt (data.table in wide format)
# Output: data.table in long format with wave variable

reshape_wide_to_long <- function(dt) {
  dt <- copy(dt)
  setDT(dt)
  
  # Patterns use ^ and $ anchors to match exact column names only
  # (e.g., matches "ghq_scale9" but not "ghq_scale9_extra")
  dflong <- melt(
    dt,
    id.vars = c("pidp", "iptwt0", "iptwt0t1", "hsu", "age", "t0mhs", "sex", "consent"),
    measure = patterns(
      ghq = "^ghq_scale9$|^ghq_scale10$|^ghq_scale11$",
      isolation = "^scisolate9$|^scisolate10$|^scisolate11$",
      mcs = "^sf12_mcs9$|^sf12_mcs10$|^sf12_mcs11$",
      pcs = "^sf12_pcs9$|^sf12_pcs10$|^sf12_pcs11$",
      life = "^life_satisf9$|^life_satisf10$|^life_satisf11$",
      health = "^health_satisf9$|^health_satisf10$|^health_satisf11$"
    ),
    variable.name = "wave"
  )
  
  # Convert wave from factor (1, 2, 3) to numeric (0, 1, 2) for modeling
  # Wave 0 = baseline (wave 9), Wave 1 = follow-up 1, Wave 2 = follow-up 2
  dflong[, wave := as.numeric(wave) - 1]
  dflong[, wave2 := wave^2]
  
  return(dflong[])
}

# Usage:
# implong <- lapply(implist, reshape_wide_to_long)

# =============================================================================
# STEP 12: Fit mixed models across imputed datasets
# =============================================================================
# Run linear mixed effects models on each imputed dataset
# Input: data_list (list of data.tables), formula (lmer formula),
#        use_consent_only (whether to subset to consenters),
#        weights_col (column name for weights)
# Output: list of lmer model objects

fit_mixed_models <- function(data_list, 
                             formula, 
                             use_consent_only = TRUE,
                             weights_col = "iptwt0t1") {
  
  model_outputs <- pblapply(data_list, function(dt) {
    # Subset to consenters if requested
    if (use_consent_only) {
      dt <- dt[consent == 1]
    }
    
    # Fit mixed model with weights
    mod <- lmer(
      formula, 
      data = dt,
      weights = dt[[weights_col]],
      control = lmerControl(
        optimizer = "optimx", 
        optCtrl = list(method = "nlminb")
      )
    )
    
    return(mod)
  })
  
  return(model_outputs)
}

# Usage:
# formula <- ghq ~ 1 + wave * t0mhs + (1 | pidp)
# formula <- ghq ~ 1 + (wave + wave2) * t0mhs + (1 | pidp)
# output <- fit_mixed_models(implong, formula)

# =============================================================================
# STEP 13: Pool results using Rubin's rules
# =============================================================================
# Pool mixed model results across imputed datasets
# Input: model_list (list of lmer model objects)
# Output: pooled results summary with confidence intervals

pool_mixed_model_results <- function(model_list) {
  # Convert to mira object for pooling
  mira_object <- as.mira(model_list)
  
  # Pool using Rubin's rules and return summary with CIs
  pooled <- pool(mira_object)
  result <- summary(pooled, conf.int = TRUE)
  
  return(result)
}

# Usage:
# pooled_results <- pool_mixed_model_results(output)
# print(pooled_results)

# =============================================================================
# ANALYSIS PIPELINE (Steps 11-13)
# =============================================================================
# Run longitudinal analysis pipeline after weighting steps

run_analysis_pipeline <- function(implist, 
                                  formula,
                                  use_consent_only = TRUE,
                                  weights_col = "iptwt0t1") {
  
  message("Step 11: Reshaping to long format...")
  implong <- lapply(implist, reshape_wide_to_long)
  
  message("Step 12: Fitting mixed models...")
  models <- fit_mixed_models(
    implong, 
    formula, 
    use_consent_only = use_consent_only,
    weights_col = weights_col
  )
  
  message("Step 13: Pooling results...")
  pooled_results <- pool_mixed_model_results(models)
  
  message("Analysis complete!")
  return(list(
    long_data = implong,
    models = models,
    pooled_results = pooled_results
  ))
}

# Usage:
# formula <- ghq ~ 1 + wave * t0mhs + (1 | pidp)
# results <- run_analysis_pipeline(implist, formula)
# print(results$pooled_results)

# =============================================================================
# STEP 14: Run multi-outcome analysis
# =============================================================================
# Run mixed model analysis for multiple outcomes at once
# Input: implong (list of data.tables in long format), 
#        outcomes (character vector of outcome variable names),
#        predictors (character string of RHS of formula without outcome),
#        use_consent_only (whether to subset to consenters),
#        weights_col (column name for weights)
# Output: named list of pooled results for each outcome

run_multi_outcome_analysis <- function(implong,
                                       outcomes = c("ghq", "mcs", "pcs", "life"),
                                       predictors = "1 + wave * t0mhs + (1 | pidp)",
                                       use_consent_only = TRUE,
                                       weights_col = "iptwt0t1") {
  
  results <- lapply(setNames(outcomes, outcomes), function(outcome) {
    message(paste("Analyzing outcome:", outcome))
    
    # Build formula dynamically
    formula <- as.formula(paste(outcome, "~", predictors))
    
    # Fit models across imputed datasets
    models <- fit_mixed_models(
      implong, 
      formula, 
      use_consent_only = use_consent_only,
      weights_col = weights_col
    )
    
    # Pool results and add outcome column
    pooled <- pool_mixed_model_results(models)
    pooled <- as.data.frame(pooled)
    pooled$outcome <- outcome
    
    return(list(
      formula = formula,
      models = models,
      pooled = pooled
    ))
  })
  
  message("All outcomes analyzed!")
  return(results)
}

# Usage:
# implong <- lapply(implist, reshape_wide_to_long)
# results <- run_multi_outcome_analysis(implong)
# results$ghq$pooled   # GHQ results
# results$mcs$pooled   # MCS results
# results$pcs$pooled   # PCS results
# results$life$pooled  # Life satisfaction results

# =============================================================================
# STEP 15: Combine and export pooled results
# =============================================================================
# Combine pooled results from multiple outcomes into a single data.frame
# Input: multi_results (output from run_multi_outcome_analysis)
# Output: combined data.frame with outcome column

combine_pooled_results <- function(multi_results) {
  combined <- do.call(rbind, lapply(names(multi_results), function(name) {
    df <- multi_results[[name]]$pooled
    df$outcome <- name
    return(df)
  }))
  
  # Reorder columns to put outcome first
  cols <- c("outcome", setdiff(names(combined), "outcome"))
  combined <- combined[, cols]
  
  return(combined)
}

# Usage:
# combined_results <- combine_pooled_results(results)
# print(combined_results)

# =============================================================================
# FULL MULTI-OUTCOME PIPELINE
# =============================================================================
# Complete pipeline from weighted data to pooled multi-outcome results

run_full_analysis <- function(implist,
                              outcomes = c("ghq", "mcs", "pcs", "life"),
                              predictors = "1 + wave * t0mhs + (1 | pidp)",
                              use_consent_only = TRUE,
                              weights_col = "iptwt0t1") {
  
  message("Step 11: Reshaping to long format...")
  implong <- lapply(implist, reshape_wide_to_long)
  
  message("Step 14: Running multi-outcome analysis...")
  multi_results <- run_multi_outcome_analysis(
    implong,
    outcomes = outcomes,
    predictors = predictors,
    use_consent_only = use_consent_only,
    weights_col = weights_col
  )
  
  message("Step 15: Combining results...")
  combined <- combine_pooled_results(multi_results)
  
  message("Analysis complete!")
  return(list(
    long_data = implong,
    by_outcome = multi_results,
    combined = combined
  ))
}

# Usage:
# # Linear model
# results <- run_full_analysis(implist)
# print(results$combined)
#
# # Quadratic model
# results_quad <- run_full_analysis(
#   implist,
#   predictors = "1 + (wave + wave2) * t0mhs + (1 | pidp)"
# )
# print(results_quad$combined)

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
