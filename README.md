# Propensity Weight Analysis

This repository contains R code for performing propensity score weighting and consent weighting on multiply imputed datasets, followed by longitudinal mixed model analysis.

## Overview

The analysis pipeline consists of 13 modular, verifiable steps:

### Weighting Pipeline (Steps 1-10)
1. **Extract imputed datasets** - Convert mice object to list of data.tables
2. **Add consent indicator** - Mark participants who consented for data linkage
3. **Filter pre-COVID** - Keep only pre-COVID observations
4. **Filter mental health sample** - Select mentally distressed (GHQ > 10) or MH service users
5. **Process HSU** - Recode HSU to binary and remove HSU cases
6. **Create interaction term** - Calculate age × GHQ interaction
7. **Fit propensity model** - Logistic regression for treatment propensity
8. **Fit consent model** - Logistic regression for consent weighting (IPCW)
9. **Calculate stabilized IPTW** - Inverse probability of treatment weights
10. **Combine weights** - Multiply IPTW by IPCW

### Analysis Pipeline (Steps 11-13)
11. **Reshape to long format** - Wide to long transformation for longitudinal analysis
12. **Fit mixed models** - Linear mixed effects models across imputed datasets
13. **Pool results** - Combine estimates using Rubin's rules

## Requirements

```r
library(data.table)
library(mice)
library(lme4)
library(pbapply)
library(optimx)
```

## Usage

### Step-by-step verification

```r
source("propensity_weight_analysis.R")

# Step 1: Extract datasets
implist <- extract_imputed_datasets(imp1)
print(paste("Number of datasets:", length(implist)))

# Step 2: Add consent indicator
implist <- lapply(implist, add_consent_indicator, consent_ids = consent)

# Step 3-10: Continue with remaining steps...
# (see propensity_weight_analysis.R for full example)
```

### Full pipeline

```r
source("propensity_weight_analysis.R")

implist <- run_full_pipeline(imp1, consent)
```

## Longitudinal Analysis (Steps 11-13)

After computing weights, reshape and fit mixed models:

### Step-by-step

```r
# Step 11: Reshape to long format
implong <- lapply(implist, reshape_wide_to_long)
print(head(implong[[1]]))

# Step 12: Fit mixed models
formula <- ghq ~ 1 + wave * t0mhs + (1 | pidp)
# Or quadratic: ghq ~ 1 + (wave + wave2) * t0mhs + (1 | pidp)
models <- fit_mixed_models(implong, formula)

# Step 13: Pool results
pooled_results <- pool_mixed_model_results(models)
print(pooled_results)
```

### Full analysis pipeline

```r
formula <- ghq ~ 1 + wave * t0mhs + (1 | pidp)
results <- run_analysis_pipeline(implist, formula)
print(results$pooled_results)
```

## Customizing Formulas

Both propensity and consent models accept custom formulas:

```r
# Custom propensity formula
custom_ps_formula <- t0mhs ~ age + sf12_pcs9 + sf12_mcs9 + lt_sick9
implist <- lapply(implist, fit_propensity_model, formula = custom_ps_formula)

# Custom consent formula
custom_consent_formula <- consent ~ age + totincome9 + imd_mean
implist <- lapply(implist, fit_consent_model, formula = custom_consent_formula)
```

## Output Variables

After running the pipeline, each dataset includes:

| Variable | Description |
|----------|-------------|
| `consent` | Consent indicator (1 = consented, 0 = not) |
| `int1` | Age × GHQ interaction term |
| `pst0` | Propensity score for T0 treatment |
| `p_consent` | Predicted probability of consent |
| `pconsent` | Marginal probability of consent |
| `ipcw` | Inverse probability of consent weight |
| `pt0` | Marginal probability of T0 treatment |
| `iptwt0` | Stabilized IPTW for T0 |
| `iptwt0t1` | Combined weight (IPTW × IPCW) |

### After reshape (long format)

| Variable | Description |
|----------|-------------|
| `wave` | Time wave (0, 1, 2) |
| `wave2` | Wave squared (for quadratic models) |
| `ghq` | GHQ score at each wave |
| `isolation` | Social isolation at each wave |
| `mcs` | SF-12 mental component score at each wave |
| `pcs` | SF-12 physical component score at each wave |
| `life` | Life satisfaction at each wave |
| `health` | Health satisfaction at each wave |