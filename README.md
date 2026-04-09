# Propensity Weight Analysis

This repository contains R code for performing propensity score weighting and consent weighting on multiply imputed datasets.

## Overview

The analysis pipeline consists of 10 modular, verifiable steps:

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

## Requirements

```r
library(data.table)
library(mice)
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