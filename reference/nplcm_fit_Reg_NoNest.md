# Fit nested partially-latent class model with regression (low-level)

Fit nested partially-latent class model with regression (low-level)

## Usage

``` r
nplcm_fit_Reg_NoNest(data_nplcm, model_options, mcmc_options)
```

## Arguments

- data_nplcm:

  Cases are on top of controls in the rows of diagnostic test results
  and the covariate matrix. This is assumed by `baker` to automatically
  write model files (`.bug`).

  - `Mobs` A list of measurements of distinct qualities (Bronze-,
    Silver, and Gold-Standard: `MBS`,`MSS`,`MGS`). The elements of the
    list should include `MBS`, `MSS`, and `MGS`. If any of the component
    is not available, please specify it as, e.g., `MGS=NULL`
    (effectively deleting `MGS` from `Mobs`).

    - `MBS` a list of data frame of bronze-standard (BrS) measurements.
      For each data frame (referred to as a 'slice'), rows are subjects,
      columns are causative agents (e.g., pathogen species). We use
      `list` here to accommodate the possibility of multiple sets of BrS
      data. They have imperfect sensitivity/specificity (e.g.
      nasopharyngeal polymerase chain reaction - NPPCR).

    - `MSS` a list of data frame of silver-standard (SS) measurements.
      Rows are subjects, columns are causative agents measured in
      specimen (e.g. blood culture). These measurements have perfect
      specificity but imperfect sensitivity.

    - `MGS` a list of data frame of gold-standard (GS) measurements.
      Rows are subject, columns are measured causative agents These
      measurements have perfect sensitivity and specificity.

  - `Y` Vector of disease status: `1` for case, `0` for control.

  - `X` Covariate matrix. A subset of columns are primary covariates in
    cause-specific- case-fraction (CSCF) functions and hence must be
    available for cases, and another subset are covariates that are
    available in the cases and the controls. The two sets of covariates
    may be identical, overlapping or completely different. In general,
    this is not the design matrix for regression models, because for
    enrollment date in a study which may have non-linear effect, basis
    expansion is often needed for approximation.

- model_options:

  A list of model options: likelihood and prior.

  `use_measurements`

  :   A vector of characters strings; can be one or more from `"BrS"`,
      `"SS"`, `"GS"`.

  `likelihood`

  :   

      cause_list

      :   The vector of causes (NB: specify);

      k_subclass

      :   The number of nested subclasses in each disease class (one of
          case classes or the control class; the same `k_subclass` is
          assumed for each class) and each slice of BrS measurements.
          `1` for conditional independence; larger than `1` for
          conditional dependence. It is only available for BrS
          measurements. It is a vector of length equal to the number of
          slices of BrS measurements;

      Eti_formula

      :   Formula for etiology regressions. You can use
          [`s_date_Eti()`](https://zhenkewu.com/baker/reference/s_date_Eti.md)
          to specify the design matrix for `R` format enrollment date;
          it will produce natural cubic spline basis. Specify `~ 1` if
          no regression is intended.

      FPR_formula

      :   formula for false positive rates (FPR) regressions; see
          [`formula()`](https://rdrr.io/r/stats/formula.html). You can
          use
          [`s_date_FPR()`](https://zhenkewu.com/baker/reference/s_date_FPR.md)
          to specify part of the design matrix for `R` format enrollment
          date; it will produce penalized-spline basis (based on
          B-splines). Specify `~ 1` if no regression is intended. (NB:
          If `effect="fixed"`,
          [`dm_Rdate_FPR()`](https://zhenkewu.com/baker/reference/dm_Rdate_FPR.md)
          will just specify a design matrix with appropriately
          standardized dates.)

  `prior`

  :   

      Eti_prior

      :   Description of etiology prior (e.g., `overall_uniform` - all
          hyperparameters are `1`; or `0_1` - all hyperparameters are
          `0.1`);

      TPR_prior

      :   Description of priors for the measurements (e.g., informative
          vs non-informative). Its length should be the same as
          `use_measurements` above. Please see examples for how to
          specify. The package can also handle multiple slices of BrS,
          SS data, so separate specification of the TPR priors are
          needed.

- mcmc_options:

  A list of Markov chain Monte Carlo (MCMC) options.

  - `debugstatus` Logical - whether to pause WinBUGS after it finishes
    model fitting; (NB: is this obsolete? Test.)

  - `n.chains` Number of MCMC chains;

  - `n.burnin` Number of burn-in iterations;

  - `n.thin` To keep every other `n.thin` samples after burn-in period;

  - `individual.pred` `TRUE` to perform individual prediction (`Icat`
    variables in the `.bug` file); `FALSE` otherwise;

  - `ppd` `TRUE` to simulate new data (`XXX.new` variables in the `.bug`
    file) from the posterior predictive distribution (ppd); `FALSE`
    otherwise;

  - `get.pEti` `TRUE` for getting posterior samples of individual
    etiologic fractions; `FALSE` otherwise. For non-regression, or
    regression models with all discrete predictors, by default this is
    `TRUE`, so no need to specify this entry. It is only relevant for
    regression models with non-discrete covariates. Because individuals
    have distinct CSCFs at their specific covariate values, it's easier
    to just store the posterior samples of the regression coefficients
    and reconstruct the pies afterwards, rather than storing them
    through `JAGS`.

  - `result.folder` Path to folder storing the results;

  - `bugsmodel.dir` Path to `.bug` model files;

  - `jags.dir` Path to where JAGS is installed; if `NULL`, this will be
    set to `jags.dir=""`.

## Value

BUGS fit results from `JAGS`.

## Details

This function prepares data, specifies hyperparameters in priors (true
positive rates and CSCFs), initializes the posterior sampling chain,
writes the model file (for JAGS or WinBUGS with slight differences in
syntax), and fits the model. Features:

- regression (not all discrete covariates);

- no nested subclasses, i.e. conditional independence of multivariate
  measurements given disease class and covariates;

- multiple BrS + multiple SS.

## See also

[write_model_NoReg](https://zhenkewu.com/baker/reference/write_model_NoReg.md)
for constructing `.bug` model file; This function then puts it in the
folder `mcmc_options$bugsmodel.dir`.

Other model fitting functions:
[`nplcm_fit_NoReg()`](https://zhenkewu.com/baker/reference/nplcm_fit_NoReg.md),
[`nplcm_fit_Reg_Nest()`](https://zhenkewu.com/baker/reference/nplcm_fit_Reg_Nest.md),
[`nplcm_fit_Reg_discrete_predictor_NoNest()`](https://zhenkewu.com/baker/reference/nplcm_fit_Reg_discrete_predictor_NoNest.md)
