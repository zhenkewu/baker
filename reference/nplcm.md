# Fit nested partially-latent class models (highest-level wrapper function)

Uses `JAGS` (OSX or Windows) operating system for Bayesian posterior
inference (see `README` file for an instruction to install `JAGS`). If
running `JAGS` on windows, please go to control panel to add the
directory to `JAGS` into ENVIRONMENTAL VARIABLE.

## Usage

``` r
nplcm(data_nplcm, model_options, mcmc_options)
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

A `JAGS` output result, fitted by function
[`R2jags::jags2()`](https://rdrr.io/pkg/R2jags/man/jags.html) from
`R2jags`. It is an object of class `nplcm` and `bugs`. Current
implemented models follow the hierarchy below:

- no regression: Fitted by at low level by
  [nplcm_fit_NoReg](https://zhenkewu.com/baker/reference/nplcm_fit_NoReg.md)

- regression: Given disease class (control or a class of cases with the
  same subset of causative agents):

  - local independence model for BrS measures: Fitted at lower level by

    - [nplcm_fit_Reg_NoNest](https://zhenkewu.com/baker/reference/nplcm_fit_Reg_NoNest.md)
      deals with the setting with two sets of covariates, one for CSCF
      regression and the other for FPR regression. The two sets of
      covariates may be identical, overlapping or non-overlapping. This
      function is called when there exists one or more than one discrete
      covariate among the union of the two covariate sets. The method
      implemented by this function directly lets FPR depend upon
      covariates. This is different from Wu and Chen (2021), which let
      the subclass weights depend upon covariates. We implemented this
      function for methods comparison.

    - [nplcm_fit_Reg_discrete_predictor_NoNest](https://zhenkewu.com/baker/reference/nplcm_fit_Reg_discrete_predictor_NoNest.md)
      deals with the setting with all discrete covariates for FPRs and
      CSCFs. The strata defined by the two sets of covariates need not
      be identical, e.g., as a result of distinct sets of covariates.
      Again, this is directly to let FPR be stratified by covariates,
      hence different from Wu and Chen (2020+) We implemented this
      function for methods comparison.

  - local dependence model for BrS measures: Fitted at lower level by
    [nplcm_fit_Reg_Nest](https://zhenkewu.com/baker/reference/nplcm_fit_Reg_Nest.md):
    This is the method introduced in Wu and Chen (2021): CSCF
    regression + case/control subclass weight regression. It does not
    provide a specialized function for the setting with all discrete
    covariates.

## Examples

``` r
# \donttest{
data(data_nplcm_noreg)
cause_list <- LETTERS[1:6]
J.BrS      <- 6
model_options_no_reg <- list(
  likelihood   = list(
    cause_list = cause_list,
    k_subclass = 2,
    Eti_formula = ~-1, # no covariate for the etiology regression
    FPR_formula = list(
      MBS1 =   ~-1)    # no covariate for the subclass weight regression
  ),
  use_measurements = c("BrS"), 
  # use bronze-standard data only for model estimation.
  prior= list(
    Eti_prior = overall_uniform(1,cause_list), 
    # Dirichlet(1,...,1) prior for the etiology.
    TPR_prior  = list(BrS = list(
      info  = "informative", # informative prior for TPRs
      input = "match_range", 
      # specify the informative prior for TPRs by specifying a plausible range.
      val = list(MBS1 = list(up =  list(rep(0.99,J.BrS)), 
                             # upper ranges: matched to 97.5% quantile of a Beta prior
                             low = list(rep(0.55,J.BrS))))
      # lower ranges: matched to 2.5% quantile of a Beta prior
    )
    )
  )
)     


set.seed(1)

run_example <- function(){
# include stratification information in file name:
thedir0    <- paste0(tempdir(),"_no_reg")

# create folders to store the model results 
dir.create(thedir0, showWarnings = FALSE)
# remove the temp folder on exit:
# add = TRUE ensures you don't overwrite other exit handlers
on.exit(unlink(thedir0, recursive = TRUE), add = TRUE) 

result_folder_no_reg <- file.path(thedir0,paste("results",collapse="_"))
thedir <- result_folder_no_reg
dir.create(thedir, showWarnings = FALSE)

# options for MCMC chains:
mcmc_options_no_reg <- list(
  debugstatus = TRUE,
  n.chains = 1,
  n.itermcmc = as.integer(200), 
  n.burnin = as.integer(100), 
  n.thin = 1,
  individual.pred = TRUE, # <- must set to TRUE! <------- NOTE! 
  ppd = FALSE,
  result.folder = thedir,
  bugsmodel.dir = thedir
)

BrS_object_1 <- make_meas_object(patho = LETTERS[1:6], 
                                 specimen = "MBS", test = "1", 
                                 quality = "BrS", cause_list = cause_list)
clean_options <- list(BrS_objects = make_list(BrS_object_1))
# place the nplcm data and cleaning options into the results folder
dput(data_nplcm_noreg,file.path(thedir,"data_nplcm.txt")) 
dput(clean_options, file.path(thedir, "data_clean_options.txt"))

rjags::load.module("glm")

nplcm_noreg <- nplcm(data_nplcm_noreg,model_options_no_reg,mcmc_options_no_reg)

}
run_example()
#> ==[baker] Results stored in: == 
#>  /tmp/RtmpRcxbeQ_no_reg/results 
# }

```
