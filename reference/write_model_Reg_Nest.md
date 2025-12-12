# Write `.bug` model file for regression model WITH nested subclasses

`write_model_Reg_Nest` automatically generates model file according to
`model_options`; This is called within
[nplcm_fit_Reg_Nest](https://zhenkewu.com/baker/reference/nplcm_fit_Reg_Nest.md).

## Usage

``` r
write_model_Reg_Nest(
  Mobs,
  prior,
  cause_list,
  Eti_formula,
  FPR_formula,
  use_measurements,
  ppd = NULL,
  use_jags = FALSE
)
```

## Arguments

- Mobs:

  Measurement data in the form of `data_nplcm`

- prior:

  Prior specification from `model_options`

- cause_list:

  A list of latent status names (crucial for building templates; see
  [`make_template()`](https://zhenkewu.com/baker/reference/make_template.md))

- Eti_formula:

  Etiology regression formula; Check
  `model_options$likelihood$Eti_formula`.

- FPR_formula:

  A list of FPR regression formula; check
  `model_options$likelihood$FPR_formula`

- use_measurements:

  "BrS", or "SS"

- ppd:

  Default is NULL; set to TRUE for posterior predictive checking

- use_jags:

  Default is FALSE; set to TRUE if want to use JAGS for model fitting.

## Value

a long character string to be written into .bug file.

## See also

[insert_bugfile_chunk_noreg_meas](https://zhenkewu.com/baker/reference/insert_bugfile_chunk_noreg_meas.md)
for inserting .bug file chunk for measurements (plug-and-play.R);
[insert_bugfile_chunk_noreg_etiology](https://zhenkewu.com/baker/reference/insert_bugfile_chunk_noreg_etiology.md)
for inserting .bug file chunk for distribution of latent status
(etiology).

Other model generation functions:
[`write_model_NoReg()`](https://zhenkewu.com/baker/reference/write_model_NoReg.md),
[`write_model_Reg_NoNest()`](https://zhenkewu.com/baker/reference/write_model_Reg_NoNest.md),
[`write_model_Reg_discrete_predictor_NoNest()`](https://zhenkewu.com/baker/reference/write_model_Reg_discrete_predictor_NoNest.md)
