# Insert measurement likelihood (nested model+regression) code chunks into .bug model file

Insert measurement likelihood (nested model+regression) code chunks into
.bug model file

## Usage

``` r
insert_bugfile_chunk_reg_nest_meas(
  Mobs,
  prior,
  cause_list,
  FPR_formula,
  use_measurements = "BrS",
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

A long character string to be inserted into .bug model file as
measurement likelihood

## See also

Called by
[write_model_Reg_NoNest](https://zhenkewu.com/baker/reference/write_model_Reg_NoNest.md)
for constructing a `.bug` file. This is usually called along with
specification of latent status regression
([insert_bugfile_chunk_reg_etiology](https://zhenkewu.com/baker/reference/insert_bugfile_chunk_reg_etiology.md)).
