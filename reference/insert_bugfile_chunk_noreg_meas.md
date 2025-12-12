# Insert measurement likelihood (without regression) code chunks into .bug model file

Insert measurement likelihood (without regression) code chunks into .bug
model file

## Usage

``` r
insert_bugfile_chunk_noreg_meas(
  k_subclass,
  Mobs,
  prior,
  cause_list,
  use_measurements = "BrS",
  ppd = NULL,
  use_jags = FALSE
)
```

## Arguments

- k_subclass:

  the number of subclasses for the slices that require conditional
  dependence modeling (only applicable to BrS data); its length is of
  the same value as the number of BrS slices.

- Mobs:

  measurement data in the form of `data_nplcm`

- prior:

  prior specification from `model_options`

- cause_list:

  a list of latent status names (crucial for building templates; see
  [`make_template()`](https://zhenkewu.com/baker/reference/make_template.md))

- use_measurements:

  "BrS", or "SS"

- ppd:

  Default is NULL; set to TRUE for posterior predictive checking

- use_jags:

  Default is FALSE; set to TRUE if want to use JAGS for model fitting.

## Value

a long character string to be inserted into .bug model file as
measurement likelihood

## See also

It is used in
[write_model_NoReg](https://zhenkewu.com/baker/reference/write_model_NoReg.md)
for constructing a .bug file along with specification of latent status
distribution
([insert_bugfile_chunk_noreg_etiology](https://zhenkewu.com/baker/reference/insert_bugfile_chunk_noreg_etiology.md))
