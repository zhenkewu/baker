# Initialize individual latent status (for `JAGS`)

Initialize individual latent status (for `JAGS`)

## Usage

``` r
init_latent_jags_multipleSS(
  MSS_list,
  cause_list,
  patho = unlist(lapply(MSS_list, colnames))
)
```

## Arguments

- MSS_list:

  A list of silver-standard measurement data, possibly with more than
  one slices; see `data_nplcm` argument in
  [`nplcm()`](https://zhenkewu.com/baker/reference/nplcm.md)

- cause_list:

  See `model_options` arguments in
  [`nplcm()`](https://zhenkewu.com/baker/reference/nplcm.md)

- patho:

  A vector of measured pathogen name for `MSS`; default is
  `colnames(MSS)`

## Value

a list of numbers, indicating categories of individual latent causes.

## Details

In `JAGS` 3.4.0, if an initial value contradicts the probabilistic
specification, e.g. `MSS_1[i,j] ~ dbern(mu_ss_1[i,j])`, where
`MSS_1[i,j]=1` but `mu_ss_1[i,j]=0`, then `JAGS` cannot understand it.
In PERCH application, this is most likely used when the specificity of
the silver-standard data is `1`. Note: this is not a problem in
`WinBUGS`.
