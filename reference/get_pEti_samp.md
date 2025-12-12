# get etiology samples by names (no regression)

get etiology samples by names (no regression)

## Usage

``` r
get_pEti_samp(res_nplcm, model_options)
```

## Arguments

- res_nplcm:

  result from model fits

- model_options:

  model specification

## Value

A list:

- `pEti_mat`: a matrix of posterior samples (iteration by cause);
  overall etiology `latent_nm`: a vector of character strings
  representing the names of the causes
