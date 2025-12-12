# insert etiology regression for latent status code chunk into .bug file; discrete predictors

insert etiology regression for latent status code chunk into .bug file;
discrete predictors

## Usage

``` r
insert_bugfile_chunk_reg_discrete_predictor_etiology(Jcause, ppd = NULL)
```

## Arguments

- Jcause:

  The number of distinct causes, i.e., categories of latent health
  status; equals `length(model_options$likelihood$cause_list)`.

- ppd:

  Default is NULL; set to TRUE for posterior predictive checking

## Value

a long character string to be inserted into .bug model file as
distribution specification for latent status
