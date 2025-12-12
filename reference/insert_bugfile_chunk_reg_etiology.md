# insert etiology regression for latent status code chunk into .bug file

insert etiology regression for latent status code chunk into .bug file

## Usage

``` r
insert_bugfile_chunk_reg_etiology(Eti_formula, Jcause, ppd = NULL)
```

## Arguments

- Eti_formula:

  Etiology regression formula; Check
  `model_options$likelihood$Eti_formula`.

- Jcause:

  The number of distinct causes, i.e., categories of latent health
  status; equals `length(model_options$likelihood$cause_list)`.

- ppd:

  Default is NULL; set to TRUE for posterior predictive checking

## Value

a long character string to be inserted into .bug model file as
distribution specification for latent status
