# make a mapping template for model fitting

`make_template` creates a mapping matrix (binary values). Each pathogen
in a measurement slice (e.g., nasal-pharyngeal PCR test) is mapped to
inform one category of latent status. All the possible categories (e.g.,
causes of pneumonia) remain the same regardless of the measurement slice
used (e.g., NPPCR or BCX).

## Usage

``` r
make_template(patho, cause_list)
```

## Arguments

- patho:

  A vector of pathogen names for a particular measurement slice. `patho`
  must be a substring of some elements in `cause_list`, e.g., "PNEU" is
  a substring of "PNEU_VT13". Also see Examples for this function.

- cause_list:

  A vector of characters; Potential categories of latent statuses.

## Value

a mapping from `patho` to `cause_list`. `NROW = length(cause_list)+1`;
`NCOL = length(patho)`. This value is crucial in model fitting to
determine which measurements are informative of a particular category of
latent status.

## Details

The first argument has to be character substrings from the second
argument. For example, the two arguments can respectively be `"A"` and
`"A_1"`, or `"A"` and `"A+B"`.The second argument can have character
strings not matched in the first argument. If so, it means some causes
of diseases are not directly measured in the current measurement slice.
For each element of `patho`, the function matches from the start of the
strings of `cause_list`. Therefore, make sure that latent statuses from
the same family (e.g., "PNEU_VT13" and "PNEU_NOVT13") need to start with
the same family name (e.g., "PNEU") followed by subcategories (e.g.,
"\_VT13" and "\_NOVT13").

## Examples

``` r
cause_list <- c("HINF","PNEU_VT13","PNEU_NOVT13","SAUR","HMPV_A_B","FLU_A",
"PARA_1","PARA_3","PARA_4","PV_EV","RHINO","RSV", "ENTRB","TB")

patho_BrS_NPPCR <- c("HINF","PNEU","SAUR","HMPV_A_B","FLU_A","PARA_1",
"PARA_3","PARA_4","PV_EV","RHINO","RSV")
make_template(patho_BrS_NPPCR,cause_list)
#>       [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10] [,11]
#>  [1,]    1    0    0    0    0    0    0    0    0     0     0
#>  [2,]    0    1    0    0    0    0    0    0    0     0     0
#>  [3,]    0    1    0    0    0    0    0    0    0     0     0
#>  [4,]    0    0    1    0    0    0    0    0    0     0     0
#>  [5,]    0    0    0    1    0    0    0    0    0     0     0
#>  [6,]    0    0    0    0    1    0    0    0    0     0     0
#>  [7,]    0    0    0    0    0    1    0    0    0     0     0
#>  [8,]    0    0    0    0    0    0    1    0    0     0     0
#>  [9,]    0    0    0    0    0    0    0    1    0     0     0
#> [10,]    0    0    0    0    0    0    0    0    1     0     0
#> [11,]    0    0    0    0    0    0    0    0    0     1     0
#> [12,]    0    0    0    0    0    0    0    0    0     0     1
#> [13,]    0    0    0    0    0    0    0    0    0     0     0
#> [14,]    0    0    0    0    0    0    0    0    0     0     0
#> [15,]    0    0    0    0    0    0    0    0    0     0     0


 cause = c("A","B1","B2","C","A+C","B+C")
 patho = c("A","B","C")
 make_template(patho,cause)
#>      [,1] [,2] [,3]
#> [1,]    1    0    0
#> [2,]    0    1    0
#> [3,]    0    1    0
#> [4,]    0    0    1
#> [5,]    1    0    1
#> [6,]    0    1    1
#> [7,]    0    0    0
 
 cause = c("A","B1","B2","C","A+C","B+C","other")
 patho = c("A","B","C")
 make_template(patho,cause)
#>      [,1] [,2] [,3]
#> [1,]    1    0    0
#> [2,]    0    1    0
#> [3,]    0    1    0
#> [4,]    0    0    1
#> [5,]    1    0    1
#> [6,]    0    1    1
#> [7,]    0    0    0
#> [8,]    0    0    0
 
 
 cause = c("A","B1","B2","X_B","Y_B","C","A+C","B+C","other")
 patho = c("A","B","C","X_B","Y_B")
 make_template(patho,cause)
#>       [,1] [,2] [,3] [,4] [,5]
#>  [1,]    1    0    0    0    0
#>  [2,]    0    1    0    0    0
#>  [3,]    0    1    0    0    0
#>  [4,]    0    0    0    1    0
#>  [5,]    0    0    0    0    1
#>  [6,]    0    0    1    0    0
#>  [7,]    1    0    1    0    0
#>  [8,]    0    1    1    0    0
#>  [9,]    0    0    0    0    0
#> [10,]    0    0    0    0    0
 
 
```
