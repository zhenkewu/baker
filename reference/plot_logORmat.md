# Visualize pairwise log odds ratios (LOR) for data that are available in both cases and controls

Visualize pairwise log odds ratios (LOR) for data that are available in
both cases and controls

## Usage

``` r
plot_logORmat(data_nplcm, pathogen_display, BrS_slice = 1, logOR_rounding = 2)
```

## Arguments

- data_nplcm:

  See
  [`assign_model()`](https://zhenkewu.com/baker/reference/assign_model.md).

- pathogen_display:

  The pathogen vector in desired order for display. It can be of larger
  length than that of `pathogen_BrS`.

- BrS_slice:

  Default is 1 - the set of BrS data to visualize.

- logOR_rounding:

  Rounding number of the log odds ratio. Default is 2.

## Value

Figure of LOR matrix and relevant s.e. and significance information.

## Details

`plot_logORmat` visualizes a matrix of pairwise log odds ratios (LOR)
for cases (upper) and controls (lower). LOR is at the top of the cell.
Below it, its standard error is in smaller type, using the same color as
the LOR. Then the estimate is divided by its standard error. We put the
actual value when the Z-statistics has an absolute value greater than
\$2\$; a plus (red) or minus (blue) if between \$1\$ and \$2\$; blank
otherwise.

## See also

Other exploratory data analysis functions:
[`get_top_pattern()`](https://zhenkewu.com/baker/reference/get_top_pattern.md),
[`show_individual()`](https://zhenkewu.com/baker/reference/show_individual.md),
[`summarize_BrS()`](https://zhenkewu.com/baker/reference/summarize_BrS.md),
[`summarize_SS()`](https://zhenkewu.com/baker/reference/summarize_SS.md),
[`visualize_season()`](https://zhenkewu.com/baker/reference/visualize_season.md)

## Examples

``` r
data(data_nplcm_noreg)
plot_logORmat(data_nplcm_noreg,names(data_nplcm_noreg$Mobs$MBS[[1]]))
#> == Visualizing pairwise log odds ratios for bronze-standard data set:  1 :  MBS1 . ==

```
