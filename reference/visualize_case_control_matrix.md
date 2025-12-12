# Visualize matrix for a quantity measured on cases and controls (a single number)

Special to case-control visualization: upper right for cases, lower left
for controls.

## Usage

``` r
visualize_case_control_matrix(
  mat,
  dim_names = ncol(mat),
  cell_metrics = "",
  folding_line = TRUE,
  axes = FALSE,
  xlab = "",
  ylab = "",
  asp = 1,
  title = ""
)
```

## Arguments

- mat:

  matrix of values: upper for cases, lower for controls;

- dim_names:

  names of the columns, from left to right. It is also the names of the
  rows, from bottom to top. Default is 1 through `ncol(mat)`;

- cell_metrics:

  the meaning of number in every cell;

- folding_line:

  Default is `TRUE` for adding dashed major diagonal line.

- axes:

  plot axes; default is `FALSE`;

- xlab:

  label for x-axis

- ylab:

  label for y-axis

- asp:

  aspect ratio; default is `1` to ensure square shape

- title:

  text for the figure

## Value

plotting function; no returned value.
