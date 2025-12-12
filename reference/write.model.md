# function to write bugs model (copied from R2WinBUGS)

function to write bugs model (copied from R2WinBUGS)

## Usage

``` r
write.model(model, con = "model.bug", digits = 5)
```

## Arguments

- model:

  R / S-PLUS function containing the BUGS model in the BUGS model
  language, for minor differences see Section Details.

- con:

  passed to writeLines which actually writes the model file

- digits:

  number of significant digits used for WinBUGS input, see formatC

## Value

write text lines to a connection; same as
[`writeLines()`](https://rdrr.io/r/base/writeLines.html)
