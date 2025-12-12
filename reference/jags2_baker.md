# Run `JAGS` from R

The jags function takes data and starting values as input. It
automatically writes a jags script, calls the model, and saves the
simulations for easy access in R. Check the R2jags::jags2 for details
about the argument.

## Usage

``` r
jags2_baker(
  data,
  inits,
  parameters.to.save,
  model.file = "model.bug",
  n.chains = 3,
  n.iter = 2000,
  n.burnin = floor(n.iter/2),
  n.thin = max(1, floor((n.iter - n.burnin)/1000)),
  DIC = TRUE,
  jags.path = "",
  working.directory = NULL,
  clearWD = TRUE,
  refresh = n.iter/50
)
```

## Arguments

- data:

  \(1\) a vector or list of the names of the data objects used by the
  model, (2) a (named) list of the data objects themselves, or (3) the
  name of a "dump" format file containing the data objects, which must
  end in ".txt", see example below for details.

- inits:

  a list with `n.chains` elements; each element of the list is itself a
  list of starting values for the `BUGS` model, *or* a function creating
  (possibly random) initial values. If inits is `NULL`, `JAGS` will
  generate initial values for parameters.

- parameters.to.save:

  character vector of the names of the parameters to save which should
  be monitored.

- model.file:

  file containing the model written in `BUGS` code. Alternatively, as in
  R2WinBUGS, `model.file` can be an R function that contains a `BUGS`
  model that is written to a temporary model file (see `tempfile`) using
  `write.model`

- n.chains:

  number of Markov chains (default: 3)

- n.iter:

  number of total iterations per chain (including burn in; default:
  2000)

- n.burnin:

  length of burn in, i.e. number of iterations to discard at the
  beginning. Default is `n.iter/2`, that is, discarding the first half
  of the simulations. If n.burnin is 0, `jags()` will run 100 iterations
  for adaption.

- n.thin:

  thinning rate. Must be a positive integer. Set `n.thin` \> 1 to save
  memory and computation time if `n.iter` is large. Default is
  `max(1, floor(n.chains * (n.iter-n.burnin) / 1000))` which will only
  thin if there are at least 2000 simulations.

- DIC:

  logical; if `TRUE` (default), compute deviance, pD, and DIC. The rule
  `pD=var(deviance) / 2` is used.

- jags.path:

  directory that contains the `JAGS` executable. The default is “”.

- working.directory:

  sets working directory during execution of this function; This should
  be the directory where model file is.

- clearWD:

  indicating whether the files `data.txt`, `inits[1:n.chains].txt`,
  `codaIndex.txt`, `jagsscript.txt`, and `CODAchain[1:nchains].txt`
  should be removed after `jags` has finished, default=TRUE.

- refresh:

  refresh frequency for progress bar, default is `n.iter/50`

## Value

Same as [`R2jags::jags()`](https://rdrr.io/pkg/R2jags/man/jags.html)

## Details

This modifies the jags2 function in R2jags package.

## See also

[`R2jags::jags()`](https://rdrr.io/pkg/R2jags/man/jags.html)
