# baker: **B**ayesian **A**nalytic **K**it for **E**tiology **R**esearch

`baker` is designed for disease etiology studies from case-control data
with multiple sources of measurements with potential errors. If you are
interested in estimating the population etiology pie (a vector of
fractions that sum to one), and the probability of each cause for a
particular individual case, try `baker`.

## Value

No returned value; documentation purpose only.

## Details

`baker` implements hierarchical Bayesian models to infer disease
etiology for multivariate binary data. We created `baker` to catalyze
effective communications between analysts and practicing clinicians that
are vital to the success of etiology studies. The `baker` package offers
modules to

- Import and tidy the PERCH data (the study that motivates the creation
  of this package),

- Transform, explore the data,

- Specify, automatically generate the model files, and fit the models
  (npLCM),

- Store and visualize posterior summaries for communicating scientific
  findings, and

- Check and compare the fitted models.

`baker` has implemented models for dependent measurements given disease
status, regression analyses of etiology, multiple imperfect
measurements, different priors for true positive rates among cases with
differential measurement characteristics, and multiple-pathogen
etiology. Scientists in Pneumonia Etiology Research for Child Health
(PERCH) study usually refer to the etiology distribution as "population
etiology pie" and "individual etiology pie" for their compositional
nature, hence the name of the package (baking the pie).

## baker functions

[`nplcm()`](https://zhenkewu.com/baker/reference/nplcm.md)

## See also

- <https://github.com/zhenkewu/baker> for the source code and
  system/software requirements to use `baker` for your data.

## Author

**Maintainer**: Zhenke Wu <zhenkewu@gmail.com>
([ORCID](https://orcid.org/0000-0001-7582-669X)) \[copyright holder\]

Authors:

- Scott Zeger <sz@jhu.edu>
  ([ORCID](https://orcid.org/0000-0001-8907-1603))

Other contributors:

- John Muschelli <muschellij2@gmail.com>
  ([ORCID](https://orcid.org/0000-0001-6469-1750)) \[contributor\]

- Irena Chen <irena@umich.edu>
  ([ORCID](https://orcid.org/0000-0002-9366-8506)) \[contributor\]
