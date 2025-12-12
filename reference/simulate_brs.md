# Simulate Bronze-Standard (BrS) Data

Simulate Bronze-Standard (BrS) Data

## Usage

``` r
simulate_brs(set_parameter, latent_samples)
```

## Arguments

- set_parameter:

  True model parameters in an npLCM specification:

  `cause_list`

  :   a vector of disease class names among cases (since the causes
      could be multi-agent (e.g., multiple pathogens may cause an
      individual case's pneumonia), so its length could be longer than
      the total number of unique causative agents)

  `etiology`

  :   a vector of proportions that sum to 100 percent

  `pathogen_BrS`

  :   a vector of putative causative agents' names measured in
      bronze-standard (BrS) data. This function simulates only one slice
      defined by ``` specimen``test``pathogen ```

  `pathogen_SS`

  :   a vector of pathogen names measured in silver-standard (SS) data.

  `meas_nm`

  :   a list of ``` specimen``test ``` names e.g.,
      `list(MBS = c("NPPCR"),MSS="BCX")` for nasopharyngeal (NP)
      specimen tested by polymerase chain reaction (PCR) - `NPPCR` and
      blood (B) tested by culture (Cx) - `BCX`

  `Lambda`

  :   controls' subclass weights \\\nu_1, \nu_2, \ldots, \nu_K\\ a
      vector of `K` probabilities that sum to 1.

  `Eta`

  :   a matrix of dimension `length(cause_list)` by `K`; each row
      represents a disease class (among cases); the values in that row
      are subclass weights \\\eta_1, \eta_2, \ldots, \eta_K\\ for that
      disease class, so needs to sum to one. In Wu et al. 2016 (JRSS-C),
      the subclass weights are the same across disease classes across
      rows. But when simulating data, one can specify rows with distinct
      subclass weights - it is a matter whether we can recover these
      parameters (possible when some cases' true disease classes are
      observed)

  `PsiBS/PsiSS`

  :   False positive rates for Bronze-Standard data and for
      Silver-Standard data. For example, the rows of `PsiBS` correspond
      to the dimension of the particular slice of BrS measures, e.g.,
      `10` for 10 causative agents measured by NPPCR; the columns
      correspond to `K` subclasses; generically, the dimension is `J` by
      `K` `PsiSS` is supposed to be a vector of all zeros (perfect
      specificity in silver-standard measures).

  `ThetaBS/ThetaSS`

  :   True positive rates \\\Theta\\ for Bronze-Standard data and for
      Silver-Standard data. Dimension is `J` by `K` (can contain `NA` if
      the total number of causative agents measured by BrS or SS exceeds
      the measured causative agents in SS. For example, in PERCH study,
      nasopharyngeal polymerase chain reaction (NPPCR; bronze-standard)
      may target 30 distinct pathogens, but blood culture (BCX;
      silver-standard) may only target a subset of the 30, so we have to
      specify `NA` in `ThetaSS`for those pathogens not targeted by BCX).

  `Nu`

  :   the number of control subjects

  `Nd`

  :   the number of case subjects

- latent_samples:

  simulated latent status for all the subjects, for use in simulating
  BrS measurements.

## Value

a data frame with first column being case-control status (case at top)
and columns of bronze-standard measurements

## See also

Other internal simulation functions:
[`simulate_latent()`](https://zhenkewu.com/baker/reference/simulate_latent.md),
[`simulate_ss()`](https://zhenkewu.com/baker/reference/simulate_ss.md)
