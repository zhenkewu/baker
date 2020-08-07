#' Hypothetical pathogens and their categories (virus or bacteria)
#'
#' This is used in simulations where the pathogen names are from the alphabet,
#' and we hope to plot etiologies grouped by virus or bacteria
#'
#' @format A matrix of two columns
#' \describe{
#'   \item{pathogen}{names of the hypothetical pathogens, A-Z}
#'   \item{pathogen_type}{category of the hypothetical pathogens, \code{B} for
#'   bacterium, \code{V} for virus, which are randomly assigned.}
#' }
"pathogen_category_simulation"
