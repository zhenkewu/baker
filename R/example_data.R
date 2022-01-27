#' Hypothetical pathogens and their categories (virus or bacteria)
#'
#' This is used in simulations where the pathogen names are from the alphabet,
#' and we hope to plot etiologies grouped by virus or bacteria
#'
#' @format A matrix of two columns
#' \describe{
#'   \item{pathogen}{names of the hypothetical pathogens, A-Z}
#'   \item{pathogen_type}{category of the hypothetical pathogens, `B` for
#'   bacterium, `V` for virus, which are randomly assigned.}
#' }
#' @usage data("pathogen_category_simulation")
#' 
#' 
#' @return No returned value; just loading data into the working space.
"pathogen_category_simulation"


#' pathogens and their categories in PERCH study (virus or bacteria)
#'
#' 231 rows indicating bacteria, virus, fungi, or other categories.
#'
#' @format A matrix of two columns
#' \describe{
#'   \item{pathogen}{names of the pathogens}
#'   \item{pathogen_type}{category of the pathogens, `B` for
#'   bacterium, `V` for virus, `F` for fungus, `O` for "not categorized"}
#' }
#' 
#' @return No returned value; just loading data into the working space.
#' 
#' @usage data("pathogen_category_perch")
"pathogen_category_perch"

## #' Simulated dataset that is structured in the format necessary for an [nplcm()] with regression 
## #'
## #' Data set for illustrating regression functionalities
## #' 
## #' @format A list containing three items 
## #' \describe{
## #'   \item{Mobs}{BrS level measurements: N = 1,200 (half cases and half controls);
## #'   one slice of BrS measurements (6 dimensional, A-F); one slice of SS measurements (2 dimensional,
## #'   A and B)}
## #'   \item{Y}{case-control status}
## #'   \item{X}{matrix of covariates (N by 4); columns: SITE (1 and 2, each with 600 subjects), 
## #'   DATE (index from 1:300), std_date (standardized DATE), ENRLDATE (actual dates)}
## #' }
## #' 
## #' @usage data("data_nplcm_reg")
## "data_nplcm_reg"


#' Simulated dataset that is structured in the format necessary for an [nplcm()] with regression 
#'
#' Data set for illustrating regression functionalities
#' 
#' @format A list containing three items 
#' \describe{
#'   \item{Mobs}{BrS level measurements: N = 1,200 (half cases and half controls);
#'   one slice of BrS measurements (6 dimensional, A-F); one slice of SS measurements (2 dimensional,
#'   A and B)}
#'   \item{Y}{case-control status}
#'   \item{X}{matrix of covariates (N by 4); columns: SITE (1 and 2, each with 600 subjects), 
#'   DATE (index from 1:300), std_date (standardized DATE), ENRLDATE (actual dates)}
#' }
#' 
#' @return No returned value; just loading data into the working space.
#' @usage data("data_nplcm_reg_nest")
"data_nplcm_reg_nest"

#' Simulated dataset that is structured in the format necessary for an [nplcm()] with regression 
#'
#' Data set for illustrating regression functionalities
#' 
#' @format A list containing three items 
#' \describe{
#'   \item{Mobs}{BrS level measurements: N = 600 (half cases and half controls);
#'   one slice of BrS measurements (6 dimensional, A-F); one slice of SS measurements (2 dimensional,
#'   A and B)}
#'   \item{Y}{case-control status}
#' }
#' 
#' @return No returned value; just loading data into the working space.
#' 
#' @usage data("data_nplcm_noreg")
"data_nplcm_noreg"



# NB: a few more data sets: data_nplcm_no_reg (no regression setting); this can be used in a lot of functions
# for visualization. Perhaps also need to add posterior samples as data set just as a way to 
# demonstrate the plotting functionality and unit test.
