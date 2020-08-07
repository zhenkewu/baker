#' Simulated dataset that is structured in the format necessary for an nplcm with regression 
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
"data_nplcm_reg"