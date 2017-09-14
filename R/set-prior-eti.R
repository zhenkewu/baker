#' specify overall uniform (symmetric Dirichlet distribution) for etiology prior
#' 
#' @param alpha any positive number, usually 1.
#' 
#' @param cause_list a list of latent status
#' 
#' @return a vector of length \code{length(cause_list)}
#' @family prior specification functions
#' @export
#'
overall_uniform <- function(alpha,cause_list) {
  if (length(alpha) > 1) {
    stop("== [baker] overall_uniform is lazy and only needs a single hyperprior for etiology; please use a positive value for 'alpha'! 
         If a non-uniform etiology prior is needed, just specify a vector of 'alpha' of equal length to 'cause_list'! ==")
  }
  rep(alpha, length(cause_list))
}


# @param n_stratum the number of strata; Default to \code{NULL} (just a single uniform
# prior); for any positive integer, the function will generate a matrix of \code{n_stratum}
#rows and \code{length(cause_list)} columns, with each row representing an etiology prior.
# 
# ,n_stratum=NULL