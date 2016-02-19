#' specify overall uniform (symmetric Dirichlet distribution) for etiology prior
#' 
#' @param alpha any positive number, usually 1s.
#' 
#' @param cause_list a list of latent status
#' 
#' @return a vector of length \code{length(cause_list)}
#' @family prior specification functions
#' @export
#'
overall_uniform <- function(alpha,cause_list) {
  if (length(alpha) == 1) {
    res <- rep(alpha, length(cause_list)); return(res)
  }
  if (length(alpha) > 1 &&
      length(alpha) < length(cause_list)) {
    stop("== not enough hyperparameters for etiology! ==")
  }
  alpha
}


# can add: taxonomy and others.
