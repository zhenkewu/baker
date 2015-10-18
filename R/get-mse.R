#' Obtain MSE, bias-squared and variance of posterior from nplcm fitted folders.
#'
#' The result is equivalent to Euclidean-type calculation after the
#' compositional vector (e.g., etiologic fraction) is centered-log-ratio (CLRB) tranformed.
#' For simulation only.
#'
#' @param DIR_NPLCM File path where Bayesian results are stored
#' @param truth True etiologic fraction vector (must sum to 1) 
#' used to generate data
#' @importFrom robCompositions cenLR
#' @return a vector of (Integrated Squared Aitchison Distance (ISAD), bias-squared, variance, truth)
#' @export
#'

get_metric <- function(DIR_NPLCM,truth){
  use_jags <- is_jags_folder(DIR_NPLCM)
  my.mse <- function(A,B){
    if (nrow(A)!=nrow(B)){
      stop("not equal number of rows!")
    }else{
      aa=cenLR(A)$x.clr
      bb=cenLR(B)$x.clr
      res <- rowSums((aa-bb)^2)
      mean(res)
    }
  }
  
  my.bias.sq <- function(A,B){
    if (nrow(A)!=nrow(B)){
      stop("not equal number of rows!")
    }else{
      aa=cenLR(A)$x.clr
      bb=cenLR(B)$x.clr
      sum((colMeans(aa-bb))^2)
    }
  }
  
  #A = as.data.frame(true_pEti)
  my.var <- function(A){
    res <- cenLR(A)$x.clr
    sum(diag(cov(res)))
  }
  
  #read in data from the folder directory:
  
  out <- nplcm_read_folder(DIR_NPLCM)
  model_options <- out$model_options
  bugs.dat      <- out$bugs.dat
  res_nplcm     <- out$res_nplcm
  #some data preparation:
  Nd <- bugs.dat$Nd
  Nu <- bugs.dat$Nu
  
  Y <-  out$Y
  
  cause_list     <- model_options$likelihood$cause_list
  Jcause         <- length(grep("^pEti\\[",colnames(res_nplcm)))
  if (Jcause != length(cause_list)){
    stop("=='model_options' and actual posterior results 'pEti' have different number of causes!==")
  }
  
  # extract and process some data and posterior samples:
  SubVarName <- rep(NA,Jcause)
  for (j in 1:Jcause){
    SubVarName[j] = paste("pEti","[",j,"]",sep="")
  }
  #get etiology fraction MCMC samples:
  pEti_mat   <- as.matrix(res_nplcm[,SubVarName])
  true_pEti  <- as.data.frame(t(replicate(nrow(pEti_mat),truth)))
  
  # Aitchison distances:
  # aDist(true_pEti,pEti_mat)
  mse     <- my.mse(true_pEti,pEti_mat)
  bias.sq <- my.bias.sq(true_pEti,pEti_mat)
  vv      <- my.var(as.data.frame(pEti_mat))
  
  res <- make_list(mse, bias.sq,vv,truth)
  res
}





