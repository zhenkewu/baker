#' Obtain Integrated Squared Aitchison Distance, Squared Bias and Variance (both on 
#' Central Log-Ratio transformed scale) that measure the discrepancy of a posterior
#' distribution of pie and a true pie. 
#'
#' The result is equivalent to Euclidean-type calculation after the
#' compositional vector (e.g., etiologic fraction) is centered-log-ratio (CLRB) transformed.
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
    sum(diag(stats::cov(res)))
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

#' Obtain direct bias that measure the discrepancy of a posterior
#' distribution of pie and a true pie. 
#'
#' @param DIR_list The list of  where Bayesian results are stored
#' @param truth True etiologic fraction vector (must sum to 1)  used to generate data;
#' Default is \code{NULL}. If a vector is supplied, then only the frist path in \code{DIR_LIST}
#' is used.
#' @param silent Default is FALSE. To suppress printing messages, set to TRUE.
#' 
#' @return a list of length two. \code{diff} is the direct differences; 
#' \code{prb} is the percent relative bias.
#' @export
#'
get_direct_bias <- function(DIR_list,truth=NULL,silent=FALSE){

  symdiff <- function( x, y) { setdiff( union(x, y), intersect(x, y))}
  
  if (is.null(truth)){
    # read from folders:
    out_list <- vector("list",length(DIR_list))
    base_nm  <- lapply(DIR_list,basename)
    names(out_list) <- base_nm
    
    pEti_samp_list  <- list()
    for (i in seq_along(DIR_list)){
      out_list[[i]]        <- nplcm_read_folder(DIR_list[[i]])
      pEti_samp_list[[i]]  <- get_pEti_samp(out_list[[i]]$res_nplcm,out_list[[i]]$model_options)
    }
    
    different_names <- symdiff(out_list[[1]]$model_options$likelihood$cause_list,out_list[[2]]$model_options$likelihood$cause_list)
    if (length(different_names) > 0 ){stop("==Two folders have different latent category names!==")}
    
    if(!silent){
      print("==The first folder in 'DIR_LIST' is used as reference for calculating relative bias. ==")
    }
    # plain difference:
    diff_mat <- pEti_samp_list[[2]]$pEti_mat - pEti_samp_list[[1]]$pEti_mat
    diff     <- colMeans(diff_mat)
    names(diff) <- pEti_samp_list[[1]]$latent_nm
    # percent relative bias:
    prb     <- colMeans(diff_mat)/colMeans(pEti_samp_list[[1]]$pEti_mat)
    names(prb) <- pEti_samp_list[[1]]$latent_nm
    
    res <- make_list(diff,prb)
    return(res)
  }
  
  if (!is.null(truth)){
    if(!silent){print("==Only the first folder in 'DIR_LIST' is used to compare with 'truth'!==")}
    # read from folders:
    out_list <- vector("list",length(DIR_list))
    base_nm  <- lapply(DIR_list,basename)
    names(out_list) <- base_nm
    
    pEti_samp_list  <- list()
    for (i in 1){
      out_list[[i]]        <- nplcm_read_folder(DIR_list[[i]])
      pEti_samp_list[[i]]  <- get_pEti_samp(out_list[[i]]$res_nplcm,out_list[[i]]$model_options)
    }
    
    if (length(out_list[[1]]$model_options$likelihood$cause_list) != length(truth)){
      stop("==Results and truth have different number of latent categories! ==")
    }
    
    # plain difference:
    diff_mat <- pEti_samp_list[[1]]$pEti_mat - t(replicate(nrow(pEti_samp_list[[1]]$pEti_mat),truth))
    diff     <- colMeans(diff_mat)
    names(diff) <- pEti_samp_list[[1]]$latent_nm
    # percent relative bias:
    prb     <- colMeans(diff_mat)/truth
    names(prb) <- pEti_samp_list[[1]]$latent_nm
    
    res <- make_list(diff,prb)
    return(res)
  }
}

#' Obtain coverage status from a result folder
#'
#' @param DIR_NPLCM Path to where Bayesian results are stored
#' @param truth True etiologic fraction vector (must sum to 1)  used to generate data.
#' 
#' @examples 
#' \dontrun{
#' DIR_NPLCM <- "~/downloads/rep_1_kfit_2/"  
#' truth     <- c(0.5,0.2,0.15,0.1,0.05)
#' get_coverage(DIR_NPLCM,truth)
#' }
#' @return A logic vector of length as \code{truth}. 1 for covered; 0 for not.
#' @export
#'
get_coverage <- function(DIR_NPLCM,truth){
    # read from folders:
    out       <- nplcm_read_folder(DIR_NPLCM)
    pEti_samp <- get_pEti_samp(out$res_nplcm,out$model_options)
    
    if (length(out$model_options$likelihood$cause_list) != length(truth)){
      stop("==Results and truth have different number of latent categories! ==")
    }
    
    # plain difference:
    diff_mat <- pEti_samp$pEti_mat - t(replicate(nrow(pEti_samp$pEti_mat),truth))
    UL <- apply(diff_mat,2,stats::quantile,0.975)
    LL <- apply(diff_mat,2,stats::quantile,0.025)
    res <- (UL>=0 & LL <=0)
    return(res)
}
