#' Calculate marginal log odds ratios
#' 
#' This only works for single-agent causes
#' 
#' @inheritParams simulate_nplcm
#' 
#' @example /inst/example/plot_associations.R
#'
#' @return a matrix of log odds ratio.
#' See the example for a figure showing pairwise odds ratios for cases (upper right, solid lines) 
#' and controls (lower left, broken lines) as the first subclass weight increases 
#' from 0 to 1. Pairwise independence is represented by the dotted horizontal lines 
#' for reference.
#' 
#' @export
#'
compute_logOR_single_cause <- function(set_parameter){
  eti  <- set_parameter$etiology
  theta <- set_parameter$ThetaBS
  psi   <- set_parameter$PsiBS
  lambda <- set_parameter$Lambda
  eta <- set_parameter$Eta
  K   <- max(sum(lambda!=0), sum(eta[1,]!=0))
  k_ind_case <- 1:K
  k_ind_ctrl <- 1:K
  if (K==1){
    k_ind_case <- which(eta[1,]!=0)  
    k_ind_ctrl <- which(lambda!=0)  
  }
  J <- length(eti)
  MAT <- matrix(NA,nrow=J, ncol = J)
  for (j in 1:(J-1)){
    for (l in (j+1):J){
      A <- 0; B <- 0; C <- 0; D <- 0
      for (c in 1:J){
        ind_j <- 0+(c==j)
        ind_l <- 0+(c==l)
        A_subclass_temp <- 0
        B_subclass_temp <- 0
        C_subclass_temp <- 0
        D_subclass_temp <- 0
        for (k in k_ind_case){
          A_subclass_temp <- A_subclass_temp + eta[j,k]*theta[j,k]^ind_j*psi[j,k]^(1-ind_j)*(theta[l,k])^(ind_l)*(psi[l,k])^(1-ind_l)
          B_subclass_temp <- B_subclass_temp + eta[j,k]*(1-theta[j,k])^ind_j*(1-psi[j,k])^(1-ind_j)*(theta[l,k])^(ind_l)*(psi[l,k])^(1-ind_l)
          C_subclass_temp <- C_subclass_temp + eta[j,k]*(1-theta[j,k])^ind_j*(1-psi[j,k])^(1-ind_j)*(1-theta[l,k])^(ind_l)*(1-psi[l,k])^(1-ind_l)
          D_subclass_temp <- D_subclass_temp + eta[j,k]*theta[j,k]^ind_j*psi[j,k]^(1-ind_j)*(1-theta[l,k])^(ind_l)*(1-psi[l,k])^(1-ind_l)
        }
        A <- A + eti[c]*A_subclass_temp  
        B <- B + eti[c]*B_subclass_temp  
        C <- C + eti[c]*C_subclass_temp  
        D <- D + eti[c]*D_subclass_temp  
      }
      
      MAT[j,l] <- log(A)-log(B)+log(C)-log(D)
    }
  }
  
  for (l in 1:(J-1)){
    for (j in (l+1):J){
      A <- 0
      B <- 0
      C <- 0
      D <- 0
      for (k in k_ind_ctrl){
        A <- A + lambda[k]*psi[j,k]*psi[l,k]
        B <- B + lambda[k]*(1-psi[j,k])*(1-psi[l,k])
        C <- C + lambda[k]*psi[j,k]*(1-psi[l,k])
        D <- D + lambda[k]*(1-psi[j,k])*psi[l,k]
      }
      MAT[j,l] <- log(A)+log(B)-log(C)-log(D)
    }  
  }
  MAT
}

## Calculate marginal log odds ratios (with other causes)
## 
## This only works for single-agent causes
## 
## @inheritParams simulate_nplcm
## 
## @return a matrix of log odds ratio
## 
## @export
##
#compute_logOR_single_and_other_cause <- function(set_parameter){
#  eti <- set_parameter$etiology
#  theta <- set_parameter$ThetaBS
#  psi <- set_parameter$PsiBS
#  lambda <- set_parameter$Lambda
#  eta <- set_parameter$Eta
#  K <- max(sum(lambda!=0), sum(eta[1,]!=0))
#  k_ind_case <- 1:K
#  k_ind_ctrl <- 1:K
#  if (K==1){
#    k_ind_case <- which(eta[1,]!=0)  
#    k_ind_ctrl <- which(lambda!=0)  
#  }
#  L <- length(eti)
#  J <- L-1
#  MAT <- matrix(NA,nrow=J, ncol = J)
#  for (j in 1:(J-1)){
#    for (l in (j+1):J){
#      A <- 0; B <- 0; C <- 0; D <- 0
#      for (c in 1:J){
#        ind_j <- 0+(c==j)
#        ind_l <- 0+(c==l)
#        A_subclass_temp <- 0
#        B_subclass_temp <- 0
#        C_subclass_temp <- 0
#        D_subclass_temp <- 0
#        for (k in k_ind_case){
#          A_subclass_temp <- A_subclass_temp + eta[j,k]*theta[j,k]^ind_j*psi[j,k]^(1-ind_j)*(theta[l,k])^(ind_l)*(psi[l,k])^(1-ind_l)
#          B_subclass_temp <- B_subclass_temp + eta[j,k]*(1-theta[j,k])^ind_j*(1-psi[j,k])^(1-ind_j)*(theta[l,k])^(ind_l)*(psi[l,k])^(1-ind_l)
#          C_subclass_temp <- C_subclass_temp + eta[j,k]*(1-theta[j,k])^ind_j*(1-psi[j,k])^(1-ind_j)*(1-theta[l,k])^(ind_l)*(1-psi[l,k])^(1-ind_l)
#          D_subclass_temp <- D_subclass_temp + eta[j,k]*theta[j,k]^ind_j*psi[j,k]^(1-ind_j)*(1-theta[l,k])^(ind_l)*(1-psi[l,k])^(1-ind_l)
#        }
#        A <- A + eti[c]*A_subclass_temp  
#        B <- B + eti[c]*B_subclass_temp  
#        C <- C + eti[c]*C_subclass_temp  
#        D <- D + eti[c]*D_subclass_temp  
#      }
#      for (c in L){
#        for (k in k_ind_case){
#          A_subclass_temp <- A_subclass_temp + eta[j,k]*psi[j,k]*psi[l,k]
#          B_subclass_temp <- B_subclass_temp + eta[j,k]*(1-psi[j,k])*psi[l,k]
#          C_subclass_temp <- C_subclass_temp + eta[j,k]*(1-psi[j,k])*(1-psi[l,k])
#          D_subclass_temp <- D_subclass_temp + eta[j,k]*psi[j,k]*(1-psi[l,k])
#        }
#        A <- A + eti[L]*A_subclass_temp  
#        B <- B + eti[L]*B_subclass_temp  
#        C <- C + eti[L]*C_subclass_temp  
#        D <- D + eti[L]*D_subclass_temp  
#      }
#      
#      MAT[j,l] <- log(A)-log(B)+log(C)-log(D)
#    }
#  }
#  
#  for (l in 1:(J-1)){
#    for (j in (l+1):J){
#      A <- 0
#      B <- 0
#      C <- 0
#      D <- 0
#      for (k in k_ind_ctrl){
#        A <- A + lambda[k]*psi[j,k]*psi[l,k]
#        B <- B + lambda[k]*(1-psi[j,k])*(1-psi[l,k])
#        C <- C + lambda[k]*psi[j,k]*(1-psi[l,k])
#        D <- D + lambda[k]*(1-psi[j,k])*psi[l,k]
#      }
#      MAT[j,l] <- log(A)+log(B)-log(C)-log(D)
#    }  
#  }
#  MAT
#}
#


#' compute positive rates for nested model with subclass mixing weights that are the same across
#' `Jcause` classes for each person (people may have different weights.)
#' 
#' The array version of this function ([compute_marg_PR_nested_reg_array]) is used in [plot_etiology_regression]
#' 
#' @param ThetaBS True positive rates for `JBrS` measures (rows) among `K` subclasses (columns)
#' @param PsiBS False positive rates; dimension same as above
#' @param pEti_mat a matrix of etiology pies for `N` subjects (rows) and `Jcause` causes (columns)
#' rows sum to ones.
#' @param subwt_mat a matrix of subclass weights for cases and controls. `N` by `K`. Rows sum
#' to ones.
#' @param case a N-vector of `1`s (cases) and `0`s (controls)
#' @param template a binary matrix with `Jcause+1` rows (`Jcause` classes of cases and `1` class of controls)
#' and `JBrS` columns for the Bronze-standard measurement (say, pick one type/slice).
#' The ones in each row indicate the measurements that will show up more frequently in cases given the cause.
#' @return a matrix of values between `0` and `1` (need not to have row sums of ones); 
#' of dimension (number of subjects, dimension of the bronze-standard measurement slice).
#' 
compute_marg_PR_nested_reg <-  function(ThetaBS,PsiBS,pEti_mat,subwt_mat,case,template){
  # ThetaBS = ThetaBS_samp[,,1];PsiBS = PsiBS_samp[,,1];
  # pEti_mat = pEti_samp[,,1];subwt_mat = subwt_samp[,,1];
  # case = data_nplcm$Y;template = templateBS
  Jcause <- ncol(pEti_mat)
  N_all  <- nrow(pEti_mat)
  JBrS   <- ncol(template)
  pEti_mat <- sweep(pEti_mat,1,case,"*") # set etiology for controls to all zeros.
  if (ncol(ThetaBS)!=ncol(PsiBS) | nrow(ThetaBS)!=nrow(PsiBS) ){
    stop("==[baker] 'ThetaBS' and 'PsiBS' must have identical # of rows and # of columns. ==")}
  K      <- ncol(ThetaBS)
  
  FPR_mat <- subwt_mat%*%t(PsiBS)
  TPR_mat <- sweep(subwt_mat,1,case,"*")%*%t(ThetaBS) # set control TPRs to be zero.
  res <- matrix(0,nrow=N_all,ncol=JBrS)
  # operate on controls:
  res[case==0,] <- FPR_mat[case==0,]
  # operate on cases:
  for (j in 1:Jcause){
    PR_mat_given_cause <- sweep(TPR_mat,2,template[j,],"*")+sweep(FPR_mat,2,1-template[j,],"*")
    # cases are either TPR by their subwt or FPR by their subwts. Controls are FPR by their subwts (TPR for them are zero).
    res <- res + sweep(PR_mat_given_cause,1,pEti_mat[,j],"*")
  }
  res
}


#' compute positive rates for nested model with subclass mixing weights that are the same across
#' `Jcause` classes for each person (people may have different weights.)
#' 
#' This is an array-version of [compute_marg_PR_nested_reg]. This is used in [plot_etiology_regression]
#' 
#' @param ThetaBS_array An array of: True positive rates for JBrS measures (rows) among K subclasses (columns)
#' @param PsiBS_array An array of: False positive rates; dimension same as above
#' @param pEti_mat_array An array of: a matrix of etiology pies for N subjects (rows) and Jcause causes (columns)
#' rows sum to ones.
#' @param subwt_mat_array An array of: a matrix of subclass weights for cases and controls. N by K. Rows sum
#' to ones.
#' @param case a N-vector of 1s (cases) and 0s (controls)
#' @param template a binary matrix with Jcause+1 rows (Jcause classes of cases and 1 class of controls)
#' and JBrS columns for the Bronze-standard measurement (say, pick one type/slice).
#' The ones in each row indicate the measurements that will show up more frequently in cases given the cause.
#' @return An array of: a matrix of values between 0 and 1 (need not to have row sums of ones); of dimension (number of subjects, dimension of the bronze-standard
#' measurement slice).
#' 
compute_marg_PR_nested_reg_array <- function(ThetaBS_array,PsiBS_array,
                                             pEti_mat_array,subwt_mat_array,case,template){
  res <- array(NA,c(dim(pEti_mat_array)[1],dim(ThetaBS_array)[1],dim(pEti_mat_array)[3]))  
  for (s in 1:(dim(ThetaBS_array)[3])){ 
    res[,,s] <- compute_marg_PR_nested_reg(
      ThetaBS_array[,,s],PsiBS_array[,,s],
      pEti_mat_array[,,s],subwt_mat_array[,,s],case,template)
  }
  res
}

