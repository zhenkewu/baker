#' Calculate marginal log odds ratios
#' 
#' @inheritParams simulate_nplcm
#' 
#' @example /inst/example/plot_associations.R
#'
#' @return A figure showing pairwise odds ratios for cases (upper right, solid lines) 
#' and controls (lower left, broken lines) as the first subclass weight increases 
#' from 0 to 1. Pairwise independence is represented by the dotted horizontal lines 
#' for reference.
#' 
#' @export
#'
compute_logOR_single_cause <- function(set_parameter){
  eti <- set_parameter$etiology
  theta <- set_parameter$ThetaBS
  psi <- set_parameter$PsiBS
  lambda <- set_parameter$Lambda
  eta <- set_parameter$Eta
  K <- max(sum(lambda!=0), sum(eta[1,]!=0))
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

#' Calculate marginal log odds ratios (with other causes)
#' 
#' @inheritParams simulate_nplcm
#' 
#' @return a matrix
#' 
#' @export
#'
compute_logOR_single_and_other_cause <- function(set_parameter){
  eti <- set_parameter$etiology
  theta <- set_parameter$ThetaBS
  psi <- set_parameter$PsiBS
  lambda <- set_parameter$Lambda
  eta <- set_parameter$Eta
  K <- max(sum(lambda!=0), sum(eta[1,]!=0))
  k_ind_case <- 1:K
  k_ind_ctrl <- 1:K
  if (K==1){
    k_ind_case <- which(eta[1,]!=0)  
    k_ind_ctrl <- which(lambda!=0)  
  }
  L <- length(eti)
  J <- L-1
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
      for (c in L){
        for (k in k_ind_case){
          A_subclass_temp <- A_subclass_temp + eta[j,k]*psi[j,k]*psi[l,k]
          B_subclass_temp <- B_subclass_temp + eta[j,k]*(1-psi[j,k])*psi[l,k]
          C_subclass_temp <- C_subclass_temp + eta[j,k]*(1-psi[j,k])*(1-psi[l,k])
          D_subclass_temp <- D_subclass_temp + eta[j,k]*psi[j,k]*(1-psi[l,k])
        }
        A <- A + eti[L]*A_subclass_temp  
        B <- B + eti[L]*B_subclass_temp  
        C <- C + eti[L]*C_subclass_temp  
        D <- D + eti[L]*D_subclass_temp  
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
