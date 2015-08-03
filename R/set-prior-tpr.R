#' Set true positive rate (TPR) prior ranges for bronze-standard data
#'
#' Current prior assignment let bacteria NPPCR to have uniform range, viral NPPCR
#' to have .5-1.0 range. The PCP (a fungus) NPPCR TPR is also set to be .5-1.0; PCP
#' has no blood culture measurements. Also, not all the bacteria have blood culture
#' measurments. One question is whether to use informative NPPCR range or non-informative
#' NPPCR range (0-100%).
#'
#'
#' @param model_options See \code{\link{nplcm}} function.
#' @param data_nplcm See \code{\link{assign_model}} function.
#' 
#' @return Parameters for the TPR priors for BrS data. It is a list of two lists (alpha and beta).
#' Alpha or beta is of the same length as the number of BrS measurement slices; each element of the list
#' is a numeric vector for alpha or beta parameters to specify Beta prior for TPRs.
#'
#' @export

set_prior_tpr_BrS <- function(model_options,data_nplcm){

  parsed_model <- assign_model(model_options,data_nplcm)
  if (parsed_model$num_slice["MBS"] == 0){stop("== No BrS data! ==")}
    
  Mobs <- data_nplcm$Mobs
  Y    <- data_nplcm$Y
  X    <- data_nplcm$X
  
  Nd <- sum(Y)
  
  if (parsed_model$num_slice["MBS"] > 0){
    
    #
    # 1. BrS data:
    #
    # mapping template (by `make_template` function):
    patho_BrS_list <- lapply(Mobs$MBS,colnames)
    template_BrS_list <-
      lapply(patho_BrS_list,make_template,model_options$likelihood$cause_list) # key.
    
    prior_BrS <- model_options$prior$TPR_prior$BrS
    res_alpha <- list()
    res_beta <- list()
    
    if (prior_BrS$info == "non-informative"){
        for (s in seq_along(Mobs$MBS)){
            res_alpha[[s]] <- res_beta[[s]] <- rep(1,ncol(Mobs$MBS[[s]]))
        }
      return(list(alpha = res_alpha,beta = res_beta))
    }
    
    if (prior_BrS$info == "informative"){
      if (prior_BrS$input == "match_range"){# begin match range:
        for (s in seq_along(Mobs$MBS)){ #iterate over slices:
            tmp_ab <- matrix(NA,nrow=2,ncol(Mobs$MBS[[s]]))
            rownames(tmp_ab) <- c("alpha","beta")
            colnames(tmp_ab) <- colnames(Mobs$MBS[[s]])
            for (j in 1:ncol(Mobs$MBS[[s]]) ){ # iterate over dimensions:
              low_tmp <- prior_BrS$val[[s]]$low[j]
              up_tmp <- prior_BrS$val[[s]]$up[j]
              tmp <- beta_parms_from_quantiles(c(low_tmp,up_tmp),p=c(0.025,.975),plot=FALSE)
              tmp_ab [,j] <- c(tmp$a,tmp$b)
            }# end iterate over dimensions.
            res_alpha[[s]] <- tmp_ab[1,]
            res_beta[[s]]  <- tmp_ab[2,]
        }# end iterate over slices.
        return(list(alpha = res_alpha,beta = res_beta))
      } # end match range.
      
      if (prior_BrS$input == "direct_beta_param"){
          return(list(alpha = lapply(prior_BrS$val,"[[","alpha"),
                      beta  = lapply(prior_BrS$val,"[[","beta")))
        }
    } # end informative.
  }
}
  
# what1 <- set_prior_tpr_BrS(m_opt,data_nplcm)
# what2 <- set_prior_tpr_BrS(m_opt2,data_nplcm)







#' Set true positive rate (TPR) prior ranges for silver-standard data
#'
#' Current prior assignment let bacteria NPPCR to have uniform range, viral NPPCR
#' to have .5-1.0 range. The PCP (a fungus) NPPCR TPR is also set to be .5-1.0; PCP
#' has no blood culture measurements. Also, not all the bacteria have blood culture
#' measurments. One question is whether to use informative NPPCR range or non-informative
#' NPPCR range (0-100%).
#'
#'
#' @param model_options See \code{\link{nplcm}} function.
#' @param data_nplcm See \code{\link{assign_model}} function.
#' 
#' @return Parameters for the TPR priors for SS data. It is a list of two lists (alpha and beta).
#' Alpha or beta is of the same length as the number of BrS measurement slices; each element of the list
#' is a numeric vector for alpha or beta parameters to specify Beta prior for TPRs.
#'
#' @export
set_prior_tpr_SS <- function(model_options,data_nplcm){
  
  parsed_model <- assign_model(model_options,data_nplcm)
  if (parsed_model$num_slice["MSS"] == 0){stop("== No SS data! ==")}
  
  Mobs <- data_nplcm$Mobs
  Y    <- data_nplcm$Y
  X    <- data_nplcm$X
  
  Nd <- sum(Y)
  
  if (parsed_model$num_slice["MSS"] > 0){
    
    #
    # 2. SS data:
    #
    # mapping template (by `make_template` function):
    patho_SS_list <- lapply(Mobs$MSS,colnames)
    template_SS_list <-
      lapply(patho_SS_list,make_template,model_options$likelihood$cause_list) # key.
    
    prior_SS <- model_options$prior$TPR_prior$SS
    
    if (!is.null(prior_SS$grp)){GSS_TPR <- length(unique(prior_SS$grp))}
    if (is.null(prior_SS$grp)){GSS_TPR <- 1}
      res_all_grp <- list()
      
      for (g in 1:GSS_TPR){ # iterate over TPR split groups:
        res_alpha <- list()
        res_beta <- list()
        
        if (prior_SS$info == "non-informative"){
          for (s in seq_along(Mobs$MSS)){
            res_alpha[[s]] <- res_beta[[s]] <- rep(1,ncol(Mobs$MSS[[s]]))
          }
          res_tmp <- (list(alpha = res_alpha,beta = res_beta))
        }
        
        if (prior_SS$info == "informative"){
          if (GSS_TPR == 1){ # ad hoc solution to non-list problem when GSS_TPR==1.
            curr_val <- list(prior_SS$val)[g]
          } else{
            curr_val <- prior_SS$val[g]
          }
          if (prior_SS$input == "match_range"){# begin match range:
            for (s in seq_along(Mobs$MSS)){ #iterate over slices:
              tmp_ab <- matrix(NA,nrow=2,ncol(Mobs$MSS[[s]]))
              rownames(tmp_ab) <- c("alpha","beta")
              colnames(tmp_ab) <- colnames(Mobs$MSS[[s]])
              for (j in 1:ncol(Mobs$MSS[[s]]) ){ # iterate over dimensions:
                low_tmp <- curr_val[[s]]$low[j]
                up_tmp <- curr_val[[s]]$up[j]
                tmp <- beta_parms_from_quantiles(c(low_tmp,up_tmp),p=c(0.025,.975),plot=FALSE)
                tmp_ab [,j] <- c(tmp$a,tmp$b)
              }# end iterate over dimensions.
              res_alpha[[s]] <- tmp_ab[1,]
              res_beta[[s]]  <- tmp_ab[2,]
            }# end iterate over slices.
            res_tmp <- (list(alpha = res_alpha,beta = res_beta))
          } # end match range.
          
          if (prior_SS$input == "direct_beta_param"){
            res_tmp <- (list(alpha = lapply(curr_val,"[[","alpha"),
                             beta  = lapply(curr_val,"[[","beta")))
          }
        } # end informative.
        res_all_grp[[g]] <- res_tmp
      } # end iteration over group split.
    
    return(res_all_grp)
  }
}
# 
# what3 <- set_prior_tpr_SS(m_opt3,data_nplcm)
# what4 <- set_prior_tpr_SS(m_opt4,data_nplcm)



