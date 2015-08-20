#' Set true positive rate (TPR) prior ranges for bronze-standard data
#' (for conditional independence models - currently also used for conditional dependence model).
#' 
#' @param slice the index of BrS measurement under consideration
#' @param model_options See \code{\link{nplcm}} function.
#' @param data_nplcm See \code{\link{assign_model}} function.
#' 
#' @return Parameters for the TPR priors for BrS data. It is a list of two lists (alpha and beta).
#' Alpha and beta is of the same length, the number of BrS measurement slices. 
#' Each element of the alpha (beta) list is a numeric vector for alpha (beta) 
#' parameters to specify Beta prior for TPRs.
#'
#' @export

set_prior_tpr_BrS_NoNest <- function(slice,model_options,data_nplcm){

  parsed_model <- assign_model(model_options,data_nplcm)
  if (parsed_model$num_slice["MBS"] == 0){stop("== No BrS data! ==")}
    
  Mobs <- data_nplcm$Mobs
  Y    <- data_nplcm$Y
  X    <- data_nplcm$X
  
  Nd <- sum(Y)
  
  if (parsed_model$num_slice["MBS"] > 0){
    
    # mapping template (by `make_template` function):
    patho_BrS_list <- lapply(Mobs$MBS,colnames)
    template_BrS_list <-
      lapply(patho_BrS_list,make_template,model_options$likelihood$cause_list) # key.
    
    prior_BrS <- model_options$prior$TPR_prior$BrS
    
    res_all_slice <- vector("list",length=parsed_model$num_slice["MBS"])
    names(res_all_slice) <- names(Mobs$MBS)
    if (prior_BrS$info == "non-informative"){
        for (s in seq_along(Mobs$MBS)){
            res_all_slice[[s]] <- list(alpha = rep(1,ncol(Mobs$MBS[[s]])),
                                       beta = rep(1,ncol(Mobs$MBS[[s]])))
        }
      return(res_all_slice[slice])
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
              tmp_ab[,j] <- c(tmp$a,tmp$b)
            }# end iterate over dimensions.
            alpha_vec <- tmp_ab[1,]; names(alpha_vec) <- colnames(Mobs$MBS[[s]])
            beta_vec  <- tmp_ab[2,]; names(beta_vec) <- colnames(Mobs$MBS[[s]])
            res_all_slice[[s]] <- list(alpha = alpha_vec,beta = beta_vec)
        }# end iterate over slices.
        return(res_all_slice[slice])
      }# end match range.
      
      # begin direct parameters for Beta:
      if (prior_BrS$input == "direct_beta_param") {
        for (s in seq_along(Mobs$MBS)){
          curr_val <- prior_BrS$val[[s]]
          tmp_alpha <- curr_val$alpha; names(tmp_alpha) <- colnames(Mobs$MBS[[s]])
          tmp_beta  <- curr_val$beta; names(tmp_beta) <- colnames(Mobs$MBS[[s]])
          res_all_slice[[s]] <- list(alpha =  tmp_alpha, beta = tmp_beta)
        }
        return(res_all_slice[slice])
      }
    } # end informative.
  }
}
  
#' Set true positive rate (TPR) prior ranges for silver-standard data. 
#'
#' @param model_options See \code{\link{nplcm}} function.
#' @param data_nplcm See \code{\link{assign_model}} function.
#' 
#' @return Parameters for the TPR priors for BrS data. It is a list of two lists (alpha and beta).
#' Alpha and beta is of the same length, the number of BrS measurement slices. 
#' Each element of the alpha (beta) list is a numeric vector for alpha (beta) 
#' parameters to specify Beta prior for TPRs.
#'
#' @export
#' 
#' 
set_prior_tpr_SS <- function(model_options,data_nplcm){
  
  parsed_model <- assign_model(model_options,data_nplcm)
  if (parsed_model$num_slice["MSS"] == 0){stop("== No SS data! ==")}
  
  Mobs <- data_nplcm$Mobs
  Y    <- data_nplcm$Y
  X    <- data_nplcm$X
  
  Nd <- sum(Y)
  
  # mapping template (by `make_template` function):
  patho_SS_list <- lapply(Mobs$MSS,colnames)
  template_SS_list <-
    lapply(patho_SS_list,make_template,model_options$likelihood$cause_list) # key.
  
  prior_SS <- model_options$prior$TPR_prior$SS
  
  GSS_TPR <- 1
  if (parsed_model$SS_grp){GSS_TPR <- length(unique(prior_SS$grp))} # if there is grouping variable,
  
  res_all_slice <- vector("list",length = parsed_model$num_slice["MSS"])
  names(res_all_slice) <- names(Mobs$MSS)
  for (s in seq_along(Mobs$MSS)){ #iterate over slices:
    res_all_grp <- vector("list",length= GSS_TPR)
    if (prior_SS$info == "non-informative"){ # <------ not correct; also need to stratify by group.
      stop("== Not implemented. Please contact maintainer for an update, or 
               comment on the github page. Thanks. ==")
    }
    
    if (prior_SS$info == "informative"){
      curr_val <- prior_SS$val[[s]]
      # begin match range:
      if (prior_SS$input == "match_range") {
        for (g in 1:GSS_TPR){
          # iterate over TPR split groups:
          tmp_ab <- matrix(NA,nrow = 2,ncol(Mobs$MSS[[s]]))
          rownames(tmp_ab) <- c("alpha","beta")
          colnames(tmp_ab) <- colnames(Mobs$MSS[[s]])
          for (j in 1:ncol(Mobs$MSS[[s]])) {
            # iterate over dimensions:
            low_tmp <- curr_val[[g]]$low[j]
            up_tmp <- curr_val[[g]]$up[j]
            tmp <-
              beta_parms_from_quantiles(c(low_tmp,up_tmp),p = c(0.025,.975),plot = FALSE)
            tmp_ab [,j] <- c(tmp$a,tmp$b)
          }# end iterate over dimensions.
          alpha_vec <- tmp_ab[1,]; names(alpha_vec) <- colnames(Mobs$MSS[[s]])
          beta_vec  <- tmp_ab[2,]; names(beta_vec) <- colnames(Mobs$MSS[[s]])
          res_all_grp[[g]] <- list(alpha = alpha_vec,beta = beta_vec)
        }# end iteration over group split.
      } # end match range.
      
      # begin direct parameters for Beta:
      if (prior_SS$input == "direct_beta_param") {
        for (g in 1:GSS_TPR){
          tmp_alpha <- curr_val[[g]]$alpha; names(tmp_alpha) <- colnames(Mobs$MSS[[s]])
          tmp_beta <- curr_val[[g]]$beta; names(tmp_beta) <- colnames(Mobs$MSS[[s]])
          res_all_grp[[g]] <- list(alpha =  tmp_alpha, beta = tmp_beta)# <--- note here we lost column names for each measurements.
        }
      }
      
      res_all_slice[[s]] <- res_all_grp
    } # end informative.
    return(res_all_slice)
  }# end iterate over slices.
}


