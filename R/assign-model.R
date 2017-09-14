#' Interpret the model specified by user
#'
#' \code{assign_model} translates options specified by a user (e.g., in 
#' \code{model_options}) into information that can be understood by \code{baker}.
#' 
#' @details \code{assign_model} will be modified to check if data are conformable
#' to specified model.
#' 
#' @param data_nplcm Data.
#' @param model_options See \code{\link{nplcm}} function.
#' @param silent Default is \code{TRUE} for no messages; \code{FALSE} otherwise.
#' @return A list of model specifications:
#' \itemize{
#'    \item \code{num_slice} A vector counting the No. of measurement slices for each
#'    level of measurement quality;
#'    \item \code{nested} Local dependence specification for modeling bronze-standard
#'    data. \code{TRUE} for nested models (conditional dependence); 
#'    \code{FALSE} for non-nested models (conditional independence). 
#'    One for each BrS slice.
#'    \item \code{regression}
#'        \itemize{
#'            \item \code{do_reg_Eti} \code{TRUE} for doing etiology regression.
#'            It means let the etiology fractions to change with covariates. 
#'            \code{FALSE} otherwise;
#'            \item \code{do_reg_FPR} \code{TRUE} for doing false positive rate 
#'            regression (for every slice of bronze-standard). It means the false
#'            positive rates, usually estimatable from controls, can vary with 
#'            covariates. \code{FALSE} otherwise.
#'            \item code{is_discrete_predictor} A list of names "Eti", and 
#'            the names for every slice of bronze-standard data. \code{TRUE}
#'            for having all discrete predictors, \code{FALSE} otherwise.
#'        }
#' }
#'
#' @family specification checking functions
#' @export
assign_model <- function(model_options,data_nplcm, silent=TRUE){
  # load options:
  likelihood       <- model_options$likelihood
  use_measurements <- model_options$use_measurements
  prior            <- model_options$prior
  
  # load data: 
  Mobs <- data_nplcm$Mobs
  Y    <- data_nplcm$Y
  X    <- data_nplcm$X
  
  nested       <- likelihood$k_subclass > 1
  
  # test the match between actual data and model_options:
  use_data_sources   <- c("MBS","MSS","MGS")[lookup_quality(use_measurements)]
  input_data_sources <-  names(Mobs)
  if (!all(use_data_sources%in%input_data_sources)){
    stop("==[baker] Please supply actual datasets as specified by 'use_measurements' in 'model_options'.\n==")
  }
  
  # get the length of each measurement quality:
  num_slice <- rep(0,3)
  names(num_slice) <- c("MBS","MSS","MGS")
  for (i in seq_along(use_data_sources)){
    num_slice[use_data_sources[i]] <- length(Mobs[[use_data_sources[i]]])
  }
  
  # specify regression for FPR: (only available for bronze-standard data. Silver-standard data automatically have FPR==0.)
  do_reg_FPR <- list() #  <---- a regression for each measurement slice?
  is_discrete_FPR_vec <- rep(NA,length(likelihood$FPR_formula))
  names(is_discrete_FPR_vec) <- names(likelihood$FPR_formula)
  for (i in seq_along(Mobs$MBS)){
    ind_tmp <-
      which(names(likelihood$FPR_formula) == names(Mobs$MBS)[i])
    form_tmp <- stats::as.formula(likelihood$FPR_formula[[ind_tmp]])
    if (!length(ind_tmp)) { # don't do regression if no regression formula is found:
      do_reg_FPR[[i]] <- FALSE
    } else{ # do regression if there is matched regression formula:
      do_reg_FPR[[i]] <-
        parse_nplcm_reg(form_tmp,data_nplcm,silent=silent)
    }
    
    is_discrete_FPR_vec[i] <- FALSE
    if (!is.null(X)){
      is_discrete_FPR_vec[i] <- (!is_intercept_only(form_tmp) & !stats::is.empty.model(form_tmp) & 
                                   is_discrete(X, form_tmp))
    }
  }
  names(do_reg_FPR) <- names(Mobs$MBS)
  
  #
  # specify regression for TPR: (every measurement slice has it.)
  #
  # do_reg_TPR <- list() #  <---- a regression for each measurement slice?
  # for (i in seq_along(Mobs$MBS)){
  #   ind_tmp <-
  #     which(names(likelihood$TPR_formula) == names(Mobs$MBS)[i])
  #   if (!length(ind_tmp)) { # don't do regression if no regression formula is found:
  #     do_reg_TPR[[i]] <- FALSE
  #   } else{ # do regression if there is matched regression formula:
  #     do_reg_TPR[[i]] <-
  #       parse_nplcm_reg(stats::as.formula(likelihood$TPR_formula[[ind_tmp]]),data_nplcm,silent=silent)
  #   }
  # }
  # 
  # names(do_reg_TPR) <- names(Mobs$MBS)
  # 
  # # if using silver-standard data:
  # if ("MSS"%in% use_data_sources){
  #   for (i in length(Mobs$MBS)+seq_along(Mobs$MSS)){
  #     ind_tmp <-
  #       which(names(likelihood$TPR_formula) == names(Mobs$MSS)[i])
  #     if (!length(ind_tmp)) { # don't do regression if no regression formula is found:
  #       do_reg_TPR[[i]] <- FALSE
  #     } else{ # do regression if there is matched regression formula:
  #       do_reg_TPR[[i]] <-
  #         parse_nplcm_reg(stats::as.formula(likelihood$TPR_formula[[ind_tmp]]),data_nplcm,silent=silent)
  #     }
  #   }
  #   names(do_reg_TPR) <- c(names(Mobs$MBS),names(Mobs$MSS))
  # }
  # specify regression for etiology:
  form_tmp   <- stats::as.formula(likelihood$Eti_formula)
  do_reg_Eti <- parse_nplcm_reg(form_tmp,data_nplcm,silent=silent)
  
  is_discrete_Eti <- FALSE
  if (!is.null(X)){
    is_discrete_Eti <- (!is_intercept_only(form_tmp) & !stats::is.empty.model(form_tmp) & 
                          is_discrete(data.frame(X,Y)[Y==1,,drop=FALSE], form_tmp))
  }
  
  is_discrete_predictor <- list(is_discrete_Eti, is_discrete_FPR_vec)
  names(is_discrete_predictor)[1] <- "Eti"
  names(is_discrete_predictor)[2] <- "FPR"
  regression <- make_list(do_reg_Eti, do_reg_FPR,is_discrete_predictor)#, do_reg_TPR)
  
  # check BrS group:
  BrS_grp   <- FALSE
  prior_BrS <- model_options$prior$TPR_prior$BrS
  GBrS_TPR  <- length(unique(prior_BrS$grp))
  grp_spec  <- (!is.null(prior_BrS$grp) && GBrS_TPR >1 )
  if (grp_spec) {
    for (s in seq_along(prior_BrS$val)){
      if (prior_BrS$input=="match_range" && 
          (length(prior_BrS$val[[s]]$up)!=GBrS_TPR | 
           length(prior_BrS$val[[s]]$low)!=GBrS_TPR) ){
        stop(paste0("==[baker] ",names(prior_BrS$val)[s])," needs ", GBrS_TPR,
             " sets of sensitivity ranges.==")
      }
    }
    BrS_grp <- TRUE
  }
  
  
  # check SS group:
  SS_grp   <- FALSE
  prior_SS <- model_options$prior$TPR_prior$SS
  grp_spec <- (!is.null(prior_SS$grp) && length(unique(prior_SS$grp)) >1 )
  if (grp_spec) {SS_grp <- TRUE}
  
  ## <-------- the following are more strict grp specifications (may cause error when running old folders):
  #   val_spec <- (num_slice["MSS"]>0 && any(lapply(prior_SS$val,length)>1))
  #   
  #   if (grp_spec && val_spec){SS_grp <- TRUE}
  #   if (grp_spec && !val_spec){stop("==Specified TPR group in 'grp' of 'model_options$prior$TPR_prior$SS',
  #                                   but either there is no SS data or the length of 'val' does not match the no. of TPR groups. ==")}
  #   if (!grp_spec && val_spec){stop("==No 'grp' specified in 'model_options$prior$TPR_prior$SS',
  #                                   but we have >1 sets of TPRs. ==")}
  
  # return results:
  make_list(num_slice, nested, regression,BrS_grp,SS_grp)
}



