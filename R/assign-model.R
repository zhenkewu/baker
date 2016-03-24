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
#'        }
#' }
#'
#' @family specification checking functions
#' 
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
  
  # specify regression for FPR:
  do_reg_FPR <- list() #  <---- a regression for each measurement slice?
  for (i in seq_along(Mobs$MBS)) {
    ind_tmp <-
      which(names(likelihood$FPR_formula) == names(Mobs$MBS)[i])
    if (!length(ind_tmp)) { # don't do regression if no regression formula is found:
      do_reg_FPR[[i]] <- FALSE
    } else{ # do regression if there is matched regression formula:
      do_reg_FPR[[i]] <-
        parse_nplcm_reg(as.formula(likelihood$FPR_formula[[ind_tmp]]),data_nplcm,silent=silent)
    }
  }
  names(do_reg_FPR) <- names(Mobs$MBS)
  
  # specify regression for etiology:
  form_tmp   <- as.formula(likelihood$Eti_formula)
  do_reg_Eti <- parse_nplcm_reg(form_tmp,data_nplcm,silent=silent)
  regression <- make_list(do_reg_Eti, do_reg_FPR)
  
  # check SS group:
  SS_grp <- FALSE
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
  make_list(num_slice, nested, regression,SS_grp)
}



