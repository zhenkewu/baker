#' Interpret the model specified by user
#'
#' \code{assign_model} recognizes the model to fit from user input.
#' 
#' @details \code{assign_model} will be modified to check if data are conformable
#' to specified model.
#' @param data_nplcm Data for model fitting.
#' @param model_options See \code{\link{nplcm}} function.
#' @param silent Default is TRUE for no messages; FALSE otherwise.
#' @return A list of information for the selected model:
#' \itemize{
#'    \item \code{num_slice} a vector counting the no. of measurement slices for every
#'    level of measurement quality
#'    \item \code{nested} TRUE for nested models in BrS data (conditional dependence); 
#'    FALSE for non-nested models (conditional independence). One for each slice.
#'    \item \code{regression}
#'        \itemize{
#'            \item \code{do_reg_Eti} TRUE for doing regression on etiology (latent status); 
#'            FALSE otherwise
#'            \item \code{do_reg_FPR} TRUE for doing regression on false positive rates 
#'            (for every slice of bronze-standard); FALSE otherwise
#'        }
#' }
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
    stop("==Please supply actual datasets as specified by 'use_measurements' in 'model_options'.==")
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

# 
# #' check if SS measurements have stratified TPRs
# #' 
# #' @param model_options See \link{nplcm}
# #' 
# #' 
# #' @return TRUE for stratified TPRs; FALSE otherwise
# #' 
# #' @export
# check_SS_grp <- function(model_options){
#   res <- FALSE
#   prior_SS <- model_options$prior$TPR_prior$SS
#   if (!is.null(prior_SS$grp) && length(unique(prior_SS$grp)) >1 ){
#     res <- TRUE
#   }
#   res
# }


# #' check the compatibility of model and parameter specification
# #' 
# #' 
# #' 
# check_spec <- function(data_nplcm, model_options, clean_options){
#   
# }
# 



