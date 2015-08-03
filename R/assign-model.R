#' Interpret the model specified by user
#'
#' \code{assign_model} recognizes the model to fit from user input.
#' 
#' @details \code{assign_model} will be modified to check if data are conformable
#' to specified model.
#' @param data_nplcm Data for model fitting.
#' @param model_options See \code{\link{nplcm}} function.
#' 
#' @return A list of information for the selected model:
#' \itemize{
#'    \item \code{num_slice} a vector counting the no. of measurement slices for every
#'    level of measurement quality
#'    \item \code{nested} TRUE for nested models (conditional dependence); 
#'    FALSE for non-nested models (conditional independence).
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

assign_model <- function(model_options,data_nplcm){
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
  do_reg_FPR <- list()
  for (i in seq_along(Mobs$MBS)) {
    ind_tmp <-
      which(names(likelihood$FPR_formula) == names(Mobs$MBS)[i])
    if (!length(ind_tmp)) {
      do_reg_FPR[[i]] <- FALSE
    } else{
      do_reg_FPR[[i]] <-
        parse_nplcm_reg(as.formula(likelihood$FPR_formula[[ind_tmp]]),data_nplcm)
    }
  }
  names(do_reg_FPR) <- names(Mobs$MBS)
  
  # specify regression for etiology:
  form_tmp <- as.formula(likelihood$Eti_formula)
  do_reg_Eti <- parse_nplcm_reg(form_tmp,data_nplcm)
  regression <- make_list(do_reg_Eti, do_reg_FPR)
  make_list(num_slice, nested, regression)
}






