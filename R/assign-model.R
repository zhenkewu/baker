#' Choose a model to fit
#'
#' \code{assign_model} recognizes the model to fit from user input.
#' 
#' @details \code{assign_model} will also inspect the actual data supplied 
#'  and check if the data conform to user's requested model. The following 
#'  features of data and user inputs are checked against each other: 
#' \enumerate{
#' \item Available types of measurement quality, 
#' i.e., gold-, silver- or bronze-standard or any combinations; 
#' \item Model for false positive rates: covariate-dependent or not;
#' \item Model for etiology: covariate-dependent or not.
#' }
#' @param data_nplcm Data for model fitting.
#' @param model_options See \code{\link{nplcm}} function.
#' 
#' @return A list of information for the selected model:
#' \itemize{
#' \item \code{measurement}
#' \itemize{
#' \item \code{quality} e.g. "BrS+SS" indicates both BrS and SS measures are 
#' available
#' \item \code{SSonly} \code{TRUE} for existence of pathogens with only SS measures;
#' otherwise, \code{FALSE};
#' \item \code{nest} \code{TRUE} for conditional dependent model; \code{FALSE}
#' for conditional independent model;
#' }
#' \item \code{reg}
#' \itemize{
#' \item \code{do_FPR_reg} \code{TRUE} for allowing FPR to be covariate-dependent; \code{FALSE} otherwise;
#' \item \code{do_Eti_reg} \code{TRUE} for allowing etiology to be covariate-dependent; \code{FALSE} otherwise;
#' }
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
  for (i in seq_along(use_measurements)){
    num_slice[i] <- num_slice[i]+length(Mobs[[i]])
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






