#' Assign model based on model specification and data structure.
#'
#' Before fitting the model, this function determines a) the available measurement
#' types, b) false positive rate and c) etiology model components based on input
#' data and model options. It also does consistency checking for the specifications 
#' in \code{model_options}, and their relations with actual inputted data.
#'
#' @param Mobs See the \code{nplcm} function.
#' @param Y See the \code{nplcm} function.
#' @param X See the \code{nplcm} function.
#' @param model_options See the \code{nplcm} function.
#' @param silent Default is \code{TRUE}: print assigned
#' model descriptions on screen; otherwise, \code{FALSE}.
#'
#' @return The correct Bayesian model type to fit. The list of model
#' type information is:
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
#' \item \code{do_FPR_reg} \code{TRUE} for FPR regression; \code{FALSE} otherwise;
#' \item \code{do_Eti_reg} \code{TRUE} for etiology regression; \code{FALSE} otherwise;
#' }
#' }
#'
#' @export

assign_model <- function(Mobs,Y,X,model_options,silent=TRUE){

  use_data_sources <- model_options$M_use
  model_meas       <- model_options$k_subclass>1

  quality_nm <- names(Mobs)
  quality_nm[which(names(Mobs)=="MBS")] <- "BrS"
  quality_nm[which(names(Mobs)=="MSS")] <- "SS"
  quality_nm[which(names(Mobs)=="MGS")] <- "GS"
  input_data_sources  <- quality_nm[!is.na(Mobs)]

  if (!all(input_data_sources%in%use_data_sources)){
    # test for data sources as specified by model_options and
    # the actual input data set:
    stop("==Please supply actual datasets as specified by 'M_use' in
         'model_options'.==")

  }else{
    if (length(model_options$TPR_prior)!=length(model_options$M_use)){
       stop("==Please input same number of TPR priors (i.e., length of 'TPR_prior') 
         as in number of measuremens levels (i.e., length of 'M_use') used 
            in 'model_options'. ==")
    }else{
      assigned_model <- list(quality= paste(use_data_sources,collapse="+"),
                             SSonly = !is.null(model_options$pathogen_SSonly_list),
                             nest   = model_meas)
                reg  <- list(do_FPR_reg = !is.null(model_options$X_reg_FPR),
                             do_Eti_reg = !is.null(model_options$X_reg_Eti))
    }
  }

  list(      measurement = assigned_model,
                     reg = reg)
}

#assign_model(Mobs,Y,X,model_options)
