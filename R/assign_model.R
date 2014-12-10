#' Assign model based on model specification and data.
#'
#' For a data set fed into \code{\link{nplcm}}, determine measurement and
#' false positive rate and etiology model components. It also does consistency
#' checking for specifications in \code{model_options}, and their relations
#' with actual inputted data.
#'
#' @param Mobs See the \code{nplcm} function.
#' @param Y See the \code{nplcm} function.
#' @param X See the \code{nplcm} function.
#' @param model_options See the \code{nplcm} function.
#' @param silent Default is \code{TRUE}, which means print on screen assigned
#' model descriptions, otherwise, \code{FALSE}.
#'
#' @return Assign the model type to be fit.
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
    stop("==Please input actual data with data sources specified by 'M_use' in
         'model_options'.==")

  }else{
    if (length(model_options$TPR_prior)!=length(model_options$M_use)){
       stop("==Please input same number of TPR priors (length of 'TPR_prior') as in number of measurement
          levels (length of 'M_use') used in 'model_options'. ==")
    }else{
      assigned_model <- list( quality= paste(use_data_sources,collapse="+"),
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
