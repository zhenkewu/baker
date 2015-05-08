#' Choose a model to fit
#'
#' \code{assign_model} recognizes the model to fit as requested by user.
#' \code{assign_model} will also inspect the actual data supplied to this function
#'  and check if the data conform to user's requested model. The following features of data and user
#' inputs are checked: 1) available measurement quality types, i.e., gold-, silver- or bronze-standard or
#' any combinations; 2) model for false positive rates: covariate-dependent or not;
#' 3) model for etiology: covariate-dependent or not.
#' 
#'
#' @param Mobs See \code{\link{nplcm}} function.
#' @param Y See \code{\link{nplcm}} function.
#' @param X See \code{\link{nplcm}} function.
#' @param model_options See \code{\link{nplcm}} function.
#' @param silent Default is \code{TRUE}: print assigned
#' model descriptions on screen; otherwise, \code{FALSE}.
#'
#' @return A list of model information:
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

assign_model <- function(Mobs,Y,X,model_options,silent=TRUE){

  use_data_sources <- model_options$M_use
  model_meas       <- model_options$k_subclass>1

  quality_nm <- names(Mobs)
  quality_nm[which(names(Mobs)=="MBS")] <- "BrS"
  quality_nm[which(names(Mobs)=="MSS")] <- "SS"
  quality_nm[which(names(Mobs)=="MGS")] <- "GS"
  input_data_sources  <- quality_nm[!is.na(Mobs) & 
                                      unlist(lapply(Mobs,function(a) !is.null(a)))]
  # the !is.null criteria is for situations where NULL is specified for MGS or MSS
  # in the list 'Mobs'.

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
