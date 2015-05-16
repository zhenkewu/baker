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
#' @param data_nplcm Data for model fitting. Details are
#' \itemize{
#' \item \code{Mobs} See \code{\link{nplcm}} function.
#' \item \code{Y} See \code{\link{nplcm}} function.
#' \item \code{X} See \code{\link{nplcm}} function.
#' \item \code{Mname} A list of pathogen names. \code{Mname_BrS} for BrS only pathogens,
#' and \code{Mname_SSonly} for pathogens with only SS measures.
#' \item \code{taxonomy} A list of pathogen taxonomy. The elements include 
#' \code{taxo_BrS} and \code{taxo_SSonly}. It can be used to represent biological
#' classifications, such as bacteria (\code{B}) versus virus (\code{V}).
#' }
#' 
#' @param model_options See \code{\link{nplcm}} function.
#' @param silent Default is \code{TRUE}: print assigned
#' model descriptions on screen; otherwise, \code{FALSE}.
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

assign_model <- function(data_nplcm,model_options,silent=TRUE){
  Mobs <- data_nplcm$Mobs
  Y    <- data_nplcm$Y
  X    <- data_nplcm$X
  
  use_data_sources <- model_options$M_use
  model_meas       <- model_options$k_subclass>1

  quality_nm <- names(Mobs)
  quality_nm[which(names(Mobs)=="MBS")] <- "BrS"
  quality_nm[which(names(Mobs)=="MSS")] <- "SS"
  quality_nm[which(names(Mobs)=="MGS")] <- "GS"
  
  input_data_sources  <- quality_nm[!is.na(Mobs) & 
                                      unlist(lapply(Mobs,function(a) !is.null(a)))]
  # the !is.null criteria is for situations where NULL is specified for MGS or MSS
  # in the list 'Mobs' (when you call 'Mobs$MGS', the output will be NULL).

  if (!all(input_data_sources%in%use_data_sources)){
    # test for data sources as specified by model_options and
    # the actual input data set:
    stop("==Please supply actual datasets as specified by 'M_use' in
         'model_options'.==")

  }else{
    if (length(model_options$TPR_prior)!=length(model_options$M_use)){
       stop("==Please input the same number of TPR priors 
            (i.e., length of 'TPR_prior') as in number of measuremens levels 
            (i.e., length of 'M_use') used in 'model_options'. ==")
    }else{
      assigned_model <- list(quality= paste(use_data_sources,collapse="+"),
                             SSonly = !is.null(model_options$pathogen_SSonly_list),
                             nest   = model_meas)
      
      
      # FPR regression:
      in_user <- !is.null(model_options$X_reg_FPR)
      in_data <- !is.null(X)
      if (!in_user && !in_data){
        do_FPR_reg <- FALSE
      } else if (in_user && !in_data){
        stop(paste0("==","Data does not have user specified covaraites for FPR regression!","=="))
      } else if (!in_user && in_data){
        do_FPR_reg <- FALSE
      } else {
        if (!all(model_options$X_reg_FPR%in%colnames(X))){
          stop(paste0("==Covariate(s), '",setdiff(model_options$X_reg_FPR,colnames(X)), 
                      "', not in data! =="))
        } else{
          do_FPR_reg <- TRUE
        }
      }
      
    
      
      # Etiology regression:
      in_user <- !is.null(model_options$X_reg_Eti)
      in_data <- !is.null(X)
      if (!in_user && !in_data){
        do_Eti_reg <- FALSE
      } else if (in_user && !in_data){
        stop(paste0("==","Data does not have user specified covaraites for Eti regression!","=="))
      } else if (!in_user && in_data){
        do_Eti_reg <- FALSE
      } else {
        if (!all(model_options$X_reg_Eti%in%colnames(X))){
          stop(paste0("==Covariate(s),",setdiff(model_options$X_reg_Eti,colnames(X)), 
                      " not in data! =="))
        } else{
          do_Eti_reg <- TRUE
        }
      }
      
      
      reg  <- list(do_FPR_reg = do_FPR_reg, do_Eti_reg = do_Eti_reg)
      list(measurement = assigned_model,reg = reg)
    }
  }
}






