#' Set the etiology prior
#'
#' Assign etiology prior to the pathogens entering the model.
#'
#'
#' @param model_options See \code{nplcm_fit} function argument list.
#'
#'
#' @return The vector of etiology prior parameters
#'
#'
#' @export

set_prior_eti <- function(model_options){
        Jcause <- length(model_options$cause_list)
        cat("==Etiology priors: ==" ,"\n",model_options$Eti_prior,"\n")
        if (model_options$Eti_prior=="overall_uniform"){
#                 if (is.null(model_options$SSonly)||model_options$SSonly==FALSE){
                  alpha    <-  rep(1,Jcause)
#                 }else {
#                   JSS_only <- length(model_options$pathogen_SSonly_list)
#                   alpha    <-  rep(1,Jcause+JSS_only)
#                 }
        }

        if (model_options$Eti_prior=="0_1"){
#           if (is.null(model_options$SSonly)||model_options$SSonly==FALSE){
            alpha    <-  rep(.1,Jcause)
#           }else {
#             JSS_only <- length(model_options$pathogen_SSonly_list)
#             alpha    <-  rep(.1,Jcause+JSS_only)
#           }
        }

        alpha
}

