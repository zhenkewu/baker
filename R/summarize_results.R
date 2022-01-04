#' `print.nplcm` summarizes the results from [nplcm()].
#'
#' @param x Output from [nplcm()].
#' @param ... Arguments passed to summary and printing methods.
#' @return Summary of object output by [nplcm()] --- need details.
#'
#' @family nplcm results
#' @export
print.nplcm <- function(x, ...){
  print(summary(x, ...), ...)
  # Return
  return(invisible(x))
}


#' `summary.nplcm` summarizes the results from [nplcm()].
#'
#' @param object Output from [nplcm()].
#' An object of class "nplcm"
#' @param res_dir file path to folder where results are stored.
#' @param type `NoReg`,`strat`, or `Reg`
#' @param ... Not used.
#' @return see [print.nplcm()]
#'
#' @family nplcm results
#' @export
summary.nplcm <- function(object,...){
  
  if (object$fitted_type=="no_reg"){
    res <- plot_panels(object$DIR_NPLCM, is_plot=FALSE)
  }
  
  if (object$fitted_type=="reg_nonest_strat" | object$fitted_type=="reg_nest_strat"){
    res <- plot_etiology_strat(object$DIR_NPLCM, strata_weights = "empirical",is_plot=FALSE)
  }
  
  if (object$fitted_type=="reg_nonest" | object$fitted_type=="reg_nest"){
    res <- plot_etiology_regression(object$DIR_NPLCM,do_plot = FALSE,plot_basis=FALSE,...)
  }
  res$fitted_type <- object$fitted_type
  class(res) <- paste0("summary.nplcm.",object$fitted_type)
  
  return(res)
  
}


#' Compact printing of [nplcm()] model fits
#'
#' `print.summary.nplcm` is a print method for class
#' `summary.nplcm.NoReg`.
#'
#' @param x output from `summary.nplcm` with `summary.nplcm.no_reg` as the output object class.
#' @param ... Not used.
#' @return see [print.nplcm()]
#'
#' @export
#' @family nplcm results
print.summary.nplcm.no_reg <- function(x,...) {
  cat("[baker] summary: model structure","\n")
  cat("           fitted type: ",x$fitted_type,"\n")
  cat("---\n")
  cat("     name measurements: ", names(x$parsed_model$num_slice),"\n")
  cat("slices of measurements: ", x$parsed_model$num_slice,"\n")
  cat("                nested: ", x$parsed_model$nested,"\n")
  cat("---\n")
  cat("            regression: ","\n")
  cat("                  etiology: ", x$parsed_model$regression$do_reg_Eti,"\n")
  cat("                  name FPR: ", names(x$parsed_model$regression$do_reg_FPR),"\n")
  cat("                       FPR: ", x$parsed_model$regression$do_reg_FPR,"\n")
  cat("---\n")
  cat("all discrete predictor: ","\n")
  cat("                  etiology: ", x$parsed_model$regression$is_discrete_predictor$Eti,"\n")
  cat("                  name FPR: ", names(x$parsed_model$regression$do_reg_FPR),"\n")
  cat("                       FPR: ", x$parsed_model$regression$is_discrete_predictor$FPR,"\n")
  
  cat("\n------- posterior summary -----------\n")
  
  print(x$res_Eti)
  # Return
  return(invisible(x))
}



#' Compact printing of [nplcm()] model fits
#'
#' `print.summary.nplcm` is a print method for class
#' `summary.nplcm.reg_nonest_strat`.
#'
#' @param x output from `summary.nplcm` with `summary.nplcm.reg_nonest_strat` as the output object class.
#' @param ... Not used.
#' @return see [print.nplcm()]
#'
#' @export
#' @family nplcm results
print.summary.nplcm.reg_nonest_strat <- function(x,...) {
  
  cat("[baker] summary: model structure","\n")
  cat("           fitted type: ",x$fitted_type,"\n")
  cat("---\n")
  cat("     name measurements: ", names(x$parsed_model$num_slice),"\n")
  cat("slices of measurements: ", x$parsed_model$num_slice,"\n")
  cat("                nested: ", x$parsed_model$nested,"\n")
  cat("---\n")
  cat("            regression: ","\n")
  cat("                  etiology: ", x$parsed_model$regression$do_reg_Eti,"\n")
  cat("                  name FPR: ", names(x$parsed_model$regression$do_reg_FPR),"\n")
  cat("                       FPR: ", x$parsed_model$regression$do_reg_FPR,"\n")
  cat("---\n")
  cat("all discrete predictor: ","\n")
  cat("                  etiology: ", x$parsed_model$regression$is_discrete_predictor$Eti,"\n")
  cat("                  name FPR: ", names(x$parsed_model$regression$do_reg_FPR),"\n")
  cat("                       FPR: ", x$parsed_model$regression$is_discrete_predictor$FPR,"\n")
  
  cat("\n------- strata definitions (by row) -----------\n")
  print(x$unique_Eti_level)
  
  cat("\n------- posterior summary -----------\n")
  
  print(x$res_list)
  
  # Return
  return(invisible(x))
}


#' Compact printing of [nplcm()] model fits
#' 
#' Same as [print.summary.nplcm.reg_nonest_strat()]
#'
#' `print.summary.nplcm` is a print method for class
#' `summary.nplcm.reg_nest_strat`.
#'
#' @param x output from `summary.nplcm` with `summary.nplcm.reg_nest_strat` as the output object class.
#' @param ... Not used.
#' @return see [print.nplcm()]
#'
#' @export
#' @family nplcm results
print.summary.nplcm.reg_nest_strat <- function(x,...) {
  
  cat("[baker] summary: model structure","\n")
  cat("           fitted type: ",x$fitted_type,"\n")
  cat("---\n")
  cat("     name measurements: ", names(x$parsed_model$num_slice),"\n")
  cat("slices of measurements: ", x$parsed_model$num_slice,"\n")
  cat("                nested: ", x$parsed_model$nested,"\n")
  cat("---\n")
  cat("            regression: ","\n")
  cat("                  etiology: ", x$parsed_model$regression$do_reg_Eti,"\n")
  cat("                  name FPR: ", names(x$parsed_model$regression$do_reg_FPR),"\n")
  cat("                       FPR: ", x$parsed_model$regression$do_reg_FPR,"\n")
  cat("---\n")
  cat("all discrete predictor: ","\n")
  cat("                  etiology: ", x$parsed_model$regression$is_discrete_predictor$Eti,"\n")
  cat("                  name FPR: ", names(x$parsed_model$regression$do_reg_FPR),"\n")
  cat("                       FPR: ", x$parsed_model$regression$is_discrete_predictor$FPR,"\n")
  
  cat("\n------- strata definitions (by row) -----------\n")
  print(x$unique_Eti_level)
  
  cat("\n------- posterior summary -----------\n")
  
  print(x$res_list)
  
  # Return
  return(invisible(x))
}

#' Compact printing of [nplcm()] model fits
#'
#' `print.summary.nplcm` is a print method for class
#' `summary.nplcm.reg_nonest`.
#'
#' @param x output from `summary.nplcm` with `summary.nplcm.reg_nonest` as the output object class.
#' @param ... Not used.
#' @return see [print.nplcm()]
#'
#' @export
#' @family nplcm results
print.summary.nplcm.reg_nonest <- function(x,...) {
  
  cat("[baker] summary: model structure","\n")
  cat("           fitted type: ",x$fitted_type,"\n")
  cat("---\n")
  cat("     name measurements: ", names(x$parsed_model$num_slice),"\n")
  cat("slices of measurements: ", x$parsed_model$num_slice,"\n")
  cat("                nested: ", x$parsed_model$nested,"\n")
  cat("---\n")
  cat("            regression: ","\n")
  cat("                  etiology: ", x$parsed_model$regression$do_reg_Eti,"\n")
  cat("                  name FPR: ", names(x$parsed_model$regression$do_reg_FPR),"\n")
  cat("                       FPR: ", x$parsed_model$regression$do_reg_FPR,"\n")
  cat("---\n")
  cat("all discrete predictor: ","\n")
  cat("                  etiology: ", x$parsed_model$regression$is_discrete_predictor$Eti,"\n")
  cat("                  name FPR: ", names(x$parsed_model$regression$do_reg_FPR),"\n")
  cat("                       FPR: ", x$parsed_model$regression$is_discrete_predictor$FPR,"\n")
  
  cat("\n------- posterior summary -----------\n")
  print(x$res)
  # Return
  return(invisible(x))
}


#' Compact printing of [nplcm()] model fits
#'
#' `print.summary.nplcm` is a print method for class
#' `summary.nplcm.reg_nest`.
#'
#' @param x output from `summary.nplcm` with `summary.nplcm.reg_nest` as the output object class.
#' @param ... Not used.
#' @return see [print.nplcm()]
#'
#' @export
#' @family nplcm results
print.summary.nplcm.reg_nest <- function(x,...) {
  
  cat("[baker] summary: model structure","\n")
  cat("           fitted type: ",x$fitted_type,"\n")
  cat("---\n")
  cat("     name measurements: ", names(x$parsed_model$num_slice),"\n")
  cat("slices of measurements: ", x$parsed_model$num_slice,"\n")
  cat("                nested: ", x$parsed_model$nested,"\n")
  cat("---\n")
  cat("            regression: ","\n")
  cat("                  etiology: ", x$parsed_model$regression$do_reg_Eti,"\n")
  cat("                  name FPR: ", names(x$parsed_model$regression$do_reg_FPR),"\n")
  cat("                       FPR: ", x$parsed_model$regression$do_reg_FPR,"\n")
  cat("---\n")
  cat("all discrete predictor: ","\n")
  cat("                  etiology: ", x$parsed_model$regression$is_discrete_predictor$Eti,"\n")
  cat("                  name FPR: ", names(x$parsed_model$regression$do_reg_FPR),"\n")
  cat("                       FPR: ", x$parsed_model$regression$is_discrete_predictor$FPR,"\n")
  
  cat("\n------- posterior summary -----------\n")
  print(x$res)
  # Return
  return(invisible(x))
}




