#' Write .bug model file for model without regression
#' 
#' \code{write_model_NoReg} automatically generates model file according to
#' \code{model_options}
#' 
#' @inheritParams insert_bugfile_chunk_noreg_meas
#' 
#' @return a long character string to be written into .bug file.
#' 
#' @seealso 
#' \link{insert_bugfile_chunk_noreg_meas} for inserting .bug file
#' chunk for measurements (plug-and-play); \link{insert_bugfile_chunk_noreg_etiology}
#' for inserting .bug file chunk for distribution of latent status (etiology).
#' 
#' @family model generation functions
#' 
#' @export
#' 
write_model_NoReg <- function(k_subclass,Mobs,prior,cause_list,
                              use_measurements,ppd=NULL,use_jags=FALSE){
  ## 1) accommodate singletons, combos, and NoA;
  ## 2) If use_jags=FALSE, check-bit to prevent case data informing FPR (see the use of cut() functions in "plug-and-play.R");
  
  chunk1 <- insert_bugfile_chunk_noreg_meas(k_subclass,Mobs,
                                            prior,cause_list,use_measurements,ppd,use_jags)
  chunk2 <- insert_bugfile_chunk_noreg_etiology(ppd)
  
  paste0("model{#BEGIN OF MODEL:\n",
         chunk1,"\n",
         chunk2,"\n",
         "}#END OF MODEL.")
}

#' Write .bug model file for regression model without nested subclasses
#' 
#' \code{write_model_Reg_NoNest} automatically generates model file according to
#' \code{model_options}
#' 
#' @inheritParams insert_bugfile_chunk_reg_nonest_meas
#' 
#' @return a long character string to be written into .bug file.
#' 
#' @seealso 
#' \link{insert_bugfile_chunk_noreg_meas} for inserting .bug file
#' chunk for measurements (plug-and-play); \link{insert_bugfile_chunk_noreg_etiology}
#' for inserting .bug file chunk for distribution of latent status (etiology).
#' 
#' @family model generation functions
#' 
#' @export
#' 
write_model_Reg_NoNest <- function(Mobs,prior,cause_list,
                                   use_measurements,ppd=NULL,use_jags=FALSE){
  chunk1 <- insert_bugfile_chunk_reg_nonest_meas(Mobs,
                                                 prior,cause_list,use_measurements,ppd,use_jags)
  chunk2 <- insert_bugfile_chunk_reg_etiology(ppd) #DONE.
  
  paste0("model{#BEGIN OF MODEL:\n",
         chunk1,"\n",
         chunk2,"\n",
         "}#END OF MODEL.")
}
