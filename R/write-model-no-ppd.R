#' Write .bug model file for no regression, conditional independence
#' 
#' 
#' @inheritParams insert_bugfile_chunk_noreg_meas
#' 
#' @return a long character string to be written as .bug file.
#' @seealso \link{insert_bugfile_chunk_noreg_meas} for inserting .bug file
#' chunk for measurements (plug-and-play); \link{insert_bugfile_chunk_noreg_etiology}
#' for inserting .bug file chunk for distribution of latent status (etiology).
#' 
#' @export
write_model_NoReg_NoNest <- function(Mobs,prior,cause_list,use_measurements){
  ## 1) accommodates singletons, combos, and NoA;
  ## 2) check-bit to prevent case data informing FPR;
  
  chunk1 <- insert_bugfile_chunk_noreg_meas(Mobs,prior,cause_list,use_measurements)
  chunk2 <- insert_bugfile_chunk_noreg_etiology()

  paste0("model{#BEGIN OF MODEL:\n",
         chunk1,"\n",
         chunk2,"\n",
         "}#END OF MODEL.")
  
}