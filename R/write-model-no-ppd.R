#' Write .bug model file for no regression, conditional independence
#' @param Mobs measurements collected by `data_nplcm`
#' @param cause_list the vector of latent status
#' @param use_measurements any single or combinations of "BrS", "SS" 
#' 
#' @return a long character string to be written as .bug file.
#' @export
write_model_NoReg_plcm <- function(Mobs,cause_list,use_measurements){
  ## 1) accommodates singletons, combos, and NoA;
  ## 2) check-bit to prevent case data informing FPR;
  
  chunk1 <- insert_bugfile_chunk_noreg_meas(Mobs,cause_list,use_measurements)
  chunk2 <- insert_bugfile_chunk_noreg_etiology()

  paste0("model{#BEGIN OF MODEL:\n",
         chunk1,"\n",
         chunk2,"\n",
         "}#END OF MODEL.")
  
}