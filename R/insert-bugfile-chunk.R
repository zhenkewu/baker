#' Insert measurement likelihood (without regression) code chunk into .bug model file
#' 
#' @param k_subclass the number of subclasses for the slices that require conditional dependence modeling; its length is
#' of the same value as the number of BrS slices.
#' @param Mobs measurement data in the form of \code{data_nplcm}
#' @param prior prior information from \code{model_options}
#' @param cause_list a list of latent status (crucial for building templates)
#' @param use_measurements "BrS", or "SS"
#' 
#' @return a long character string to be inserted into target .bug model file
#' 
#' @seealso \link{write_model_NoReg} for constructing a .bug file
#' 
#' @export
insert_bugfile_chunk_noreg_meas <-
  function(k_subclass,Mobs,prior,cause_list,use_measurements = "BrS") {
    if (!("BrS" %in% use_measurements) && !("SS" %in% use_measurements)){stop("==No BrS or SS measurements specified in the model! ==")}
    for (s in seq_along(Mobs$MBS)){
      if (k_subclass[s]>1 && ncol(Mobs$MBS[[s]])==1){stop("==cannot do nested modeling for BrS measurements with only one column! ==")}  
    }
    
    # generate file:
    if ("BrS" %in% use_measurements){
      chunk_BrS_case <- ""
      chunk_BrS_ctrl <- ""
      chunk_BrS_subclass <- ""
      chunk_BrS_param <- ""
      for (s in seq_along(Mobs$MBS)){# begin iterate over slices:
        k_curr <- k_subclass[s]
        if (k_curr == 1 ){
          chunk_BrS_case  <- paste0(chunk_BrS_case,  add_meas_BrS_case_NoNest_Slice(s,Mobs,cause_list)$plug)
          chunk_BrS_ctrl  <- paste0(chunk_BrS_ctrl,  add_meas_BrS_ctrl_NoNest_Slice(s,Mobs,cause_list)$plug)
          chunk_BrS_param <- paste0(chunk_BrS_param, add_meas_BrS_param_NoNest_Slice(s,Mobs,cause_list)$plug)
        }
        
        if ((k_curr > 1 )){
          chunk_BrS_case  <- paste0(chunk_BrS_case,  add_meas_BrS_case_Nest_Slice(s,Mobs,cause_list)$plug)
          chunk_BrS_ctrl  <- paste0(chunk_BrS_ctrl,  add_meas_BrS_ctrl_Nest_Slice(s,Mobs,cause_list)$plug)
          chunk_BrS_subclass  <- paste0(chunk_BrS_subclass,  add_meas_BrS_subclass_Nest_Slice(s,Mobs,cause_list)$plug)
          chunk_BrS_param <- paste0(chunk_BrS_param, add_meas_BrS_param_Nest_Slice(s,Mobs,cause_list)$plug)
        }
      }# end iterate over slices.
    }

    if ("BrS" %in% use_measurements & !("SS" %in% use_measurements)) {
      chunk <- paste0(
        "# BrS measurements:
        for (i in 1:Nd){
        ",chunk_BrS_case,"
        }
        for (i in (Nd+1):(Nd+Nu)){
        ",chunk_BrS_ctrl,"
        }
 
        ",chunk_BrS_subclass,"

        # bronze-standard measurement characteristics:
        ",chunk_BrS_param
      )
    }
    
    if (!("BrS" %in% use_measurements) & ("SS" %in% use_measurements)) {
      nslice_SS <- length(Mobs$MSS)
      chunk <- paste0(
        "# SS measurements:
        for (i in 1:Nd){
        ",add_meas_SS_case(nslice_SS,Mobs,prior,cause_list)$plug,"
        }
        
        # silver-standard measurement characteristics:
        ",add_meas_SS_param(nslice_SS,Mobs,prior,cause_list)$plug
      )
    }
    
    if ("BrS" %in% use_measurements & ("SS" %in% use_measurements)) {
      nslice_SS <- length(Mobs$MSS)
      chunk <- paste0(
        "# BrS and SS measurements:
        for (i in 1:Nd){
        ",chunk_BrS_case,"
        ",add_meas_SS_case(nslice_SS,Mobs,prior,cause_list)$plug,"
        }
        for (i in (Nd+1):(Nd+Nu)){
        ",chunk_BrS_ctrl,"
        }
        
        ",chunk_BrS_subclass,"

        # measurement characteristics:
        ",chunk_BrS_param,
          add_meas_SS_param(nslice_SS,Mobs,prior,cause_list)$plug
      )
    }
    paste0(chunk,"\n")
}



#' insert etiology code chunks into .bug file
#' 
#' @return a long character string to be inserted into target .bug model file
#' @export

insert_bugfile_chunk_noreg_etiology <- function(){

  chunk_etiology <- paste0("
  # etiology priors
  for (i in 1:Nd){
    Icat[i] ~ dcat(pEti[1:Jcause])
    
  }
  pEti[1:Jcause]~ddirch(alpha[])")
  
  paste0(chunk_etiology,"\n")
}
