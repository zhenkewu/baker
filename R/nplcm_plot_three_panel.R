#' Plot three-panel figures for nested partially-latent model results
#'
#' Visualize the model outputs for communicating how data inform the
#' DN: 1. current implementation: nplcm, BrS and SS.
#' "Jfull" here is not the same as in other functions: it refers to the number of
#' pathogens even if there are pathogens with only silver-standard data;
#' in other functions, "Jfull" refers to the number of pathogens that have BrS data.
#' 2. Missing data for BrS or SS are dropped when calculating observed measurement
#' prevalences
#'
#' @param DIR_NPLCM File path to the folder containing posterior samples
#'
#' @param ss_upperlimit The upper limit of horizontal bar for the silver-standard
#' subpanel (the middle panel). The default value is .25.
#'
#' @param eti_upperlimit The upper limit of horizontal bar for the etiology
#' posterior subpanel (the rightmost panel). The default value is .4
#' @importFrom coda read.coda
#' @importFrom binom binom.confint
#' @return None
#'
#' @export
nplcm_plot_three_panel <- function(DIR_NPLCM,
                                   ss_upperlimit=1,
                                   eti_upperlimit=.4){#BEGIN function
  
  # read NPLCM outputs:
  out           <- nplcm_read_folder(DIR_NPLCM)
  # organize ouputs:
  Mobs          <- out$Mobs
  Y             <- out$Y
  model_options <- out$model_options
  clean_options <- out$clean_options
  res_nplcm     <- out$res_nplcm
  bugs.dat      <- out$bugs.dat
  rm(out)
  
  # Determine which three-panel plot to draw:
  parsing <- assign_model(Mobs,Y,X,model_options)
  # X not needed in the three-panel plot, but because 'assign_model' was designed
  # to distinguish models even with X, so we have to stick to the useage of 
  # assign_model.
    
  if (!any(unlist(parsing$reg))){
    # if no stratification or regression:
    if (parsing$measurement$quality=="BrS+SS"){
      if (!parsing$measurement$SSonly){
        if (!parsing$measurement$nest){
          print("== BrS+SS; no SSonly; subclass number: K = 1. ==")
        }else{
          print("== BrS+SS; no SSonly; subclass number: K > 1. ==")
        }
        #
        # Plot:
        # 
        
        # put layout in place for three panels:
        layout(matrix(c(1,2,3),1,3,byrow = TRUE),
               widths=c(3,2,3),heights=c(8))

        # BrS panel:
        nplcm_plot_BrS_panel(Mobs$MBS,model_options,clean_options,res_nplcm,
                             bugs.dat,
                             top_BrS = 1.3,prior_shape = "interval")

        # SS panel:
        nplcm_plot_SS_panel(Mobs$MSS,model_options,clean_options,res_nplcm,
                            bugs.dat,
                            top_SS = 0.3)
  
        # Etiology panel:
        nplcm_plot_pie_panel(model_options,res_nplcm,bugs.dat,
                             top_pie = 1)
          
      } else{
        if (!parsing$measurement$nest){
          stop("== BrS+SS; SSonly; subclass number: K = 1: not done.  ==")

        }else{
          stop("== BrS+SS; SSonly; subclass number: K > 1: not done.  ==")
        }
      }
    }else if (parsing$measurement$quality=="BrS"){
      if (!parsing$measurement$SSonly){
        if (!parsing$measurement$nest){
          stop("== BrS; no SSonly; subclass number: K = 1: not done.  ==")
        }else{
          stop("== BrS; no SSonly; subclass number: K > 1: not done.  ==")
        }
      }
    }
  } else{
    stop("== Three panel plot not implemented for stratification or regression
         settings. Please check back later for updates. Thanks. ==")
  }

}# END function
