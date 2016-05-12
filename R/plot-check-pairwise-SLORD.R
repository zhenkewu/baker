#' Posterior predictive checking for nested partially latent class models - 
#' pairwise log odds ratio (only for bronze-standard data)
#' 
#' At each MCMC iteration, we generate a new data set based on the model and 
#' parameter values at that iteration. The sample size of the new data set equals
#' that of the actual data set, i.e. the same number of cases and controls.
#' 
#' @param DIR_NPLCM File path to the folder that stores results from npLCM fit.
#' @param slice Default is 1, for the first slice of BrS data.
#' @return A figure of posterior predicted log odds ratio compared with the observed 
#' log odds ratio for the BrS data. The function generates this figure in your working directory automatically.
#' @family visualization functions
#' @family model checking functions
#' @export
#' 

plot_check_pairwise_SLORD <- function(DIR_NPLCM,slice = 1){
  old_par <- graphics::par(graphics::par("mfrow", "mar"))
  on.exit(graphics::par(old_par))
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
  
  curr_BrS_object <-clean_options$BrS_objects[[slice]]
  pathogen_BrS_list <- curr_BrS_object$patho
  
  #
  # test: (pathogen_display can be separately specified)
  #
  pathogen_display <- pathogen_BrS_list
  index_display    <- my_reorder(pathogen_display,pathogen_BrS_list)
  pathogen_name    <- pathogen_BrS_list[index_display]
  JBrS     <- length(pathogen_name)
  
  Nd <- bugs.dat$Nd
  Nu <- bugs.dat$Nu
  
  curr_MBS <- bugs.dat[[grep(paste0("^MBS_",slice,"$"),names(bugs.dat))]]
  
  MBS.case <- curr_MBS[Y==1,index_display]
  MBS.ctrl <- curr_MBS[Y==0,index_display]
  
  is_jags <- is_jags_folder(DIR_NPLCM)
  
  MBS.new <- res_nplcm[,grep(paste0("^MBS.new_",slice,"\\["),colnames(res_nplcm))]
  
  # to handle JAGS outputs' transposition:
  if (is_jags){
    swap_col <- c(sapply(1:(Nd+Nu),function(s) s+((1:JBrS)-1)*(Nd+Nu)))
    MBS.new  <- MBS.new[,swap_col]
  }
  
  Niter   <- nrow(MBS.new)
  
  
  # Figure 1: Standardized Log OR differences; Upper for cases, lower
  #           for controls:
  
  # observed:
  out <- logOR(MBS.case,MBS.ctrl)
  logORmat       <- out$logOR
  logORmat.se    <- out$logOR.se
  
  # posterior predicted:
  logORmat.ppd     <- array(NA,c(JBrS,JBrS,Niter))
  logORmat.se.ppd  <- array(NA,c(JBrS,JBrS,Niter))
  std.logORmat.ppd <- array(NA,c(JBrS,JBrS,Niter))
  
  for (iter in 1:Niter){
    mat_case <- matrix(MBS.new[iter,1:(JBrS*Nd)],ncol=JBrS,byrow=TRUE)[,index_display]
    mat_ctrl <- matrix(MBS.new[iter,-(1:(JBrS*Nd))],ncol=JBrS,byrow=TRUE)[,index_display]
    
    out.tmp    <- logOR(mat_case,mat_ctrl)
    logORmat.ppd[,,iter]    <- out.tmp$logOR
    logORmat.se.ppd[,,iter] <- out.tmp$logOR.se
  }
  
  misfit.vis.mat <- matrix(NA,nrow=JBrS,ncol=JBrS)
  mean.logORmat.ppd <- apply(logORmat.ppd,c(1,2),
                             function(v) mean(v[!is.na(v)]))
  sd.logORmat.ppd <- apply(logORmat.ppd,c(1,2),
                           function(v) stats::sd(v[!is.na(v)]))
  misfit.vis.mat <- (logORmat-mean.logORmat.ppd)/sd.logORmat.ppd
  
  visualize_case_control_matrix(misfit.vis.mat,pathogen_name,cell_metrics="(obs-mean)/stats::sd")
  
  cat("==A figure is generated for model checking: pairwise standardized log odds ratio difference (SLORD). 
          Stored in ",DIR_NPLCM," ==")
}
