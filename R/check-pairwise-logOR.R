#' Posterior predictive checking for nested partially latent class models - 
#' pairwise log odds ratio
#' 
#' At each MCMC iteration, we generate a new data set based on the model and 
#' parameter values at that iteration. The sample size of the new data set equals
#' that of the actual data set, i.e. the same number of cases and controls.
#' 
#' @param DIR_NPLCM File path to the folder that stores results from npLCM fit.
#' 
#' @return A figure of posterior predicted log odds ratio compared with the observed 
#' log odds ratio for the BrS data. The function generates this figure in your working directory automatically.
#' 
#' @export
#' 

check_pairwise_logOR <- function(DIR_NPLCM){
      #
      # read in the data:
      #  
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
      pathogen_BrS_list <- model_options$pathogen_BrS_list
  
      #
      # test: (pathogen_display can be separately specified)
      #
      pathogen_display <- pathogen_BrS_list
      index_display    <- my_reorder(pathogen_display,pathogen_BrS_list)
      pathogen_name    <- pathogen_BrS_list[index_display]
      
  
      JBrS  <- bugs.dat$JBrS
      Nd <- bugs.dat$Nd
      Nu <- bugs.dat$Nu
      
      MBS.case <- bugs.dat$MBS[Y==1,index_display]
      MBS.ctrl <- bugs.dat$MBS[Y==0,index_display]
      
      MBS.new <- res_nplcm[,grep("MBS.new",colnames(res_nplcm))]
      Niter   <- nrow(MBS.new)
      
      
      #
      # Figure 1: Standardized Log OR differences; Upper for cases, lower
      #           for controls:
      # 
      
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
                        function(v) sd(v[!is.na(v)]))
      misfit.vis.mat <- (logORmat-mean.logORmat.ppd)/sd.logORmat.ppd
      op <- par()
      pdf(paste0(DIR_NPLCM,"//predicted_vs_observed_logOR_standardized.pdf"),
          width=15,height=12)
      visualize_case_control_matrix(misfit.vis.mat,pathogen_name,cell_metrics="(obs-mean)/sd")
      dev.off()
      
      cat("==A figure is generated for model checking: pairwise log odds ratio. 
          Stored in ",DIR_NPLCM," ==")
      #
      # Figure 2: Boxplot of posterior predictive distribution and observed log OR:
      #
}

# @param pathogen_display The pathogen vector in desired order for display.
# It can be of larger length than that of \code{pathogen_BrS}.