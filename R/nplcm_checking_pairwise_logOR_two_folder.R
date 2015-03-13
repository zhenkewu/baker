#' Posterior predictive checking for nested partially latent class models - 
#' pairwise log odds ratio (parallel comparisons for two folders)
#' 
#' At each MCMC iteration, we generate a new data set based on the model and 
#' parameter values at that iteration. The sample size of the new data set equals
#' that of the actual data set, i.e. the same number of cases and controls.
#' 
#' @param DIR_NPLCM File path to the folder that stores results from npLCM fit.
#' 
#' @importFrom coda read.coda
#' 
#' @return A figure of posterior predicted log odds ratio compared with the observed 
#' log odds ratio for the BrS data. The function generates this figure in your working directory automatically.
#' 
#' @export
#' 

nplcm_checking_pairwise_logOR_two_folder <- function(DIR_NPLCM1,DIR_NPLCM2){
  
  obs_logOR_list   <- list()
  obs_logORse_list   <- list()
  pred_list <- list()
  Niter_list <- list()
  misfit.vis.mat_list <- list()
  count <- 1
  for (DIR_NPLCM in c(DIR_NPLCM1,DIR_NPLCM2)){
      #
      # read in the data:
      #  
      # remember that the data.txt file in the WinBUGS working folder is transposed:
      bugs.dat <- dget(paste(DIR_NPLCM,"data.txt",sep="/"))
      for (bugs.variable.name in names(bugs.dat)){
        if (!is.null(dim(bugs.dat[[bugs.variable.name]]))) {
          dim(bugs.dat[[bugs.variable.name]]) <- rev(dim(bugs.dat[[bugs.variable.name]]))
          bugs.dat[[bugs.variable.name]] <- aperm(bugs.dat[[bugs.variable.name]])
        }
        assign(bugs.variable.name, bugs.dat[[bugs.variable.name]])
      }
      
      model_options     <- dget(paste(DIR_NPLCM,"model_options.txt",sep="/"))
      pathogen_BrS_list <- model_options$pathogen_BrS_list
      
      #
      # test: (pathogen_display can be separately specified)
      #
      pathogen_display <- pathogen_BrS_list
      index_display    <- my_reorder(pathogen_display,pathogen_BrS_list)
      pathogen_name    <- pathogen_BrS_list[index_display]
      
      ## reading nplcm outputs:
      res_nplcm <- read.coda(paste(DIR_NPLCM,"coda1.txt",sep="/"),
                             paste(DIR_NPLCM,"codaIndex.txt",sep="/"),
                             quiet=TRUE)
      
      JBrS  <- bugs.dat$JBrS
      Nd <- bugs.dat$Nd
      Nu <- bugs.dat$Nu
      Y  <- c(rep(1,Nd),rep(0,Nu))
      
      MBS.case <- bugs.dat$MBS[Y==1,index_display]
      MBS.ctrl <- bugs.dat$MBS[Y==0,index_display]
      
      MBS.new <- res_nplcm[,grep("MBS.new",colnames(res_nplcm))]
      Niter   <- nrow(MBS.new)
      
    
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
      
      
      obs_logOR_list[[count]] <- logORmat
      obs_logORse_list[[count]] <- logORmat.se
      pred_list[[count]] <- logORmat.ppd
      Niter_list[[count]]  <- Niter
      misfit.vis.mat_list[[count]] <- misfit.vis.mat
        count <- count +1
    }# END loop over folders storing results.

  
  #
  # organize variables for plotting purposes:
  #
  
  # vector to store observed log odds ratios
  vec.logOR.obs.case <- rep(NA,JBrS*(JBrS-1)/2)
  vec.logOR.obs.ctrl <- rep(NA,JBrS*(JBrS-1)/2)
  
  # vector to store size of the dots (related to standard error of log odds ratios):
  vec.case.size      <- rep(NA,JBrS*(JBrS-1)/2)
  vec.ctrl.size      <- rep(NA,JBrS*(JBrS-1)/2)
  
  # posterior predictive log odds ratios; folder 1
  mat.logOR.ppd.case.1 = matrix(NA,nrow=Niter_list[[1]],ncol=JBrS*(JBrS-1)/2)
  mat.logOR.ppd.ctrl.1 = matrix(NA,nrow=Niter_list[[1]],ncol=JBrS*(JBrS-1)/2)
  
  # posterior predictive log odds ratios; folder 2
  mat.logOR.ppd.case.2 = matrix(NA,nrow=Niter_list[[2]],ncol=JBrS*(JBrS-1)/2)
  mat.logOR.ppd.ctrl.2 = matrix(NA,nrow=Niter_list[[2]],ncol=JBrS*(JBrS-1)/2)
  
  pair.count = 0
  plot.ind.case.1  = matrix(FALSE,nrow=2,ncol=choose(JBrS,2))
  plot.ind.ctrl.1  = matrix(FALSE,nrow=2,ncol=choose(JBrS,2))
  
  plot.ind.case.2  = matrix(FALSE,nrow=2,ncol=choose(JBrS,2))
  plot.ind.ctrl.2  = matrix(FALSE,nrow=2,ncol=choose(JBrS,2))
  for (j1 in 1:(JBrS-1)){
    for (j2 in (j1+1):JBrS){
      pair.count = pair.count+1
      
      vec.case.size[pair.count] = 1/(obs_logORse_list[[1]][j1,j2])^2
      vec.ctrl.size[pair.count] = 1/(obs_logORse_list[[1]][j2,j1])^2
      
      names(vec.logOR.obs.case)[pair.count] = paste0("(",j1,",",j2,")")
      names(vec.logOR.obs.ctrl)[pair.count] = paste0("(",j1,",",j2,")")
      
      vec.logOR.obs.case[pair.count] = obs_logOR_list[[1]][j1,j2]
      vec.logOR.obs.ctrl[pair.count] = obs_logOR_list[[1]][j2,j1]
      
      mat.logOR.ppd.case.1[,pair.count] = pred_list[[1]][j1,j2,]
      mat.logOR.ppd.ctrl.1[,pair.count] = pred_list[[1]][j2,j1,]
      
      mat.logOR.ppd.case.2[,pair.count] = pred_list[[2]][j1,j2,]
      mat.logOR.ppd.ctrl.2[,pair.count] = pred_list[[2]][j2,j1,]
      
      
      tt = misfit.vis.mat_list[[1]][j1,j2]
      if (abs(tt)>2 & is.finite(tt) & !is.na(tt) ){
        plot.ind.case.1[1,pair.count]=TRUE
        plot.ind.case.1[2,pair.count]=ifelse(tt>2,TRUE,FALSE)
      }
      tt = misfit.vis.mat_list[[1]][j2,j1]
      if (abs(tt)>2 & is.finite(tt) & !is.na(tt) ){
        plot.ind.ctrl.1[1,pair.count]=TRUE
        plot.ind.ctrl.1[2,pair.count]=ifelse(tt>2,TRUE,FALSE)
      }
      tt = misfit.vis.mat_list[[2]][j1,j2]
      if (abs(tt)>2 & is.finite(tt) & !is.na(tt) ){
        plot.ind.case.2[1,pair.count]=TRUE
        plot.ind.case.2[2,pair.count]=ifelse(tt>2,TRUE,FALSE)
      }
      tt = misfit.vis.mat_list[[2]][j2,j1]
      if (abs(tt)>2 & is.finite(tt) & !is.na(tt) ){
        plot.ind.ctrl.2[1,pair.count]=TRUE
        plot.ind.ctrl.2[2,pair.count]=ifelse(tt>2,TRUE,FALSE)
      }
    }
  }
  ind.both.case = (plot.ind.case.1[1,] | plot.ind.case.2[1,])
  ind.both.ctrl = (plot.ind.ctrl.1[1,] | plot.ind.ctrl.2[1,])
  
  ind.axis.show = (ind.both.case | ind.both.ctrl)
  
  mat.logOR.ppd.case.1[,!ind.both.case]=NA
  mat.logOR.ppd.ctrl.1[,!ind.both.ctrl]=NA
  mat.logOR.ppd.case.2[,!ind.both.case]=NA
  mat.logOR.ppd.ctrl.2[,!ind.both.ctrl]=NA
  
  vec.case.size = vec.case.size/max(vec.case.size[!is.na(vec.case.size)])*4
  vec.ctrl.size = vec.ctrl.size/max(vec.ctrl.size[!is.na(vec.ctrl.size)])*4
  
  
  #
  # start plotting:
  #
  at.loc = 1:choose(JBrS,2)
  boxwex = 0.4
  loc.gap = boxwex/2
  
  
  pdf(paste0("parallel_ppd_stat_excess_2.pdf"),
      width=10,height=10)
  op<-par()
  par(mfrow=c(2,1))
  
  #
  # cases:
  #
  par(mar = c(0,6,2,2))
  
  vec.logOR.obs.case.excess = vec.logOR.obs.case
  vec.logOR.obs.case.excess[!ind.both.case]=NA
  
  boxplot(mat.logOR.ppd.case.1,xaxt="n",at=at.loc-loc.gap,boxwex=boxwex,
          ylim=c(-4.5,4.5),outline=FALSE)
  boxplot(mat.logOR.ppd.case.2,xaxt="n",add=TRUE,at=at.loc+loc.gap,boxwex=boxwex,
          col="dodgerblue2",outline=FALSE)
  mtext("Log Odds Ratio",side = 2,at=-4,line=3,cex=2)
  
  abline(h=0,col="gray",lwd=2)
  
  points(1:choose(JBrS,2),vec.logOR.obs.case.excess,xaxt="n",xlab="pathogen pairs",
         ylab="log odds ratio",col="red",pch=1,cex=vec.case.size,lwd=2)
  points(1:choose(JBrS,2),vec.logOR.obs.case.excess,xaxt="n",xlab="pathogen pairs",
         ylab="log odds ratio",col="red",pch=20,cex=1)
  
  name_pairs <- names(vec.logOR.obs.case)
  name_pairs[!ind.axis.show] <- NA
  axis(1,at = 1:choose(JBrS,2),labels=name_pairs,las=2)
  legend("top","case",bty="n",cex=3)
  
  #
  # controls:
  #
  par(mar = c(4,6,0,2))
  
  vec.logOR.obs.ctrl.excess = vec.logOR.obs.ctrl
  vec.logOR.obs.ctrl.excess[!ind.both.ctrl]=NA
  
  boxplot(mat.logOR.ppd.ctrl.1,xaxt="n",at=at.loc-loc.gap,boxwex=boxwex,
          ylim=c(-4.5,4.5),outline=FALSE)
  boxplot(mat.logOR.ppd.ctrl.2,xaxt="n",add=TRUE,at=at.loc+loc.gap,boxwex=boxwex,
          col = "dodgerblue2",outline=FALSE)
  abline(h=0,col="gray",lwd=2)
  points(1:choose(JBrS,2),vec.logOR.obs.ctrl.excess,xaxt="n",xlab="pathogen pairs",
         ylab="log odds ratio",col="red",pch=1,cex=vec.ctrl.size,lwd=2)
  points(1:choose(JBrS,2),vec.logOR.obs.ctrl.excess,xaxt="n",xlab="pathogen pairs",
         ylab="log odds ratio",col="red",pch=20,cex=1)
  #axis(1,at = 1:choose(JBrS,2),labels=names(vec.logOR.obs.ctrl),las=2)
  
  legend("bottom","control",bty="n",cex=3)
  #   arrows(14, -2, 19-boxwex,0 , col = "black",length=0.1,angle=20,lwd=2)
  #   text(14,-2.5,"pLCM",cex=2)
  #   arrows(24, -2, 19+boxwex,0 , col = "dodgerblue2",length=0.1,angle=20,lwd=2)
  #   text(24,-2.5,"npLCM",cex=2,col="dodgerblue2")
  
  dev.off()
  
  cat("== A figure is generated for model checking: parallel pairwise log 
       odds ratio for two model fits.","\n",
       "The two model fits are extracted from:", "\n",
      "(1)",DIR_NPLCM1,"\n",
      "(2)",DIR_NPLCM2,"\n",
      "The figure is stored in the current working directory ==")

}
