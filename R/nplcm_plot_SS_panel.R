#' Plot silver-standard (SS) panel
#' 
#' Now only works for singleton etiologies. Current the prior shape can only be
#' represented by intervals.
#' 
#' @param MSS Matrix of silver-standard measurements. Rows for subjects (cases at
#' the top, controls at the bottom), columns for pathogen-specimen combination.
#' \code{MSS} has fewer columns than MBS and those pathogens with both BrS and SS
#' measurements should be arranged in the first several columns. This should
#' have been done by \code{\link{perch_data_clean}}.
#' @param model_options See \code{\link{nplcm}}
#' @param clean_options See \code{\link{perch_data_clean}}
#' @param res_nplcm See \code{\link{nplcm_read_folder}}
#' @param bugs.dat Data input for the model fitting.
#' @param top_SS Default is \code{0.3}. Numerical value to specify the rightmost limit 
#' on the horizontal axis for the SS panel.
#' @param cexval Default is 1 - size of text of the BrS percentages.
#' @param srtval Default is 0 - the direction of the text for the BrS percentages.
#' 
#' @importFrom binom binom.confint
#' 
#' @export

nplcm_plot_SS_panel <- function(MSS,model_options,clean_options,res_nplcm,
                                bugs.dat,
                                 top_SS = 0.3,
                                 cexval = 1,
                                 srtval = 0){
  #
  # now only deal with singleton etiologies:
  # 
  
  # total no. of causes:
  Jcause     <- length(model_options$cause_list)
  # extract and process some data and posterior samples:
  SubVarName <- rep(NA,Jcause)
  for (j in 1:Jcause){
    SubVarName[j] = paste("pEti","[",j,"]",sep="")
  }
  
  # get etiology fraction MCMC samples:
  pEti_mat   <- res_nplcm[,SubVarName]
  pEti_mean  <- colMeans(pEti_mat)
  pEti_mean0 <- pEti_mean
  
  # order the causes by posterior mean:
  ord <- order(pEti_mean)
  
  pEti_mean_ord <- pEti_mean[ord]
  pEti_mat_ord  <- pEti_mat[,ord]
  
  # quantiles for etiology: outer is 97.5% CI, inner is 50% CI
  pEti_q1   <- apply(pEti_mat,2,quantile,probs=0.025)[ord]
  pEti_q2   <- apply(pEti_mat,2,quantile,probs=0.975)[ord]
  pEti_innerq1   <- apply(pEti_mat,2,quantile,probs=0.25)[ord]
  pEti_innerq2   <- apply(pEti_mat,2,quantile,probs=0.75)[ord]
  
  # complete list of pathogens:
  if (is.null(model_options$SSonly) || model_options$SSonly==FALSE){
    pathogen_list     <- model_options$pathogen_BrS_list
  } else{
    pathogen_list     <- c(model_options$pathogen_BrS_list,
                           model_options$pathogen_SSonly_list)
  }
  
  Jfull               <- length(pathogen_list)
  
  if (Jfull != Jcause){
    stop("== The number of causes is different from the total number of pathogens.
         The multiple-cause visualization, including NoA, is being developed.
         Please contact developer. Thanks. ==")
  }
  JBrS                <- length(model_options$pathogen_BrS_list)
  pathogens_ord       <- pathogen_list[ord]
  
  Nd <- bugs.dat$Nd
  Nu <- bugs.dat$Nu
  
  #
  # prepare information for SS panel:
  #
  SS_index  <- which(colMeans(is.na(MSS))<.9)
  JSS       <- length(SS_index)
  
  ind.SS = rep(NA,JSS) # tells where the the SS row should go.
  for (j in 1:JSS){
    ind.SS[j] = which(ord==j)
  }
  if (!is.null(model_options$SSonly) && model_options$SSonly==TRUE){
    SSonlydat = bugs.dat$MSS.only
    JSS_only   = Jfull-JBrS
    ind.SSonly = rep(NA,JSS_only)
    for (j in 1:(JSS_only)){
      ind.SSonly[j] = which(ord == j+JBrS)
    }
    MSS = cbind(bugs.dat$MSS,SSonlydat)
    SS_index = which(colMeans(is.na(MSS))<.9)
    JSS      = length(SS_index)
    
    ind.SS = rep(NA,JSS)
    for (j in 1:JSS){
      ind.SS[j] = which(ord==ifelse(j<=JSS-JSS_only,j,
                                    j-JSS+JSS_only+JBrS))
    }
  }
  
  if (is.null(clean_options$allow_missing)||
        clean_options$allow_missing==FALSE){
    tmpSS.case = binom.confint(colSums(MSS[,1:JSS]), nrow(MSS),
                               conf.level = 0.95, methods = "ac")
  }else{
    ind_MSS_not_na <- which(rowSums(is.na(MSS[,1:JSS]))==0)
    tmpSS.case = binom.confint(colSums(MSS[ind_MSS_not_na,1:JSS]),
                               length(ind_MSS_not_na),
                               conf.level = 0.95, methods = "ac")
  }
  
  SScomp = rbind(round(tmpSS.case$mean,5),rep(NA,JSS))
  SScomp_q1 = rbind(tmpSS.case[,c("lower")],rep(NA,JSS))
  SScomp_q2 = rbind(tmpSS.case[,c("upper")],rep(NA,JSS))
  
  
  theta_matSS = (res_nplcm[,grep("thetaSS",colnames(res_nplcm))])
  theta_meanSS = colMeans(theta_matSS)
  
  theta_matSS_q1=apply(theta_matSS,2,quantile,0.025)
  theta_matSS_q2=apply(theta_matSS,2,quantile,0.975)
  
  fittedmean_SS_pos = sapply(1:JSS, function(s)
    mean(pEti_mat_ord[,ind.SS[s]]*theta_matSS[,s]))
  # note that here we used ind.SS[] to map back to the compelte
  # vector of pathogens.
  
  #
  # plotting SS panel:
  #
  par(mar=c(5.1,0,4.1,0))
  
  plotat = seq(0.5,Jfull+0.5,by=1/4)[-(c(1,(1:Jfull)*4+1))]
  #plotat.short = plotat[1:length(c(rbind(thetameanG,Gcomp)))]
  plotat.calc = function(j) {c(3*j-2,3*j-1,3*j)}
  plotat.short = rep(NA,JSS*3)
  for (j in 1:JSS){
    plotat.short[c(3*j-2,3*j-1,3*j)] = plotat[plotat.calc(ind.SS[j])]
  }
 
  
  # collect prior information on SS TPRs:
  alphaS <- bugs.dat$alphaS
  betaS  <- bugs.dat$betaS
  if (!is.null(model_options$SSonly) && model_options$SSonly==TRUE){
    alphaS.only <- bugs.dat$alphaS.only
    betaS.only  <- bugs.dat$betaS.only
    alphaS      <- c(bugs.dat$alphaS,alphaS.only)
    betaS       <- c(bugs.dat$betaS,betaS.only)
  }

  #
  # plotting:
  #
  plot(c(rbind(theta_meanSS,SScomp)),plotat.short,yaxt="n",xlim=c(0,top_SS),
       ylim=c(0.5,Jfull+.5),#xaxt="n",
       ylab="",xlab="probability",
       pch = c(rbind(rep(20,Jfull),rep(20,Jfull),rep(20,Jfull))),
       col=c(rbind(rbind(rep(1,Jfull),rep("blue",Jfull),rep(1,Jfull)))),
       cex = c(rbind(rep(1,Jfull),rep(2,Jfull),rep(2,Jfull))))
  
  points(c(rbind(fittedmean_SS_pos,matrix("",nrow=2,ncol=JSS))),plotat.short,
         yaxt="n",xlim=c(0,top_SS),
         ylim=c(0.5,Jfull+.5),xaxt="n",
         ylab="",#xlab="Gold Positive Rate",
         pch = c(rbind(rep(2,Jfull),rep(NA,Jfull),rep(NA,Jfull))),
         col=c(rbind(rbind(rep(1,Jfull),rep(1,Jfull),rep(1,Jfull)))),
         cex = c(rbind(rep(1,Jfull),rep(2,Jfull),rep(2,Jfull))))
  
  points(c(rbind(theta_matSS_q2,SScomp_q2)),plotat.short,
         pch=c(rbind(rep("|",Jfull),rep("|",Jfull),rep("|",Jfull))),
         cex=c(rbind(rep(1,Jfull),rep(1,Jfull),rep(1,Jfull))),
         col=c(rbind(rep(1,Jfull),rep(1,Jfull),rep(1,Jfull))))
  points(c(rbind(theta_matSS_q1,SScomp_q1)),plotat.short,
         pch=c(rbind(rep("|",Jfull),rep("|",Jfull),rep("|",Jfull))),
         cex=c(rbind(rep(1,Jfull),rep(1,Jfull),rep(1,Jfull))),
         col=c(rbind(rep(1,Jfull),rep(1,Jfull),rep(1,Jfull))))
  
  #inner 25%-75%
  theta_matSS_innerq1=apply(theta_matSS,2,quantile,0.25)
  theta_matSS_innerq2=apply(theta_matSS,2,quantile,0.75)
  points(c(rbind(theta_matSS_innerq1,SScomp_q1)),plotat.short,
         pch=c(rbind(rep("[",Jfull),rep("|",Jfull),rep("|",Jfull))),
         cex=c(rbind(rep(1,Jfull),rep(1,Jfull),rep(1,Jfull))),
         col=c(rbind(rep(1,Jfull),rep(1,Jfull),rep(1,Jfull))))
  points(c(rbind(theta_matSS_innerq2,SScomp_q1)),plotat.short,
         pch=c(rbind(rep("]",Jfull),rep("|",Jfull),rep("|",Jfull))),
         cex=c(rbind(rep(1,Jfull),rep(1,Jfull),rep(1,Jfull))),
         col=c(rbind(rep(1,Jfull),rep(1,Jfull),rep(1,Jfull))))
  counter = 0
  for (s in 1:length(plotat.short)){
    segments(y0=plotat.short[s],x0=c(rbind(theta_matSS_innerq1,SScomp_q1))[s],
             y1=plotat.short[s],x1=c(rbind(theta_matSS_innerq2,SScomp_q2))[s],
             col="black",
             lwd=1)
  }
  
  # row separation lines
  abline(h=seq(1.5,Jfull-.5,by=1)[ind.SS],lty=2,lwd=0.5,col="blue")
  abline(h=seq(1.5,Jfull-.5,by=1)[ind.SS]-1,lty=2,lwd=0.5,col="blue")
  
  
  counter = 0
  for (s in 1:length(plotat.short)){
    segments(y0=plotat.short[s],x0=c(rbind(theta_matSS_q1,SScomp_q1))[s],
             y1=plotat.short[s],x1=c(rbind(theta_matSS_q2,SScomp_q2))[s],col="black",
             lty=ifelse((s-1)%%3<2,1,1))
    if ((s-1)%%3>=1){
      counter=counter+1
      text(c(SScomp)[counter],plotat.short[s]+0.125,
           paste0(round(100*c(SScomp),1)[counter],"%"),srt=srtval,cex=cexval)
    }
  }
  
  for (s in 1:JSS){
    text(theta_meanSS[s],plotat.short[3*s-2]+.125,paste(round(100*theta_meanSS[s],2),"%"))
    
    # put prior shapes on gold sensitivity
    tmp = rbeta(10000,alphaS[s],betaS[s])
    points(quantile(tmp,0.025),ind.SS[s]-.45,pch="|")
    points(quantile(tmp,0.975),ind.SS[s]-.45,pch="|")
    points(quantile(tmp,0.25),ind.SS[s]-.45,pch="[")
    points(quantile(tmp,0.75),ind.SS[s]-.45,pch="]")
    segments(quantile(tmp,0.025),ind.SS[s]-.45,
             quantile(tmp,0.975),ind.SS[s]-.45,lty=1)
    
    #boxplot(tmp,at = ind.SS[s]-0.45, boxwex=1/8 ,col="gray",
    #        add=TRUE,horizontal=TRUE,outline=FALSE,xaxt="n")
    
  }
  
  mtext(expression(underline("SS")),line=1,cex=1.8)

}
