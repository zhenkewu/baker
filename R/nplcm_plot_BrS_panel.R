#' Plot bronze-standard (BrS) panel
#' 
#' Now only works for singleton etiologies.
#' 
#' @param MBS Matrix of bronze-standard measurements. Rows for subjects (cases at
#' the top, controls at the bottom), columns for pathogen-specimen combination.
#' @param model_options See \code{\link{nplcm}}
#' @param clean_options See \code{\link{clean_perch_data}}
#' @param res_nplcm See \code{\link{nplcm_read_folder}}
#' @param bugs.dat Data input for the model fitting.
#' @param top_BrS Numerical value to specify the rightmost limit 
#' on the horizontal axis for the BrS panel.
#' @param prior_shape \code{interval} or \code{boxplot} - for how to represent
#' prior/posteriors of the TPR/FPRs of measurements.
#' @param cexval Default is 1 - size of text of the BrS percentages.
#' @param srtval Default is 0 - the direction of the text for the BrS percentages.
#' 
#' @importFrom binom binom.confint
#' 
#' @export

nplcm_plot_BrS_panel <- function(MBS,model_options,clean_options,res_nplcm,
                                 bugs.dat,
                                 top_BrS = 1.3,prior_shape = "boxplot",
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
  
  
  colnames(MBS) <-model_options$pathogen_BrS_list
  if (is.null(clean_options$allow_missing)||
        clean_options$allow_missing==FALSE){
    tmp.case <- binom.confint(colSums(MBS[1:Nu,]),
                              Nd, conf.level = 0.95, methods = "ac")[ord,]
    tmp.ctrl <- binom.confint(colSums(MBS[-(1:Nd),]),
                              Nu, conf.level = 0.95, methods = "ac")[ord,]
  }else{
    ind_case_BrS_dat_not_na <- which(rowSums(is.na(MBS[1:Nd,]))==0)
    ind_ctrl_BrS_dat_not_na <- which(rowSums(is.na(MBS[-(1:Nd),]))==0)
    
    tmp.case <- binom.confint(colSums(MBS[(1:Nd)[ind_case_BrS_dat_not_na],]),
                              Nd, conf.level = 0.95, methods = "ac")[ord,]
    tmp.ctrl <- binom.confint(colSums(MBS[(Nd+(1:Nu))[ind_ctrl_BrS_dat_not_na],]),
                              Nu, conf.level = 0.95, methods = "ac")[ord,]
  }

  # case and control positive rate, lower and upper limit
  Bcomp    <- rbind(round(tmp.case$mean,5),round(tmp.ctrl$mean,5))
  Bcomp_q1 <- rbind(tmp.case[,c("lower")],tmp.ctrl[,c("lower")])
  Bcomp_q2 <- rbind(tmp.case[,c("upper")],tmp.ctrl[,c("upper")])
  
  
  # sensitivity (TPR) and 1-specificity (FPR) information
  if (model_options$k_subclass==1){
    #
    # plcm (just obtain TPR and FPR estimates):
    #
    if (is.null(model_options$SSonly)|| model_options$SSonly==FALSE){
      theta_mat <- (res_nplcm[,grep("thetaBS",colnames(res_nplcm))])[,ord]
    }else {
      theta_mat <- cbind(res_nplcm[,grep("thetaBS",colnames(res_nplcm))],
                         matrix(NA,nrow=nrow(res_nplcm),ncol=Jfull-JBrS))[,ord]
    }
    theta_mean <- colMeans(theta_mat)
    
    #posterior distribution of FPR:
    if (is.null(model_options$SSonly)||model_options$SSonly==FALSE){
      psi_mat  <- (res_nplcm[,grep("psiBS",colnames(res_nplcm))])[,ord]
    }else{
      psi_mat  <- cbind(res_nplcm[,grep("psiBS",colnames(res_nplcm))],
                        matrix(NA,nrow=nrow(res_nplcm),ncol=Jfull-JBrS))[,ord]
    }
    psi_mean   <- colMeans(psi_mat)
    
    ## model fitted postive rate for each pathogen
    fittedmean_case <- sapply(1:Jcause,
                              function(s) mean(pEti_mat_ord[,s]*theta_mat[,s]+
                                                 (1-pEti_mat_ord[,s])*psi_mat[,s]))
    fittedmean_control <- psi_mean
  
  } else {
    #
    # nplcm (need to marginalize over subclasses to calculate marginal TPR):
    #
    if (is.null(model_options$SSonly)|| model_options$SSonly==FALSE){
      theta_mat  <- (res_nplcm[,grep("ThetaBS.marg",colnames(res_nplcm))])[,ord]
    }else{
      theta_mat <- cbind(res_nplcm[,grep("ThetaBS.marg",colnames(res_nplcm))],
                         matrix(NA,nrow=nrow(res_nplcm),ncol=Jfull-JBrS))[,ord]
    }
    theta_mean <- colMeans(theta_mat)
    
    #posterior distribution of FPR
    if (is.null(model_options$SSonly)||model_options$SSonly==FALSE){
       psi_mat  <- (res_nplcm[,grep("PsiBS.marg",colnames(res_nplcm))])[,ord]
    }else{
       psi_mat  <- cbind(res_nplcm[,grep("PsiBS.marg",colnames(res_nplcm))],
                        matrix(NA,nrow=nrow(res_nplcm),ncol=Jfull-JBrS))[,ord]
    }
    psi_mean <- colMeans(psi_mat)
    
    psi_mat.case = array(NA,c(nrow(res_nplcm),Jfull,Jcause),
                         dimnames = list(NULL,1:Jfull,1:Jcause))
    
    psi_mat.case.tmp = (res_nplcm[,grep("PsiBS.case",colnames(res_nplcm))])
    
    for (j in 1:JBrS){
      for (s in 1:Jcause){
        tmp.nm  <- paste0("PsiBS.case","[",j,",",s,"]")
        ind.tmp <- which(colnames(psi_mat.case.tmp)==tmp.nm)
        psi_mat.case[,j,s] = psi_mat.case.tmp[,ind.tmp]
      }
    }
    
    psi_mat.case.ord = psi_mat.case[,ord,ord]
    
    # model fitted postive rate for each pathogen
    fittedmean_case    <- sapply(1:Jfull,function(s) 
                            mean(pEti_mat_ord[,s]*theta_mat[,s]+
                              rowSums(pEti_mat_ord[,-s]*psi_mat.case.ord[,s,-s])))
    fittedmean_control <- psi_mean
    
  }
  
  #
  # start plotting:
  #
  
  op <- par(mar=c(5.1,4.1,4.1,0))
  
  plotat <- seq(0.5,Jfull+0.5,by=1/4)[-(c(1,(1:Jfull)*4+1))]
  #plot case control positive rates, and fitted case rates
  plot(c(rbind(fittedmean_case,Bcomp)),plotat,yaxt="n",
       xlim=c(0,top_BrS),
       ylim=c(0.5,Jfull+.5),xaxt="n",
       ylab="",xlab="probability",
       pch = c(rbind(rep(2,Jfull),rep(20,Jfull),rep(20,Jfull))),
       col=c(rbind(rbind(rep(1,Jfull),rep("dodgerblue2",Jfull),rep("dodgerblue2",Jfull)))),
       cex = c(rbind(rep(1,Jfull),rep(2,Jfull),rep(2,Jfull))))
  
  #add axis labels on the left:
  axis(2,at = plotat,labels=rep(c("","case","ctrl"),Jfull),las=2)
  #add ticks from 0 to 1 for x-bar:
  axis(1,at = c(0,0.2,0.4,0.6,0.8,1),labels= c(0,0.2,0.4,0.6,0.8,1),las=1)
  
  #plot TPR posterior mean, and upper CI bound for case/control rates:
  points(c(rbind(theta_mean,Bcomp_q2)),plotat,
         pch=c(rbind(rep("+",Jfull),rep("|",Jfull),rep("|",Jfull))),
         cex=c(rbind(rep(2,Jfull),rep(1,Jfull),rep(1,Jfull))),
         col=c(rbind(rep(1,Jfull),rep(1,Jfull),rep(1,Jfull))))
  #plot FPR posterior mean, and lower CI bound for case/control rates:
  points(c(rbind(psi_mean,Bcomp_q1)),plotat,
         pch=c(rbind(rep("*",Jfull),rep("|",Jfull),rep("|",Jfull))),
         cex=c(rbind(rep(2,Jfull),rep(1,Jfull),rep(1,Jfull))),
         col=c(rbind(rep(1,Jfull),rep(1,Jfull),rep(1,Jfull))))
  abline(h=seq(1.5,Jfull-.5,by=1),lty=2,lwd=0.5,col="blue")
  abline(v=1,lty=2,lwd=.5)
  
  ## conditional odds ratios
  COR = function(brs.data,nd,pathogens){
    y = rep(c(1,0),times=c(nd,nrow(brs.data)-nd))
    #X = matrix(NA,nrow=nrow(brs.data),ncol=length(pathogens))
    x.nm   = paste(pathogens)
    #colnames(X) = x.nm
    #for (j in 1:length(pathogens)){
    #  X[,j] = brs.data[,x.nm[j]]
    #}
    dat.reg = as.data.frame(cbind(y,brs.data))
    form = as.formula(paste0("y~",paste(x.nm,collapse="+")))
    
    fit = glm(form,data=dat.reg,family=binomial,na.action="na.omit")
    
    if (sum(is.na(coef(fit)))==0 & sum(diag(vcov(fit))>100)==0){
      res0 = cbind(exp(suppressMessages(confint(fit))),exp(fit$coef))[-1,]
      res = list(ORinterval = res0,label = "conditional OR")
    }else{
      print("Conditional OR not calculatable. Switch to mariginal OR.")
      res0 = matrix(NA,nrow=length(pathogens),ncol=3)
      res0 = data.frame(res0)
      for (l in 1:nrow(res0)){
        tb <- table(dat.reg$y,dat.reg[,x.nm[l]])
        form_tmp <- as.formula(paste0("y~",x.nm[l]))
        fit_tmp  <- glm(form_tmp,data=dat.reg,family=binomial,na.action = "na.omit")
        if (length(vcov(fit_tmp))>1 && vcov(fit_tmp)[2,2]<100 && ncol(tb)==2){
          #print(l)
          res0[l,] <- cbind(exp(suppressMessages(confint(fit_tmp))),exp(fit_tmp$coef))[-1,]
        }
      }
      res = list(ORinterval=res0,label="marginal OR")
    }
  }
  
  tmp0 <- COR(MBS,Nd,pathogens_ord[ord<=JBrS])
  tmp  <- tmp0$ORinterval
  
  #plot conditional odds ratio on the right:
  incre = 0
  for (s in (1:Jfull)){
    if (s %in% which(ord<=JBrS)){
      incre = incre + 1
      L=round(tmp[incre,1],1)
      C=round(tmp[incre,3],1)
      R=round(tmp[incre,2],1)
      text(top_BrS-0.12,s+1/(2*Jfull),C,cex=1.5)
      text(top_BrS-0.12,s-.2,paste(c(L,"   ",R),collapse=" "),cex=1.2)
    }
  }
  legend("topright",tmp0$label,bty="n")
  
  counter = 0
  #each row, connect FPR and TPR posterior means
  for (s in 1:(3*Jfull)){
    segments(y0=plotat[s],x0=c(rbind(psi_mean,Bcomp_q1))[s],
             y1=plotat[s],x1=c(rbind(theta_mean,Bcomp_q2))[s],
             lty=ifelse((s-1)%%3<1,4,1))
    if ((s-1)%%3>=1){
      counter=counter+1
      tmp.hpos = ifelse(c(Bcomp_q2)[counter]+0.15>0.95,c(Bcomp_q1)[counter]-0.2,c(Bcomp_q2)[counter]+0.15 )
      text(tmp.hpos,plotat[s],paste0(round(100*c(Bcomp),1)[counter],"%"),
           srt=srtval,cex=cexval)
    }
  }
  
  for (s in 1:(Jfull)){
    segments(y0=plotat[3*s-1],x0=c(rbind(fittedmean_case,Bcomp))[3*s-1],
             y1=plotat[3*s],x1=c(rbind(fittedmean_case,Bcomp))[3*s],col="dodgerblue2",lwd=2)
  }
  # put prior shapes on bronze-standard sensitivity
  
  alphaB         <- bugs.dat$alphaB
  betaB          <- bugs.dat$betaB
  
  for (s in 1:Jfull){
    if (s %in% which(ord<=JBrS)){
      if (prior_shape == "interval"){
          # prior of TPR:
          tmp = rbeta(10000,alphaB[ord[s]],betaB[ord[s]])
          
          points(quantile(tmp,0.025),s-.45,pch="|")
          points(quantile(tmp,0.975),s-.45,pch="|")
          points(quantile(tmp,0.25),s-.45,pch="[")
          points(quantile(tmp,0.75),s-.45,pch="]")
          segments(quantile(tmp,0.025),s-.45,
                   quantile(tmp,0.975),s-.45,lty=1)
          
          # posterior of TPR:
          tmp.post = as.matrix(theta_mat)[,s]
          points(quantile(tmp.post,0.025),s-.35,pch="|")
          points(quantile(tmp.post,0.975),s-.35,pch="|")
          points(quantile(tmp.post,0.25),s-.35,pch="[")
          points(quantile(tmp.post,0.75),s-.35,pch="]")
          segments(quantile(tmp.post,0.025),s-.35,
                   quantile(tmp.post,0.975),s-.35,lty=1)
      } else if (prior_shape == "boxplot"){
          tmp = rbeta(10000,alphaB[ord[s]],betaB[ord[s]])
          boxplot(tmp,at = s-0.45, boxwex=1/10 , col="gray",
                  add=TRUE,horizontal=TRUE,outline=FALSE,xaxt="n")
          tmp.post = as.matrix(theta_mat)[,s]
          boxplot(tmp.post,at = s-0.35,boxwex=1/10,add=TRUE,
                  horizontal=TRUE,outline=FALSE,xaxt="n")
      }
    }
  }
  axis(2,at=(1:Jfull)-.45,labels=rep("",Jfull),las=2,cex.axis=.5)
  axis(2,at=(1:Jfull)-.35,labels=rep("",Jfull),las=2,cex.axis=.5)
  
  mtext(expression(underline("BrS")),line=1,cex=1.8)

    
}