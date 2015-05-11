#' Plot etiology (pie) panel
#' 
#' Now only works for singleton etiologies.
#' 
#' @param model_options See \code{\link{nplcm}}
#' @param res_nplcm See \code{\link{nplcm_read_folder}}
#' @param bugs.dat Data input for the model fitting.
#' @param top_pie Numerical value to specify the rightmost limit 
#' on the horizontal axis for the pie panel.
#' 
#' @importFrom binom binom.confint
#' 
#' @export

nplcm_plot_pie_panel <- function(model_options,res_nplcm,
                                 bugs.dat,
                                 top_pie = 1){
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
  # plot etiology information:
  #
  dotcolor = "black"
  #op <- par(mar=c(5.1,6,4.1,1.1))
  op <- par(mar=c(5.1,0,4.1,10))
  plot(c(pEti_mean_ord),1:(Jfull),
       yaxt="n",#xaxt="n",
       xlim=c(0,top_pie),ylim=c(0.5,Jfull+0.5),col=c("black"),
       ylab="",xlab="probability",
       pch=c(20),cex=2)
  axis(4,at=1:Jfull,labels=paste(paste(pathogens_ord,ord,sep=" ("),")",sep=""),las=2,cex.axis=1.5)
  abline(h=seq(1.5,Jfull-.5,by=1),lty=2,lwd=0.5,col="blue")
  #draw axis within plot:
  for (s in 1:(Jfull-1)){
    axis(1, seq(0,1,by=.2), lwd=0,lwd.ticks=1,#labels=rep("",length(seq(0,1,by=.2))),
         pos = seq(1.5,Jfull-.5,by=1)[s], cex.axis = 0.8,lty=2,col="blue")
    # axis(1, seq(0,1,by=.2), lwd=0,lwd.ticks=0,#labels=rep("",length(seq(0,1,by=.2))),
    #      pos = seq(1.5,Jfull-.5,by=1)[s]+0.3, cex.axis = 0.8,lty=2,col="blue")
  }
  points(c(pEti_q1),1:(Jfull),pch="|",cex=1)
  points(c(pEti_q2),1:(Jfull),pch="|",cex=1)
  points(c(pEti_innerq1),1:(Jfull),pch="[",cex=1)
  points(c(pEti_innerq2),1:(Jfull),pch="]",cex=1)
  
  mtext(expression(underline(hat(pi))),line=1,cex=1.8)
  #mtext(c(expression(bold("--")),":prior","-",":posterior"),col=c("gray","black","black","black"),
  #      adj=c(0,0.1,0.3,0.4),line=.8,cex=.8,lwd=2)
  legend("topright",c("prior","posterior"),lty=c(2,1),col=c("gray","black"),
         lwd = 4,horiz=TRUE,cex=1,bty="n")
  pgrid = seq(0,1,by=0.01)
  
  alpha <- bugs.dat$alpha#eti_prior_set(model_options)
  
  for (s in 1:(Jfull)){
    segments(y0=s,x0=c(pEti_q1)[s],y1=s,x1=c(pEti_q2)[s],col=dotcolor)
    segments(y0=s,x0=c(pEti_innerq1)[s],y1=s,x1=c(pEti_innerq2)[s],col = dotcolor,lwd=2)
    text(.8,s,paste0("=",paste0(round(100*c(pEti_mean_ord),1)[s],"%")),srt=0,cex=2)
    text(.65,s,bquote(hat(pi)[.(ord[s])]),srt=0,cex=2)
    tmp.density = dbeta(pgrid,alpha[ord[s]],sum(alpha[-ord[s]]))
    points(pgrid,tmp.density/(3*max(tmp.density))+s-0.45,type="l",col="gray",lwd=4,lty=2)
    ##posterior density
    tmp.post.density = density(pEti_mat_ord[,s],from=0,to=1)
    tmp.x = tmp.post.density$x
    tmp.y = tmp.post.density$y
    points(tmp.x,tmp.y/(3*max(tmp.y))+s-0.45,col="black",lwd=4,type="l")
    
  }
  
}