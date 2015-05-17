#' Posterior predictive checking for the nested partially class models - 
#' frequent patterns in the BrS data.
#' 
#' At each MCMC iteration, we generate a new data set based on the model and 
#' parameter values at that iteration. The sample size of the new data set equals
#' that of the actual data set, i.e. the same number of cases and controls.
#' 
#' @param DIR_NPLCM File path to the folder that stores results from npLCM fit.
#' @param npat.ctrl Number of the most common BrS measurement pattern among controls.
#' Default is 10.
#' @param npat.case Number of the most common BrS measurement pattern among cases.
#' Default is 10.
#' @param cex.pattern Size of patterns; Default is 1.
#' 
#' @return A figure of posterior predicted frequencies compared with the observed 
#' frequencies of the most common patterns for the BrS data. The function generates
#' this figure in your working directory automatically.
#' @export
#' 

check_common_pattern <- function(DIR_NPLCM,npat.case=10,npat.ctrl = 10,
                                          cex.pattern = 1){
    
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
  
 
  #
  # can be modified to actual visualization order:
  #
  pathogen_BrS_list <- model_options$pathogen_BrS_list
  
  JBrS  <- bugs.dat$JBrS
  Nd    <- bugs.dat$Nd
  Nu    <- bugs.dat$Nu
  #
  # test: (pathogen_display can be separately specified):
  #
  pathogen_display  <- pathogen_BrS_list
  index_display <- my_reorder(pathogen_display,pathogen_BrS_list)
  pathogen_name <- pathogen_BrS_list[index_display]
  
  MBS.case <- bugs.dat$MBS[Y==1,index_display]
  MBS.ctrl <- bugs.dat$MBS[Y==0,index_display]
  
  # get posterior predicted datasets:
  MBS.new <- res_nplcm[,grep("MBS.new",colnames(res_nplcm))]
  new.tmp <- MBS.new
  
  #
  # Cases:
  #
  nonzero.pat <- apply(MBS.case,1,paste,collapse = "" )
  casepat     <- sort(table(nonzero.pat),decreasing=TRUE)
  
  npat.case <- min(npat.case,length(casepat))
  casepat.high      <- casepat[1:npat.case]/Nd
  casepat.high.name <- names(casepat.high)
  
  # Rows: iterations; 
  # Columns: frequencies of npat most common patterns+the rest:
  ctpat.case <- matrix(NA,nrow=nrow(new.tmp),ncol=npat.case+1)
  
  for (iter in 1:nrow(new.tmp)){
    mat             <- matrix(new.tmp[iter,1:(JBrS*Nd)],ncol=JBrS,byrow=TRUE)[,index_display]
    nonzero.pat.tmp <- apply(mat , 1 , paste , collapse = "" )
    casepat.tmp     <- table(nonzero.pat.tmp)
    
    ppd.pat.ct = (sapply(casepat.high.name,function(x) {
      indtmp = names(casepat.tmp)==x
      if (sum(indtmp)==0){
        res = 0
      } else{
        res = casepat.tmp[indtmp]
      }     
      res
    }
    ))
    
    ctpat.case[iter,-ncol(ctpat.case)] = ppd.pat.ct/Nd
    ctpat.case[iter,ncol(ctpat.case)]  = 1-sum(ppd.pat.ct/Nd)
  }
  
  #
  # Controls:
  #
  nonzero.pat <- apply(MBS.ctrl,1,paste,collapse = "" )
  ctrlpat     <- sort(table(nonzero.pat),decreasing=TRUE)
  
  npat.ctrl <- min(npat.ctrl,length(ctrlpat))
  ctrlpat.high      <- ctrlpat[1:npat.ctrl]/Nu
  ctrlpat.high.name <- names(ctrlpat.high)
  
  # Rows: iterations; 
  # Columns: frequencies of npat most common patterns+the rest:
  ctpat.ctrl <- matrix(NA,nrow=nrow(new.tmp),ncol=npat.ctrl+1)
  
  for (iter in 1:nrow(new.tmp)){
    mat             <- matrix(new.tmp[iter,-(1:(JBrS*Nd))],ncol=JBrS,byrow=TRUE)[,index_display]
    nonzero.pat.tmp <- apply(mat,1,paste,collapse = "" )
    ctrlpat.tmp     <- table(nonzero.pat.tmp)
    
    ppd.pat.ct = (sapply(ctrlpat.high.name,function(x) {
      indtmp = names(ctrlpat.tmp)==x
      if (sum(indtmp)==0){
        res = 0
      } else{
        res = ctrlpat.tmp[indtmp]
      }     
      res
    }
    ))
    
    ctpat.ctrl[iter,-ncol(ctpat.ctrl)] = ppd.pat.ct/Nu
    ctpat.ctrl[iter,ncol(ctpat.ctrl)]  =  1-sum(ppd.pat.ct/Nu)
  }
  
  
  #
  # start plotting:
  #
  #
  # start plotting:
  #
  pdf(paste0(DIR_NPLCM,"//frequent_pattern_check.pdf"),
      width=20,height=10)
  nf  <- layout(matrix(c(1,2),nrow=1,ncol=2),
                widths=c(10,10),heights=c(6,6))
  top <- max(ctpat.case,ctpat.ctrl)+0.05
  #layout.show(nf)
  
  #
  # cases:
  #
  boxwex =0.2 
  loc.gap = 0#boxwex/1.9
  par(mar=c(5.1, 4.1, 4.1, 2.1)+c(JBrS,1,0,-2.1),xpd=TRUE)
  boxplot(ctpat.case,at=1:(npat.case+1)-loc.gap,xlab="",
          ylab="frequency",cex.lab=2,xaxt="n",cex.main=2,ylim=c(0,top),
          boxwex=boxwex,yaxt="n",outline=FALSE,bty="n")
  
  casepat.high.name <- NA2dot(casepat.high.name)
  
  #mtext("pattern (ordered by observed frequency)",1,line=8)
  axis(1,at=1:(length(casepat.high.name)+1),labels=c(casepat.high.name,"other"),
       las=2,cex.axis=cex.pattern)
  axis(2,at = seq(0.1,top,by=0.1),labels=seq(0.1,top,by=0.1),cex.axis=2)
  ci95.case = apply(ctpat.case,2,quantile,c(0.025,0.975))
  matplot(1:(npat.case+1)-loc.gap,t(ci95.case),add=TRUE,col="black",pch="*",cex=2)
  
  points(c(casepat.high,1-sum(casepat.high)),col="red",pch=1,cex=1,lwd=2)
  cumpat = round(cumsum(c(casepat.high,1-sum(casepat.high)))*100,1)
  
  mtext("case",side = 3,line=-5,cex=4)
  mtext("BrS pattern",side=1,line=JBrS+3,cex=2)
  # # cumulative frequencies:
  # plot(1:(npat.case+1),rep(0,npat.case+1)+cumpat/100,type="o",ylim=c(0,1),
  #      ylab="Cumulative",xaxt="n",xlab="n")
  # for (s in 1:(length(casepat.high)+1)){
  #   text(s,top+0.07,paste0(cumpat[s],"%"),srt=45,cex=2)
  # }
  
  #
  # controls:
  #
  boxwex =0.2 
  loc.gap = 0#boxwex/1.9
  par(mar=c(5.1, 4.1, 4.1, 2.1)+c(JBrS,-4.1,0,0),xpd=TRUE)
  boxplot(ctpat.ctrl,at=1:(npat.ctrl+1)-loc.gap,xlab="",
          ylab="",xaxt="n",cex.main=2,ylim=c(0,top),
          boxwex=boxwex,yaxt="n",cex.lab=2,outline=FALSE,
          yaxt="n")
  
  ctrlpat.high.name <- NA2dot(ctrlpat.high.name)
  # mtext("pattern (ordered by observed frequency)",1,line=8)
  axis(1,at=1:(length(ctrlpat.high.name)+1),labels=c(ctrlpat.high.name,"other"),
       las=2,cex.axis=cex.pattern)
  #axis(2,at = seq(0.1,top,by=0.1),labels=seq(0.1,top,by=0.1),cex.axis=2)
  ci95.ctrl = apply(ctpat.ctrl,2,quantile,c(0.025,0.975))
  matplot(1:(npat.ctrl+1)-loc.gap,t(ci95.ctrl),add=TRUE,col="black",pch="*",cex=2)
  
  points(c(ctrlpat.high,1-sum(ctrlpat.high)),col="red",pch=1,cex=1,lwd=2)
  cumpat = round(cumsum(c(ctrlpat.high,1-sum(ctrlpat.high)))*100,1)
  
  mtext("control",side = 3,line=-5,cex=4)
  mtext("BrS pattern",side=1,line=JBrS+3,cex=2)
  legend("topright",c("observed frequency","97.5% posterior predictive quantile",
                      "2.5% posterior predictive quantile"),
         col=c("red","black","black"),pch=c(1,8,8),bty="n")
  
  # par(mar=c(0,0,0,2.1))
  # plot(1:(npat.ctrl+1),rep(0,npat.ctrl+1)+cumpat/100,type="o",ylim=c(0,1),
  #      xaxt="n",xlab="n",ylab="",yaxt="n")
  # for (s in 1:(length(ctrlpat.high)+1)){
  #   text(s,top+0.07,paste0(cumpat[s],"%"),srt=45,cex=2)
  # }
  dev.off()
  
  
  cat("==A figure is generated for model checking: frequent BrS measurement
      patterns. Stored in ",DIR_NPLCM," ==")
}