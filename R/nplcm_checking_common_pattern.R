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
#' @importFrom coda read.coda
#' @return A figure of posterior predicted frequencies compared with the observed 
#' frequencies of the most common patterns for the BrS data. The function generates
#' this figure in your working directory automatically.
#' @export
#' 

nplcm_checking_common_pattern <- function(DIR_NPLCM,npat.ctrl = 10,npat.case=10){
  
#   DIR_NPLCM = result.folder
#   npat.ctrl = 20
#   npat.case = 20
  
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
  
  ## reading nplcm outputs:
  res_nplcm <- read.coda(paste(DIR_NPLCM,"coda1.txt",sep="/"),
                         paste(DIR_NPLCM,"codaIndex.txt",sep="/"),
                         quiet=TRUE)
  
  JBrS  <- bugs.dat$JBrS
  Nd <- bugs.dat$Nd
  Nu <- bugs.dat$Nu
  Y  <- c(rep(1,Nd),rep(0,Nu))
  
  MBS.case <- bugs.dat$MBS[Y==1,]
  MBS.ctrl <- bugs.dat$MBS[Y==0,]
  
  #
  # Controls:
  #
  nonzero.pat <- apply(MBS.ctrl,1,paste,collapse = "" )
  ctrlpat     <- sort(table(nonzero.pat),decreasing=TRUE)
  
  ctrlpat.high      <- ctrlpat[1:npat.ctrl]/Nu
  ctrlpat.high.name <- names(ctrlpat.high)
  
  # get posterior predicted datasets:
  MBS.new <- res_nplcm[,grep("MBS.new",colnames(res_nplcm))]
  new.tmp <- MBS.new
  
  # Rows: iterations; 
  # Columns: frequencies of npat most common patterns+the rest:
  ctpat.ctrl <- matrix(NA,nrow=nrow(new.tmp),ncol=npat.ctrl+1)
  
  for (iter in 1:nrow(new.tmp)){
    mat             <- matrix(new.tmp[iter,-(1:(JBrS*Nd))],ncol=JBrS,byrow=TRUE)
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
  # Cases:
  #
  nonzero.pat <- apply(MBS.case,1,paste,collapse = "" )
  casepat     <- sort(table(nonzero.pat),decreasing=TRUE)
  
  casepat.high      <- casepat[1:npat.case]/Nd
  casepat.high.name <- names(casepat.high)
  
  ctpat.case <- matrix(NA,nrow=nrow(new.tmp),ncol=npat.case+1)
  
  for (iter in 1:nrow(new.tmp)){
    mat             <- matrix(new.tmp[iter,1:(JBrS*Nd)],ncol=JBrS,byrow=TRUE)
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
  # start plotting:
  #
  op<-par()
  pdf(paste0(result.folder,"//frequent_pattern_check.pdf"),
      width=16,height=16)
  nf  <- layout(matrix(c(1,2),nrow=2),widths=c(20,20),heights=c(8,8))
  top <- .5
  #layout.show(nf)
  ## controls:
  boxwex =0.2 
  loc.gap = boxwex/1.9
  par(mai=op$mai+c(JBrS/11+.5,1,0,0),xpd=TRUE)
  boxplot(ctpat.ctrl,at=1:(npat.ctrl+1)-loc.gap,xlab="",
          ylab="frequency",xaxt="n",cex.main=2,ylim=c(0,top),
          boxwex=boxwex,yaxt="n",cex.lab=2,outline=FALSE)
  
  # mtext("pattern (ordered by observed frequency)",1,line=8)
  axis(1,at=1:(length(ctrlpat.high.name)+1),labels=c(ctrlpat.high.name,"other"),
       las=2,cex.axis=2)
  axis(2,at = seq(0.1,top,by=0.1),labels=seq(0.1,top,by=0.1),cex.axis=2)
  ci95.ctrl = apply(ctpat.ctrl,2,quantile,c(0.025,0.975))
  matplot(1:(npat.ctrl+1)-loc.gap,t(ci95.ctrl),add=TRUE,col="black",pch="*",cex=2)
  
  points(c(ctrlpat.high,1-sum(ctrlpat.high)),col="red",pch=1,cex=1,lwd=2)
  cumpat = round(cumsum(c(ctrlpat.high,1-sum(ctrlpat.high)))*100,1)
  for (s in 1:(length(ctrlpat.high)+1)){
    text(s,top+0.07,paste0(cumpat[s],"%"),srt=45,cex=2)
  }
  #legend("topright",c("observed frequency","2.5% posterior predictive quantile",
  #                    "97.5% posterior predictive quantile"),
  #col=c("red","blue","blue"),lty=c(1,1,2),pch=c(15,3,3))
  
  ## cases:
  boxwex =0.2 
  loc.gap = boxwex/1.9
  par(mai=op$mai+c(JBrS/11+.5,1,0,0),xpd=TRUE)
  boxplot(ctpat.case,at=1:(npat.case+1)-loc.gap,xlab="",
          ylab="frequency",cex.lab=2,xaxt="n",cex.main=2,ylim=c(0,top),
          boxwex=boxwex,yaxt="n",outline=FALSE)
  
  #mtext("pattern (ordered by observed frequency)",1,line=8)
  axis(1,at=1:(length(casepat.high.name)+1),labels=c(casepat.high.name,"other"),
       las=2,cex.axis=2)
  axis(2,at = seq(0.1,top,by=0.1),labels=seq(0.1,top,by=0.1),cex.axis=2)
  ci95.case = apply(ctpat.case,2,quantile,c(0.025,0.975))
  matplot(1:(npat.case+1)-loc.gap,t(ci95.case),add=TRUE,col="black",pch="*",cex=2)
  
  points(c(casepat.high,1-sum(casepat.high)),col="red",pch=1,cex=1,lwd=2)
  cumpat = round(cumsum(c(casepat.high,1-sum(casepat.high)))*100,1)
  for (s in 1:(length(casepat.high)+1)){
    text(s,top+0.07,paste0(cumpat[s],"%"),srt=45,cex=2)
  }
  dev.off()
  
  cat("==A figure is generated for model checking: frequent BrS measurement
      patterns. Stored in ",DIR_NPLCM," ==")
}