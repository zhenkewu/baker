#' Individual prediction from nested partially-latent class models
#'
#' DN:1.currently only works for both BrS and SS data, and SS only data. Any causes.
#'
#' @param DIR_NPLCM The file path to the folder that stores results from npLCM fit.
#' @param npat The number of most frequent patterns among cases to display.
#' @importFrom coda read.coda
#' @return A figure with individual diagnosis plots. The Bronze-Standard (BrS)
#'  measurement pattern is listed at the top of the barplots. Below it is the 
#'  percentage of the observed frequency of that pattern among cases.
#'  A second vector of 1/0 pattern is shown in red. It represents the 
#' Silver-Standard(SS) measurements. The percentage below it is the fraction of
#'  cases having this SS pattern among cases with the BrS pattern on the top.
#'  The whole list of causes are shown at the bottom.
#'
#' @export


nplcm_plot_individual_diagnosis <- function(DIR_NPLCM,npat=16){#BEGIN of function:

  # remember that the data.txt file in winbugs working folder is transposed:
  bugs.dat <- dget(paste(DIR_NPLCM,"data.txt",sep="/"))
  for (bugs.variable.name in names(bugs.dat)) {
    if (!is.null(dim(bugs.dat[[bugs.variable.name]]))) {
      dim(bugs.dat[[bugs.variable.name]]) <- rev(dim(bugs.dat[[bugs.variable.name]]))
      bugs.dat[[bugs.variable.name]] <- aperm(bugs.dat[[bugs.variable.name]])
    }
    assign(bugs.variable.name, bugs.dat[[bugs.variable.name]])
  }
  # change to model_options
  model_options  <- dget(paste(DIR_NPLCM,"model_options.txt",sep="/"))
  pathogen_list <- c(model_options$pathogen_BrS_list,
                     model_options$pathogen_SSonly_list)
  
  VarName ="Icat"
  
  MBS   = bugs.dat$MBS
  MSS   = bugs.dat$MSS
  JBrS  = bugs.dat$JBrS
  JSSonly = bugs.dat$JSSonly
  MSS.only = bugs.dat$MSS.only
  
  Jcasue = bugs.dat$Jcause
  
  Nd = bugs.dat$Nd
  Nu = bugs.dat$Nu
  Y = c(rep(1,Nd),rep(0,Nu))
  
  MBS.case = bugs.dat$MBS[Y==1,]
  SS_index  <- which(colMeans(is.na(bugs.dat$MSS))<.9)
  JSS       <- length(SS_index)
  MSS.case = cbind(MSS[,SS_index,drop=FALSE],MSS.only)
  ## cases:
  nonzero.pat=apply(MBS.case , 1 , paste , collapse = "" )
  casepat = sort(table(nonzero.pat),decreasing=TRUE)
  
  casepat.high = rev(sort(casepat))[1:npat]/Nd
  casepat.high.name = names(casepat.high)
  
  nonzero.pat=apply(MBS.case , 1 , paste , collapse = "" )
  ## look for individuals who have the interesting patterns:
  ind.interest.list = lapply(casepat.high.name,
                             function(x){res = which(nonzero.pat==x)})
  names(ind.interest.list) = casepat.high.name
  
  
  ## reading nplcm outputs:
  res_nplcm <- read.coda(paste(DIR_NPLCM,"coda1.txt",sep="/"),
                         paste(DIR_NPLCM,"codaIndex.txt",sep="/"),
                         quiet=TRUE)
  
  # prepare results for individual diagnosis:
  d1_nplcm  = nrow(res_nplcm)
  d2_nplcm  = Nd
  Icat_nplcm = array(NA,c(d1_nplcm,d2_nplcm))
  for (i in 1:d2_nplcm){
    SubVarName = paste(VarName,"[",i,"]",sep="")
    Icat_nplcm[,i] = res_nplcm[,SubVarName]
  }
  
  
  for (i in 1:length(casepat.high.name)){
    ind.tmp = ind.interest.list[[i]]
    
    #check all the patterns for cases in ind.tmp:
    nonzero.pat_MSS <- apply(MSS.case[ind.tmp,,drop=FALSE],1,paste,collapse="")
    MSS_pat <- sort(table(nonzero.pat_MSS),decreasing=TRUE)
    
    MSS_pat_high = rev(sort(MSS_pat))/(length(ind.tmp))
    MSS_pat_high.name = names(MSS_pat_high)
    
    ind.interest.list.MSS = lapply(MSS_pat_high.name,
                                   function(x){res = which(nonzero.pat_MSS==x)})
    # if you want to have patients with the same BrS pattern in the same plot,,
    # then delete the comment sign below:
    # par(mfrow=c(1,length(ind.interest.list.MSS)))
    for (j in 1:length(ind.interest.list.MSS)){
      #prediction for each MSS pattern within a MBS pattern:
      #nplcm:
      pred_nplcm_tmp = matrix(NA,nrow=length(ind.interest.list.MSS[[j]]),ncol=Jcause)
      count.flag = 0
      for (ii in ind.interest.list.MSS[[j]]){
        count.flag = count.flag+1
        pred_nplcm_tmp[count.flag,] = sapply(1:Jcause,
                                             function(j){sum(Icat_nplcm[,ind.tmp[ii]]==j)})/d1_nplcm
      }
      pred_nplcm = colMeans(pred_nplcm_tmp) # d1_nplcm for MCMC samples.
      
      
      #start plotting:
      par(mar=c(5.1, 4.1, 4.1, 2.1)+c(3,0,3,0),xpd=TRUE)
      barplot(rbind(pred_nplcm),beside=TRUE,ylim=c(0,1.1),
              col=c("dodgerblue2"),yaxt="n")
      axis(side = 2,at=seq(0,1,by=0.1),labels=seq(0,1,by=0.1),cex.axis=2,las=2)
      #mtext(casepat.high.name[i],side=3,line=-2,cex=3)
      #BrS:
      casepat.high.name.str <- unlist(strsplit(casepat.high.name[i],split=""))
      
      if (!is.null(MSS.only)){
        # if has SSonly measurements:
        for (j2 in 1:(JBrS+JSSonly)){
          points(2*(j2)-.5,1,pch=casepat.high.name.str[j2],cex=1.5)
        }
        #mtext(paste0(round(length(ind.tmp)/Nd,3)*100,"%"),side=3,line=-4,cex=2)
        text(2*((JBrS+JSSonly)/2)-.5,0.9,paste0(round(length(ind.tmp)/Nd,3)*100,"%"),
             cex=1.5)
        #SS:
        MSS_pat_high.name.str <- unlist(strsplit(MSS_pat_high.name[j],split=""))
        for (j2 in 1:JBrS){
          points(2*(j2)-.5,0.8,pch=MSS_pat_high.name.str[j2],col="red",cex=1.5)
        }
        for (j2 in 1:JSSonly){
          points(2*(j2+JBrS)-.5,0.8,pch=MSS_pat_high.name.str[j2+JSS],col="red",cex=1.5)
        }
        text(2*((JBrS+JSSonly)/2)-.5,.7,paste0(round(length(ind.interest.list.MSS[[j]])/length(ind.tmp),3)*100,"%"),col="red",
             cex=1.5)
        axis(3,at=2*(1:(JBrS+JSSonly))-.5,labels=pathogen_list,
             las=2,line=-0.5)
      }else{
        # if no SSonly measurements:
        for (j2 in 1:(JBrS)){
          points(2*(j2)-.5,1,pch=casepat.high.name.str[j2],cex=1.5)
        }
        #mtext(paste0(round(length(ind.tmp)/Nd,3)*100,"%"),side=3,line=-4,cex=2)
        text(2*((JBrS)/2)-.5,0.9,paste0(round(length(ind.tmp)/Nd,3)*100,"%"),
             cex=1.5)
        #SS:
        MSS_pat_high.name.str <- unlist(strsplit(MSS_pat_high.name[j],split=""))
        for (j2 in 1:JBrS){
          points(2*(j2)-.5,0.8,pch=MSS_pat_high.name.str[j2],col="red",cex=1.5)
        }
        text(2*((JBrS)/2)-.5,.7,paste0(round(length(ind.interest.list.MSS[[j]])/length(ind.tmp),3)*100,"%"),col="red",
             cex=1.5)
        axis(3,at=2*(1:(JBrS))-.5,labels=pathogen_list,
             las=2,line=-0.5)
      }
      
      # add labels to the top and bottom symbols:
      text(-2,1.1,"Pathogen:",adj=0)
      axis(1,at=2*(1:Jcause)-.5,labels=model_options$cause_list,
           las=2,line=0.5)
      
      text(-2,-0.1,"Cause:",adj=0)
    }
  }
}#END of function








