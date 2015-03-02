#' combine results from two folder, specially made for data visualization:
#'
#'This function creates visualizations that are convenient for comparing the
#'individual prediction from two models: pLCM and npLCM.
#'
#' @param DIR_NPLCM The folder path to results from npLCM fit
#' @param DIR_PLCM The folder path to results from pLCM fit
#' @param npat The number of bronze-standard patterns to display individual diagnosis
#' plots
#'
#' @return Individual plots comparing two nested models: pLCM and npLCM
#'
#' @importFrom coda read.coda
#' @export
#'

combined_visualization <- function(DIR_NPLCM,DIR_PLCM,npat=16){

      ## common to both directories:
      # the two folder should use the same data:
      # !!! remember that the data.txt file in winbugs working folder is transposed:
#       bugs2jags(paste(DIR_NPLCM,"data.txt",sep="/"),paste(DIR_NPLCM,"data_R_format.R",sep="/")) #Creates file "epildata.R"
#
#       source(paste(DIR_NPLCM,"data_R_format.R",sep="/"))

      bugs.dat <- dget(paste(DIR_NPLCM,"data.txt",sep="/"))
      for (bugs.variable.name in names(bugs.dat)) {
        if (!is.null(dim(bugs.dat[[bugs.variable.name]]))) {
          dim(bugs.dat[[bugs.variable.name]]) <- rev(dim(bugs.dat[[bugs.variable.name]]))
          bugs.dat[[bugs.variable.name]] <- aperm(bugs.dat[[bugs.variable.name]])
        }
        assign(bugs.variable.name, bugs.dat[[bugs.variable.name]])
      }


      model_options  <- dget(paste(DIR_NPLCM,"model_options.txt",sep="/"))
      pathogen_list <- model_options$pathogen_list
      VarName ="Icat"

      MBS   = bugs.dat$MBS
      MSS   = bugs.dat$MSS
      Jfull = ncol(MBS)

      Nd = bugs.dat$Nd
      Nu = bugs.dat$Nu
      Y = c(rep(1,Nd),rep(0,Nu))

      MBS.case = bugs.dat$MBS[Y==1,]
      SS_index  <- which(colMeans(is.na(bugs.dat$MSS))<.9)
      JSS       <- length(SS_index)
      MSS.case = MSS[,SS_index]
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
                           paste(DIR_NPLCM,"codaIndex.txt",sep="/"),quiet=TRUE)
    ## reading plcm outputs:
    res_plcm <- read.coda(paste(DIR_PLCM,"coda1.txt",sep="/"),
                         paste(DIR_PLCM,"codaIndex.txt",sep="/"),quiet=TRUE)

    # prepare results for individual diagnosis:
    d1_nplcm  = nrow(res_nplcm)
    d2_nplcm  = Nd
    Icat_nplcm = array(NA,c(d1_nplcm,d2_nplcm))
    for (i in 1:d2_nplcm){
      SubVarName = paste(VarName,"[",i,"]",sep="")
      Icat_nplcm[,i] = res_nplcm[,SubVarName]
    }


    d1_plcm  = nrow(res_plcm)
    d2_plcm  = Nd
    Icat_plcm = array(NA,c(d1_plcm,d2_plcm))
    for (i in 1:d2_plcm){
      SubVarName = paste(VarName,"[",i,"]",sep="")
      Icat_plcm[,i] = res_plcm[,SubVarName]
    }



    for (i in 1:length(casepat.high.name)){
      ind.tmp = ind.interest.list[[i]]

      #check all the patterns for cases in ind.tmp:
      nonzero.pat_MSS <- apply(MSS.case[ind.tmp,],1,paste,collapse="")
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
        pred_nplcm_tmp = matrix(NA,nrow=length(ind.interest.list.MSS[[j]]),ncol=Jfull)
        count.flag = 0
        for (ii in ind.interest.list.MSS[[j]]){
          count.flag = count.flag+1
          pred_nplcm_tmp[count.flag,] = sapply(1:Jfull,function(j){sum(Icat_nplcm[,ind.tmp[ii]]==j)})/d1_nplcm
        }
        pred_nplcm = colMeans(pred_nplcm_tmp) # d1_nplcm for MCMC samples.



        #prediction for each MSS pattern within a MBS pattern:
        #plcm:
        pred_plcm_tmp = matrix(NA,nrow=length(ind.interest.list.MSS[[j]]),ncol=Jfull)
        count.flag = 0
        for (ii in ind.interest.list.MSS[[j]]){
          count.flag = count.flag+1
          pred_plcm_tmp[count.flag,] = sapply(1:Jfull,function(j){sum(Icat_plcm[,ind.tmp[ii]]==j)})/d1_plcm
        }
        pred_plcm = colMeans(pred_plcm_tmp) # d1_plcm for MCMC samples.

        #start plotting:
        barplot(rbind(pred_plcm,pred_nplcm),beside=TRUE,ylim=c(0,1.1),
                col=c("white","dodgerblue2"),yaxt="n")
        axis(side = 2,at=seq(0,1,by=0.1),labels=seq(0,1,by=0.1),cex.axis=2,las=2)
        #mtext(casepat.high.name[i],side=3,line=-2,cex=3)
        #BrS:
        casepat.high.name.str <- unlist(strsplit(casepat.high.name[i],split=""))
        for (j2 in 1:Jfull){
          points(3*(j2)-1,1,pch=casepat.high.name.str[j2],cex=1.5)
        }
        #mtext(paste0(round(length(ind.tmp)/Nd,3)*100,"%"),side=3,line=-4,cex=2)
        text(3*(Jfull/2)-1,0.9,paste0(round(length(ind.tmp)/Nd,3)*100,"%"),
             cex=1.5)
        #SS:
        MSS_pat_high.name.str <- unlist(strsplit(MSS_pat_high.name[j],split=""))
        for (j2 in 1:Jfull){
          points(3*(j2)-1,0.8,pch=MSS_pat_high.name.str[j2],col="red",cex=1.5)
        }
        #mtext(paste0(MSS_pat_high.name[j],paste(rep("-",Jfull-JSS),collapse="")),
        #      side=3,line=-6,cex=3,col="red")
        text(3*(Jfull/2)-.1,.7,paste0(round(length(ind.interest.list.MSS[[j]])/length(ind.tmp),3)*100,"%"),col="red",
             cex=1.5)
        #mtext(paste0(round(length(ind.interest.list.MSS[[j]])/length(ind.tmp),3)*100,"%"),
        #      side=3,line=-8,cex=2,col="red")
        #                  for (j2 in 1:Jfull){
        #                    text(2*(j2)-.5,0.5,model_options$pathogen_list[j2],srt=90,cex=.8)
        #                  }
        axis(3,at=3*(1:Jfull)-1,labels=pathogen_list,
             las=2,line=-1)
      }
    }
}
