#' Plot three-panel figures for nested partially-latent model results
#'
#' Visualize the model outputs for communicating how data inform the
#' DN: 1. current implementation: nplcm, BrS and SS.
#' "Jfull" here is not the same as in other functions: it refers to the number of
#' pathogens even if there are pathogens with only silver-standard data;
#' in other functions, "Jfull" refers to the number of pathogens that have BrS data.
#' 2. Missing data for BrS or SS are dropped when calculating fractionsa
#'
#' @param DIR_NPLCM directory of results
#'
#' @param ss_upperlimit The upper limit of horizontal bar for the silver-standard
#' subpanel (the middle panel). The default value is .25.
#'
#' @param eti_upperlimit The upper limit of horizontal bar for the etiology
#' posterior subpanel (the rightmost panel). The default value is .4
#' @importFrom coda read.coda
#' @importFrom binom binom.confint
#' @return None
#'
#' @export
nplcm_plot_three_panel <- function(DIR_NPLCM,
                                   ss_upperlimit=1,
                                   eti_upperlimit=.4){#BEGIN function

            #read in data from the folder directory:
            bugs.dat <- dget(paste(DIR_NPLCM,"data.txt",sep="/"))
            for (bugs.variable.name in names(bugs.dat)) {
                  if (!is.null(dim(bugs.dat[[bugs.variable.name]]))) {
                    dim(bugs.dat[[bugs.variable.name]]) <- rev(dim(bugs.dat[[bugs.variable.name]]))
                    bugs.dat[[bugs.variable.name]] <- aperm(bugs.dat[[bugs.variable.name]])
                  }
                  assign(bugs.variable.name, bugs.dat[[bugs.variable.name]])
            }

            model_options  <- dget(paste(DIR_NPLCM,"model_options.txt",sep="/"))
            #compatibility checking:
            if (length(model_options$M_use)!=length(model_options$TPR_prior)){
              stop("The number of measurement source(s) is different from
                         the number of TPR prior option!
                         Make them equal, and match with order!")
            }
            clean_options <- dget(paste(DIR_NPLCM,"data_clean_options.txt",sep="/"))

            #some data preparation:
            Nd <- bugs.dat$Nd
            Nu <- bugs.dat$Nu

            Y = c(rep(1,Nd),rep(0,Nu))

            model_data_source <- rep(NA,3)
            names(model_data_source) <- c("MBS","MSS","MGS")
            model_data_source[1] <- c("no","yes")["BrS"%in%model_options$M_use+1]
            model_data_source[2] <- c("no","yes")["SS"%in%model_options$M_use+1]
            model_data_source[3] <- c("no","yes")["GS"%in%model_options$M_use+1]

            pathogen_list     <- model_options$pathogen_list

            Jfull_BrS         <- length(pathogen_list)

            #reading nplcm outputs:
            res_nplcm <- read.coda(paste(DIR_NPLCM,"coda1.txt",sep="/"),
                                   paste(DIR_NPLCM,"codaIndex.txt",sep="/"),
                                   quiet=TRUE)

            Jfull <- length(grep("pEti",colnames(res_nplcm)))

            # extract and process some data and posterior samples:
            SubVarName <- rep(NA,Jfull)
            for (j in 1:Jfull){
                  SubVarName[j] = paste("pEti","[",j,"]",sep="")
            }

            #get etiology fraction MCMC samples:
            pEti_mat   <- res_nplcm[,SubVarName]
            pEti_mean  <- colMeans(pEti_mat)
            pEti_mean0 <- pEti_mean

            # order the pathogens by posterior mean:
            ord <- order(pEti_mean)
            if (is.null(model_options$SSonly) || model_options$SSonly==FALSE){
                  pathogens_ord <- model_options$pathogen_list[ord]
            } else{
                  pathogens_ord <- c(model_options$pathogen_list,
                                     model_options$pathogen_SSonly_list)[ord]
            }
            pEti_mean_ord <- pEti_mean[ord]
            pEti_mat_ord  <- pEti_mat[,ord]

            # quantiles for etiology: outer is 97.5% CI, inner is 50% CI

            pEti_q1   <- apply(pEti_mat,2,quantile,probs=0.025)[ord]
            pEti_q2   <- apply(pEti_mat,2,quantile,probs=0.975)[ord]
            pEti_innerq1   <- apply(pEti_mat,2,quantile,probs=0.25)[ord]
            pEti_innerq2   <- apply(pEti_mat,2,quantile,probs=0.75)[ord]

            ncase <- sum(Y==1)
            nctrl <- sum(Y==0)

            MBS <- bugs.dat$MBS
            colnames(MBS) <-model_options$pathogen_list
            if (is.null(clean_options$allow_missing)||
                  clean_options$allow_missing==FALSE){
                  tmp.case <- binom.confint(colSums(MBS[1:ncase,]),
                                            ncase, conf.level = 0.95, methods = "ac")[ord,]
                  tmp.ctrl <- binom.confint(colSums(MBS[-(1:ncase),]),
                                            nctrl, conf.level = 0.95, methods = "ac")[ord,]
            }else{
                  ind_case_BrS_dat_not_na <- which(rowSums(is.na(MBS[1:ncase,]))==0)
                  ind_ctrl_BrS_dat_not_na <- which(rowSums(is.na(MBS[-(1:ncase),]))==0)

                  tmp.case <- binom.confint(colSums(MBS[(1:ncase)[ind_case_BrS_dat_not_na],]),
                                            ncase, conf.level = 0.95, methods = "ac")[ord,]
                  tmp.ctrl <- binom.confint(colSums(MBS[(ncase+(1:nctrl))[ind_ctrl_BrS_dat_not_na],]),
                                            nctrl, conf.level = 0.95, methods = "ac")[ord,]
            }

              # case and control positive rate, lower and upper limit
              Bcomp   <- rbind(round(tmp.case$mean,5),round(tmp.ctrl$mean,5))
              Bcomp_q1 <- rbind(tmp.case[,c("lower")],tmp.ctrl[,c("lower")])
              Bcomp_q2 <- rbind(tmp.case[,c("upper")],tmp.ctrl[,c("upper")])

            #posterior distribution of TPR and FPR:

            if (model_options$k_subclass>1){
                  # nplcm (need to marginalize over subclasses to calculate marginal TPR):
                  theta_mat  <- (res_nplcm[,grep("ThetaBS.marg",colnames(res_nplcm))])[,ord]
                  theta_mean <- colMeans(theta_mat)

                  #posterior distribution of FPR
                  psi_mat  <- (res_nplcm[,grep("PsiBS.marg",colnames(res_nplcm))])[,ord]
                  psi_mean <- colMeans(psi_mat)

                  psi_mat.case = array(NA,c(nrow(res_nplcm),Jfull,Jfull),
                                        dimnames = list(NULL,1:Jfull,1:Jfull))

                  psi_mat.case.tmp = (res_nplcm[,grep("PsiBS.case",colnames(res_nplcm))])
                  for (j in 1:Jfull){
                    for (s in 1:Jfull){
                      tmp.nm = paste0("PsiBS.case","[",j,",",s,"]")
                      ind.tmp = which(colnames(psi_mat.case.tmp)==tmp.nm)
                      psi_mat.case[,j,s] = psi_mat.case.tmp[,ind.tmp]
                    }
                  }

                  psi_mat.case.ord = psi_mat.case[,ord,ord]

                  # model fitted postive rate for each pathogen
                  fittedmean_case <- sapply(1:Jfull,function(s) mean(pEti_mat_ord[,s]*theta_mat[,s]+
                                                   rowSums(pEti_mat_ord[,-s]*psi_mat.case.ord[,s,-s])))
                  fittedmean_control <- psi_mean
            } else {
                  # plcm (just obtain TPR and FPR estimates):
                  if (is.null(model_options$SSonly)|| model_options$SSonly==FALSE){
                     theta_mat <- (res_nplcm[,grep("thetaBS",colnames(res_nplcm))])[,ord]
                  }else {
                     theta_mat <- cbind(res_nplcm[,grep("thetaBS",colnames(res_nplcm))],
                                       matrix(NA,nrow=nrow(res_nplcm),ncol=Jfull-Jfull_BrS))[,ord]
                  }
                  theta_mean <- colMeans(theta_mat)

                  #posterior distribution of FPR
                  if (is.null(model_options$SSonly)||model_options$SSonly==FALSE){
                      psi_mat  <- (res_nplcm[,grep("psiBS",colnames(res_nplcm))])[,ord]
                  }else{
                      psi_mat  <- cbind(res_nplcm[,grep("psiBS",colnames(res_nplcm))],
                                        matrix(NA,nrow=nrow(res_nplcm),ncol=Jfull-Jfull_BrS))[,ord]
                  }
                  psi_mean <- colMeans(psi_mat)

                  ## model fitted postive rate for each pathogen
                  fittedmean_case <- sapply(1:Jfull,
                                            function(s) mean(pEti_mat_ord[,s]*theta_mat[,s]+
                                                               (1-pEti_mat_ord[,s])*psi_mat[,s]))
                  fittedmean_control <- psi_mean

            }


            ## start plotting:
            # if only BrS data is used:
            if (model_data_source[1]=="yes" &
                    model_data_source[2]=="no" &
                       model_data_source[3]=="no"){
                        layout(matrix(c(1,2),1,2,byrow = TRUE),
                                              widths=c(4,4),heights=c(8))
            }

            # if both BrS and SS data are used:
            if (model_data_source[1]=="yes" &
                    model_data_source[2]=="yes" &
                       model_data_source[3]=="no"){

              layout(matrix(c(1,2,3),1,3,byrow = TRUE),
                     widths=c(3,2,3),heights=c(8))
            }
            cexval=1
            srtval=0
            ## BEGIN subplot 1: BrS Info -----------------------------------------------
            top2 = 1.3
            op<-par(mar=c(5.1,4.1,4.1,0))
            plotat = seq(0.5,Jfull+0.5,by=1/4)[-(c(1,(1:Jfull)*4+1))]
            #plot case control positive rates, and fitted case rates
            plot(c(rbind(fittedmean_case,Bcomp)),plotat,yaxt="n",xlim=c(0,top2),
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

            tmp0 = COR(MBS,ncase,pathogens_ord[ord<=Jfull_BrS])
            tmp = tmp0$ORinterval

            #plot conditional odds ratio on the right:
            incre = 0
            for (s in (1:Jfull)){
                  if (s %in% which(ord<=Jfull_BrS)){
                      incre = incre + 1
                      L=round(tmp[incre,1],1)
                      C=round(tmp[incre,3],1)
                      R=round(tmp[incre,2],1)
                      text(top2-0.12,s+1/(2*Jfull),C,cex=1.5)
                      text(top2-0.12,s-.2,paste(c(L,"   ",R),collapse=" "),cex=1.2)
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
                text(tmp.hpos,plotat[s],paste0(round(100*c(Bcomp),1)[counter],"%"),srt=srtval,cex=cexval)
              }
            }

            for (s in 1:(Jfull)){
              segments(y0=plotat[3*s-1],x0=c(rbind(fittedmean_case,Bcomp))[3*s-1],
                       y1=plotat[3*s],x1=c(rbind(fittedmean_case,Bcomp))[3*s],col="dodgerblue2",lwd=2)
            }
            # put prior shapes on bronze-standard sensitivity

            TPR_prior_list <-  bugs.dat#TPR_prior_set(model_options,Mobs,Y)

            # if both BrS and SS data are used:
            if (model_data_source[1]=="yes" &
                  model_data_source[2]=="no" &
                    model_data_source[3]=="no"){
                       alphaB <- TPR_prior_list$alphaB
                       betaB <- TPR_prior_list$betaB
            }
            # if both BrS and SS data are used:
            if (model_data_source[1]=="yes" &
                  model_data_source[2]=="yes" &
                   model_data_source[3]=="no"){
                     alphaB <- TPR_prior_list$alphaB
                     betaB <- TPR_prior_list$betaB
                     #alphaS <- TPR_prior_list$alphaS
                     #betaS <- TPR_prior_list$betaS
                     if (!is.null(model_options$SSonly) && model_options$SSonly==TRUE){
                          alphaS.only <- TPR_prior_list$alphaS.only
                          betaS.only <- TPR_prior_list$betaS.only
                          alphaS <- c(TPR_prior_list$alphaS,alphaS.only)
                          betaS <- c(TPR_prior_list$betaS,betaS.only)
                     }
                    used_cat <- TPR_prior_list$used_cat
            }


            for (s in 1:Jfull){
              if (s %in% which(ord<=Jfull_BrS)){
                  tmp = rbeta(10000,alphaB[ord[s]],betaB[ord[s]])
                  boxplot(tmp,at = s-0.45, boxwex=1/10 , col="gray",
                          add=TRUE,horizontal=TRUE,outline=FALSE,xaxt="n")
                  tmp.post = as.matrix(theta_mat)[,s]
                  boxplot(tmp.post,at = s-0.35,boxwex=1/10,add=TRUE,
                          horizontal=TRUE,outline=FALSE,xaxt="n")
              }
            }
            axis(2,at=(1:Jfull)-.45,labels=rep("",Jfull),las=2,cex.axis=.5)
            axis(2,at=(1:Jfull)-.35,labels=rep("",Jfull),las=2,cex.axis=.5)

            mtext(expression(underline("BrS")),line=1,cex=1.8)
            par(op)
            ## END subplot 1: BrS Info--------------------------------------------------

            # if both BrS and SS data are used:
            if (model_data_source[1]=="yes" &
                  model_data_source[2]=="yes" &
                    model_data_source[3]=="no"){

                    ## BEGIN subplot 2: SS Info -----------------------------------------------

                    ## check DN about how to warn user about missing some BCX:

                    SSdat = bugs.dat$MSS

                    SS_index  <- which(colMeans(is.na(SSdat))<.9)
                    JSS       <- length(SS_index)

                    ind.SS = rep(NA,JSS) # tells where the the SS row should go.
                    for (j in 1:JSS){
                      ind.SS[j] = which(ord==j)
                    }
                    if (!is.null(model_options$SSonly) && model_options$SSonly==TRUE){
                                SSonlydat = bugs.dat$MSS.only
                                JSS_only   = Jfull-Jfull_BrS
                                ind.SSonly = rep(NA,JSS_only)
                                for (j in 1:(JSS_only)){
                                  ind.SSonly[j] = which(ord == j+Jfull_BrS)
                                }
                                SSdat = cbind(bugs.dat$MSS,SSonlydat)
                                SS_index = which(colMeans(is.na(SSdat))<.9)
                                JSS      = length(SS_index)

                                ind.SS = rep(NA,JSS)
                                for (j in 1:JSS){
                                  ind.SS[j] = which(ord==ifelse(j<=JSS-JSS_only,j,
                                                                j-JSS+JSS_only+Jfull_BrS))
                                }
                    }

                    if (is.null(clean_options$allow_missing)||
                            clean_options$allow_missing==FALSE){
                          tmpSS.case = binom.confint(colSums(SSdat[,1:JSS]), nrow(SSdat),
                                                    conf.level = 0.95, methods = "ac")
                    }else{
                          ind_SSdat_not_na <- which(rowSums(is.na(SSdat[,1:JSS]))==0)
                          tmpSS.case = binom.confint(colSums(SSdat[ind_SSdat_not_na,1:JSS]),
                                                     length(ind_SSdat_not_na),
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

                    par(mar=c(5.1,0,4.1,0))

                    plotat = seq(0.5,Jfull+0.5,by=1/4)[-(c(1,(1:Jfull)*4+1))]
                    #plotat.short = plotat[1:length(c(rbind(thetameanG,Gcomp)))]
                    plotat.calc = function(j) {c(3*j-2,3*j-1,3*j)}
                    plotat.short = rep(NA,JSS*3)
                    for (j in 1:JSS){
                      plotat.short[c(3*j-2,3*j-1,3*j)] = plotat[plotat.calc(ind.SS[j])]
                    }


                    top3 = ss_upperlimit
                    plot(c(rbind(theta_meanSS,SScomp)),plotat.short,yaxt="n",xlim=c(0,top3),
                         ylim=c(0.5,Jfull+.5),#xaxt="n",
                         ylab="",xlab="probability",
                         pch = c(rbind(rep(20,Jfull),rep(20,Jfull),rep(20,Jfull))),
                         col=c(rbind(rbind(rep(1,Jfull),rep("blue",Jfull),rep(1,Jfull)))),
                         cex = c(rbind(rep(1,Jfull),rep(2,Jfull),rep(2,Jfull))))
                    #axis(1,at = seq(0,top3,len=10),labels= seq(0,top3,len=10),las=1)
                    points(c(rbind(fittedmean_SS_pos,matrix("",nrow=2,ncol=JSS))),plotat.short,
                             yaxt="n",xlim=c(0,top3),
                           ylim=c(0.5,Jfull+.5),xaxt="n",
                           ylab="",#xlab="Gold Positive Rate",
                           pch = c(rbind(rep(2,Jfull),rep(NA,Jfull),rep(NA,Jfull))),
                           col=c(rbind(rbind(rep(1,Jfull),rep(1,Jfull),rep(1,Jfull)))),
                           cex = c(rbind(rep(1,Jfull),rep(2,Jfull),rep(2,Jfull))))

                    #axis(2,at = plotat.short,labels=rep(c("TPR-G","case","ctrl"),JSS),las=2)
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
                    par(op)
                    ## END subplot 2: SS Info--------------------------------------------------
            }

            ## BEGIN subplot 3: etiology Info -------------------------------------------
            top=eti_upperlimit#max(pq2BG)
            dotcolor = "black"
            #op <- par(mar=c(5.1,6,4.1,1.1))
            op <- par(mar=c(5.1,0,4.1,9))
            plot(c(pEti_mean_ord),1:(Jfull),
                 yaxt="n",#xaxt="n",
                 xlim=c(0,top),ylim=c(0.5,Jfull+0.5),col=c("black"),
                 ylab="",xlab="probability",
                 pch=c(20),cex=2)
            axis(4,at=1:Jfull,labels=pathogens_ord,las=2,cex.axis=1.5)
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
              segments(y0=s,x0=c(pEti_innerq1)[s],y1=s,x1=c(pEti_innerq2)[s],col=dotcolor,lwd=2)
              #text(pmeanBG[s],s+0.25,paste0(round(100*c(pmeanBG),2)[s],"%"),srt=srtval,cex=cexval)
              text(.35,s+0.25,paste0("=",paste0(round(100*c(pEti_mean_ord),1)[s],"%")),srt=srtval,cex=2)
              text(.27,s+0.25,bquote(hat(pi)[.(ord[s])]),srt=srtval,cex=2)
              #text(top-0.05,s,pathogens[s],cex=.8,srt=srtval)
              tmp.density = dbeta(pgrid,alpha[ord[s]],sum(alpha[-ord[s]]))
              points(pgrid,tmp.density/(3*max(tmp.density))+s-0.45,type="l",col="gray",lwd=4,lty=2)
              ##posterior density
              tmp.post.density = density(pEti_mat_ord[,s],from=0,to=1)
              tmp.x = tmp.post.density$x
              tmp.y = tmp.post.density$y
              points(tmp.x,tmp.y/(3*max(tmp.y))+s-0.45,col="black",lwd=4,type="l")

            }

            par(op)
            ## END subplot 3: etiology Info----------------------------------------------

}# END function
