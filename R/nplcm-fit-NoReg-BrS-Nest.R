if(getRversion() >= "2.15.1") utils::globalVariables(c("equals","Icat",
                                                       "Z","ThetaBS","Icat.new",
                                                       "Z.new","PsiBS","inprod2"))


#' Fit nested partially-latent class model (low-level)
#'
#' Features:
#' \itemize{
#' \item no regression;
#' \item bronze- (BrS) measurements;
#' \item conditional dependence;
#' \item all pathogens have BrS measurements.
#' }
#'
#' @inheritParams nplcm
#' @return WinBUGS fit results.
#'
#' @export

nplcm_fit_NoReg_BrS_Nest <-
  function(data_nplcm,model_options,mcmc_options){
    # Record the settings of current analysis:
    cat("==Results stored in: ==","\n",mcmc_options$result.folder)
    #model_options:
    dput(model_options,file.path(mcmc_options$result.folder,"model_options.txt"))
    #mcmc_options:
    dput(mcmc_options,file.path(mcmc_options$result.folder,"mcmc_options.txt"))
    
    Mobs <- data_nplcm$Mobs
    Y    <- data_nplcm$Y
    
    # define generic function to call WinBUGS:
    call.bugs <- function(data, inits, parameters,m.file,
                          bugsmodel.dir = mcmc_options$bugsmodel.dir,
                          winbugs.dir   = mcmc_options$winbugs.dir,
                          nitermcmc     = mcmc_options$n.itermcmc,
                          nburnin       = mcmc_options$n.burnin,
                          nthin         = mcmc_options$n.thin,
                          nchains       = mcmc_options$n.chains,
                          dic = FALSE,
                          is.debug = mcmc_options$debugstatus,
                          workd= mcmc_options$result.folder,...) {

      m.file <- file.path(bugsmodel.dir, m.file);
      f.tmp <- function() {
        ##winbugs
        gs <- R2WinBUGS::bugs(data, inits, parameters,
                   model.file = m.file,
                   working.directory=workd,
                   n.chains = nchains,
                   n.iter   = nitermcmc,
                   n.burnin = nburnin,
                   n.thin   = nthin,
                   bugs.directory=winbugs.dir,
                   DIC=dic,
                   debug=is.debug,...);

        gs;
      }

      bugs.try  <- try(rst.bugs <- f.tmp(), silent=FALSE);
      if (class(bugs.try) == "try-error") {
        rst.bugs <- NULL;
      }
      rst.bugs;
    }

    #-------------------------------------------------------------------#
    # prepare data:
    #parsing <- assign_model(data_nplcm,model_options)
    Nd <- sum(Y==1)
    Nu <- sum(Y==0)

    cat("==True positive rate (TPR) prior(s) for ==\n",
        model_options$M_use,"\n",
        " is(are respectively): \n",
        model_options$TPR_prior,"\n")

    cause_list              <- model_options$cause_list
    pathogen_BrS_list       <- model_options$pathogen_BrS_list
    pathogen_SSonly_list    <- model_options$pathogen_SSonly_list

    # get the count of pathogens:
    # number of all BrS available pathogens:
    JBrS         <- length(pathogen_BrS_list)
    # number of all SS only pathogens:
    JSSonly      <- length(pathogen_SSonly_list)
    # number of all causes possible: singletons, combos, NoA, i.e.
    # the number of rows in the template:
    Jcause      <- length(cause_list)

    template <- rbind(as.matrix(rbind(symb2I(c(cause_list),
                                             c(pathogen_BrS_list,pathogen_SSonly_list)))),
                      rep(0,JBrS+JSSonly)) # last row for controls.

    # fit model :
    # plcm - BrS + SS and SSonly:
    MBS.case <- Mobs$MBS[Y==1,]
    MBS.ctrl <- Mobs$MBS[Y==0,]
    MBS      <- as.matrix(rbind(MBS.case,MBS.ctrl))

    #MSS.case <- Mobs$MSS[Y==1,1:JBrS]
    #MSS.case <- as.matrix(MSS.case)

    #SS_index <- which(colMeans(is.na(MSS.case))<0.9)#.9 is arbitrary; any number <1 will work.
    #JSS      <- length(SS_index)
    #MSS      <- MSS.case[,SS_index]

    # set priors:
    alpha          <- set_prior_eti(model_options)

    TPR_prior_list <- set_prior_tpr(model_options,data_nplcm)
    alphaB      <- TPR_prior_list$alphaB
    betaB       <- TPR_prior_list$betaB
    #alphaS      <- TPR_prior_list$alphaS
    #betaS       <- TPR_prior_list$betaS
#     if (parsing$measurement$SSonly){
#       MSS.only.case <- Mobs$MSS[Y==1,(1:JSSonly)+JBrS]
#       MSS.only <- as.matrix(MSS.only.case)
#       alphaS.only <- TPR_prior_list$alphaS.only
#       betaS.only  <- TPR_prior_list$betaS.only
#     }


    K        <- model_options$k_subclass

    mybugs <- function(...){
      inits      <- function(){list(pEti = rep(1/Jcause,Jcause),
                                    r0 = c(rep(.5,K-1),NA),
                                    r1 = cbind(matrix(rep(.5,Jcause*(K-1)),
                                                      nrow=Jcause,ncol=K-1),
                                               rep(NA,Jcause)),
                                    alphadp0 = 1)};
      data       <- c("Nd","Nu","JBrS","Jcause",
                      "alpha","template","K",
                      #"JSS","MSS",
                      #"MSS.only","JSSonly","alphaS.only","betaS.only",
                      "MBS","alphaB","betaB");

      if (mcmc_options$individual.pred==FALSE & mcmc_options$ppd==TRUE){
        parameters <- c("pEti","Lambda","Eta","alphadp0","MBS.new",
                        "ThetaBS.marg","PsiBS.marg","PsiBS.case",
                        "ThetaBS","PsiBS")

      } else if(mcmc_options$individual.pred==TRUE & mcmc_options$ppd==TRUE){
        parameters <- c("pEti","Lambda","Eta","alphadp0","Icat","MBS.new",
                        "ThetaBS.marg","PsiBS.marg","PsiBS.case",
                        "ThetaBS","PsiBS")

      } else if (mcmc_options$individual.pred==TRUE & mcmc_options$ppd==FALSE){
        parameters <- c("pEti","Lambda","Eta","alphadp0","Icat",
                        "ThetaBS.marg","PsiBS.marg","PsiBS.case",
                        "ThetaBS","PsiBS")

      } else if (mcmc_options$individual.pred==FALSE & mcmc_options$ppd==FALSE){
        parameters <- c("pEti","Lambda","Eta","alphadp0",
                        "ThetaBS.marg","PsiBS.marg","PsiBS.case",
                        "ThetaBS","PsiBS")

      }
      rst.bugs   <- call.bugs(data, inits, parameters,...);
      rst.bugs
    }
    
    
    #-----------------BEGINNING OF MODELS----------------------------------#
    
    #
    # write the .bug files into mcmc_options$bugsmodel.dir; could later set it equal to result.folder.
    #
    
    # the one with posterior predictive distribution:
    model_NoReg_BrS_nplcm_ppd <- function(){
      #BEGIN MODEL
      ## 0) conditional dependence;
      ## 1) accommodates singletons, combos, and NoA;
      ## 2) check-bit to prevent case data informing FPR;
      
      ## case BrS measurements:
      for (i in 1:Nd){
        for (j in 1:JBrS){
          ind[i,j] <- equals(1,template[Icat[i],j])
          MBS[i,j]~dbern(mu_bs.bound[i,j])
          mu_bs.bound[i,j]<-max(0.000001,min(0.999999,mu_bs[i,j]))
          mu_bs[i,j]<-PR_BS[i,j,Z[i]]
          
          for (s in 1:K){
            PR_BS[i,j,s]<-PsiBS.cut[j,s]*(1-ind[i,j])+ThetaBS[j,s]*ind[i,j]
          }
          
          #predictive checking:
          ind.new[i,j] <- equals(1,template[Icat.new[i],j])
          MBS.new[i,j]~dbern(mu_bs.bound.new[i,j])
          mu_bs.bound.new[i,j]<-max(0.000001,min(0.999999,mu_bs.new[i,j]))
          mu_bs.new[i,j]<-PR_BS.new[i,j,Z.new[i]]
          
          for (s in 1:K){
            PR_BS.new[i,j,s]<-PsiBS.cut[j,s]*(1-ind.new[i,j])+ThetaBS[j,s]*ind.new[i,j]
          }
        }
      }
      
      ## cut the feedback from case model to FPR:
      for (j in 1:JBrS){
        for (s in 1:K){
          PsiBS.cut[j,s]<-cut(PsiBS[j,s])
        }
      }
      
      ## control BrS measurements
      for (i in (Nd+1):(Nd+Nu)){
        for (j in 1:JBrS){
          MBS[i,j]~dbern(mu_bs.bound[i,j])
          mu_bs.bound[i,j] <-max(0.000001,min(0.999999,mu_bs[i,j]))
          mu_bs[i,j]<-PsiBS[j,Z[i]]
          
          #predictive checking:
          MBS.new[i,j]~dbern(mu_bs.bound.new[i,j])
          mu_bs.bound.new[i,j] <-max(0.000001,min(0.999999,mu_bs.new[i,j]))
          mu_bs.new[i,j]<-PsiBS[j,Z.new[i]]
        }
      }
      
      ## Case SS measure on pathogens that also have BrS:
      #for (i in 1:Nd){
      #for (j in 1:JSS){
      #  MSS[i,j] ~ dbern(mu_ss[i,j])
      #  mu_ss[i,j]<-ind[i,j]*thetaSS[j]+(1-ind[i,j])*psiSS[j]
      #}
      
      #for (j in 1:JSSonly){
      #	  ind[i,j+JBrS]<- equals(1,template[Icat[i],j+JBrS])
      #	  MSS.only[i,j] ~ dbern(mu_ss.only[i,j])
      #	  mu_ss.only[i,j]<-ind[i,j+JBrS]*thetaSS.only[j]
      #}
      #}
      
      for (i in 1:Nd){
        Z[i] ~ dcat(Eta[Icat[i],1:K])
        Icat[i] ~ dcat(pEti[1:Jcause])
        
        #predictive checking:
        Z.new[i] ~ dcat(Eta[Icat.new[i],1:K])
        Icat.new[i] ~ dcat(pEti[1:Jcause])
      }
      
      ######################
      ##etiology prior
      ######################
      pEti[1:Jcause]~ddirch(alpha[])
      for (i in (Nd+1):(Nd+Nu)){
        Z[i]~dcat(Lambda[1:K])
        
        #predictive checking:
        Z.new[i]~dcat(Lambda[1:K])
      }
      
      ####################################
      ### stick-breaking specification 2
      ####################################
      Lambda0[1]<-r0[1]
      r0[K]<-1
      for(j in 2:K) {Lambda0[j]<-r0[j]*(1-r0[j-1])*Lambda0[j-1]/r0[j-1]}
      for(k in 1:K-1){
        r0[k]~dbeta(1,alphadp0)%_%I(0.000001,0.999999)
      }
      
      for (k in 1:K-1){Lambda[k]<-max(0.000001,min(0.999999,Lambda0[k]))}
      Lambda[K]<-1-sum(Lambda[1:(K-1)])
      
      for (s in 1:Jcause){
        Eta0[s,1]<-r1[s,1]
        r1[s,K]<-1
        for(j in 2:K) {Eta0[s,j]<-r1[s,j]*(1-r1[s,j-1])*Eta0[s,j-1]/r1[s,j-1]}
        for(k in 1:K-1){
          r1[s,k]~dbeta(1,alphadp0)%_%I(0.000001,0.999999)
        }
      }
      
      for (s in 1:Jcause){
        for (k in 1:K-1){Eta[s,k]<-max(0.000001,min(0.999999,Eta0[s,k]))}
        Eta[s,K]<-1-sum(Eta[s,1:(K-1)])
      }
      
      alphadp0~dgamma(.25,.25)%_%I(0.001,20)
      
      #########################
      ## priors on TPR and FPR:
      #########################
      
      for (j in 1:JBrS){
        for (s in 1:K){
          PsiBS[j,s]~dbeta(1,1)
          #ThetaBS[j,s]~dbeta(1,1)
          ThetaBS[j,s]~dbeta(alphaB[j],betaB[j])
        }
        ThetaBS.marg[j]<-inprod2(ThetaBS[j,1:K],Eta[j,1:K])
        PsiBS.marg[j]<-inprod2(PsiBS[j,1:K],Lambda[1:K])
        
        for (l in 1:Jcause){
          PsiBS.case[j,l]<-inprod2(PsiBS[j,1:K],Eta[l,1:K])
        }
      }
      
      # silver-standard measurement characteristics:
      #for (j in 1:JSS){
      #    thetaSS[j]~dbeta(alphaS[j],betaS[j])
      #    psiSS[j]<-0
      #}
      
      #for (j in 1:JSSonly){
      #	thetaSS.only[j]~dbeta(alphaS.only[j],betaS.only[j])
      #}
      
      
      
    }
    # the one without ppd:
    model_NoReg_BrS_nplcm <- function(){
      #BEGIN MODEL
      ## 0) conditional dependence;
      ## 1) accommodates singletons, combos, and NoA;
      ## 2) check-bit to prevent case data informing FPR;
      
      ## case BrS measurements:
      for (i in 1:Nd){
        for (j in 1:JBrS){
          ind[i,j] <- equals(1,template[Icat[i],j])
          MBS[i,j]~dbern(mu_bs.bound[i,j])
          mu_bs.bound[i,j]<-max(0.000001,min(0.999999,mu_bs[i,j]))
          mu_bs[i,j]<-PR_BS[i,j,Z[i]]
          
          for (s in 1:K){
            PR_BS[i,j,s]<-PsiBS.cut[j,s]*(1-ind[i,j])+ThetaBS[j,s]*ind[i,j]
          }
        }
      }
      
      ## cut the feedback from case model to FPR:
      for (j in 1:JBrS){
        for (s in 1:K){
          PsiBS.cut[j,s]<-cut(PsiBS[j,s])
        }
      }
      
      ## control BrS measurements
      for (i in (Nd+1):(Nd+Nu)){
        for (j in 1:JBrS){
          MBS[i,j]~dbern(mu_bs.bound[i,j])
          mu_bs.bound[i,j] <-max(0.000001,min(0.999999,mu_bs[i,j]))
          mu_bs[i,j]<-PsiBS[j,Z[i]]
        }
      }
      
      ## Case SS measure on pathogens that also have BrS:
      #for (i in 1:Nd){
      #for (j in 1:JSS){
      #  MSS[i,j] ~ dbern(mu_ss[i,j])
      #  mu_ss[i,j]<-ind[i,j]*thetaSS[j]+(1-ind[i,j])*psiSS[j]
      #}
      
      #for (j in 1:JSSonly){
      #	  ind[i,j+JBrS]<- equals(1,template[Icat[i],j+JBrS])
      #	  MSS.only[i,j] ~ dbern(mu_ss.only[i,j])
      #	  mu_ss.only[i,j]<-ind[i,j+JBrS]*thetaSS.only[j]
      #}
      #}
      
      for (i in 1:Nd){
        Z[i] ~ dcat(Eta[Icat[i],1:K])
        Icat[i] ~ dcat(pEti[1:Jcause])
      }
      
      ######################
      ##etiology prior
      ######################
      pEti[1:Jcause]~ddirch(alpha[])
      for (i in (Nd+1):(Nd+Nu)){
        Z[i]~dcat(Lambda[1:K])
      }
      
      ####################################
      ### stick-breaking specification 2
      ####################################
      Lambda0[1]<-r0[1]
      r0[K]<-1
      for(j in 2:K) {Lambda0[j]<-r0[j]*(1-r0[j-1])*Lambda0[j-1]/r0[j-1]}
      for(k in 1:K-1){
        r0[k]~dbeta(1,alphadp0)%_%I(0.000001,0.999999)
      }
      
      for (k in 1:K-1){Lambda[k]<-max(0.000001,min(0.999999,Lambda0[k]))}
      Lambda[K]<-1-sum(Lambda[1:(K-1)])
      
      for (s in 1:Jcause){
        Eta0[s,1]<-r1[s,1]
        r1[s,K]<-1
        for(j in 2:K) {Eta0[s,j]<-r1[s,j]*(1-r1[s,j-1])*Eta0[s,j-1]/r1[s,j-1]}
        for(k in 1:K-1){
          r1[s,k]~dbeta(1,alphadp0)%_%I(0.000001,0.999999)
        }
      }
      
      for (s in 1:Jcause){
        for (k in 1:K-1){Eta[s,k]<-max(0.000001,min(0.999999,Eta0[s,k]))}
        Eta[s,K]<-1-sum(Eta[s,1:(K-1)])
      }
      
      alphadp0~dgamma(.25,.25)%_%I(0.001,20)
      
      #########################
      ## priors on TPR and FPR:
      #########################
      
      for (j in 1:JBrS){
        for (s in 1:K){
          PsiBS[j,s]~dbeta(1,1)
          #ThetaBS[j,s]~dbeta(1,1)
          ThetaBS[j,s]~dbeta(alphaB[j],betaB[j])
        }
        ThetaBS.marg[j]<-inprod2(ThetaBS[j,1:K],Eta[j,1:K])
        PsiBS.marg[j]<-inprod2(PsiBS[j,1:K],Lambda[1:K])
        
        for (l in 1:Jcause){
          PsiBS.case[j,l]<-inprod2(PsiBS[j,1:K],Eta[l,1:K])
          # for calculating predicted positive rate from cases:
          # ThetaBS.case[j,l]<-inprod2(ThetaBS[j,1:K],Eta[l,1:K])
        }
      }
      
      # silver-standard measurement characteristics:
      #for (j in 1:JSS){
      #    thetaSS[j]~dbeta(alphaS[j],betaS[j])
      #    psiSS[j]<-0
      #}
      
      #for (j in 1:JSSonly){
      #	thetaSS.only[j]~dbeta(alphaS.only[j],betaS.only[j])
      #}
    
    }
    
    #-----------------END OF MODELS----------------------------------#
    
    #
    # run the model:
    #
    if (mcmc_options$ppd==TRUE){
      model_func         <- model_NoReg_BrS_nplcm_ppd
      model_bugfile_name <- "model_NoReg_BrS_nplcm_ppd.bug"
      #file.show(filename)
      # gs <- mybugs("model_NoReg_BrS_plcm_ppd.bug")
    } else {
      model_func         <- model_NoReg_BrS_nplcm
      model_bugfile_name <- "model_NoReg_BrS_nplcm.bug"
      #gs <- mybugs("model_NoReg_BrS_plcm.bug")
    }
    
    filename <- file.path(mcmc_options$bugsmodel.dir, model_bugfile_name)
    R2WinBUGS::write.model(model_func, filename)
    gs <- mybugs(model_bugfile_name)
}
