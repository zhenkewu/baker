if(getRversion() >= "2.15.1") utils::globalVariables(c("equals","Icat","thetaBS",
                                                       "Icat.new","psiBS"))



#' Fit nested partially-latent class model (low-level)
#'
#' Features:
#' \itemize{
#' \item no regression;
#' \item bronze- (BrS) measurements;
#' \item conditional independence;
#' \item all pathogens have BrS measurements.
#' }
#'
#' @inheritParams nplcm
#' @return WinBUGS fit results.
#'

nplcm_fit_NoReg_BrS_NoNest<-
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
    parsing <- assign_model(data_nplcm,model_options)
    Nd      <- sum(Y==1)
    Nu      <- sum(Y==0)

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

    if (parsing$measurement$SSonly){
      MSS.only.case <- Mobs$MSS[Y==1,(1:JSSonly)+JBrS]
      MSS.only <- as.matrix(MSS.only.case)
      alphaS.only <- TPR_prior_list$alphaS.only
      betaS.only  <- TPR_prior_list$betaS.only
    }

    mybugs <- function(...){
      inits      <- function(){list(thetaBS = rbeta(JBrS,1,1),
                                    psiBS   = rbeta(JBrS,1,1))};
      data       <- c("Nd","Nu","JBrS","Jcause","alpha",
                      "template",
                      "MBS",
                      "alphaB","betaB");

      if (mcmc_options$individual.pred==FALSE & mcmc_options$ppd==TRUE){
        parameters <- c("thetaBS","psiBS","pEti","MBS.new");

      } else if(mcmc_options$individual.pred==TRUE & mcmc_options$ppd==TRUE){
        parameters <- c("thetaBS","psiBS","pEti","Icat","MBS.new")

      } else if (mcmc_options$individual.pred==TRUE & mcmc_options$ppd==FALSE){
        parameters <- c("thetaBS","psiBS","pEti","Icat")

      } else if (mcmc_options$individual.pred==FALSE & mcmc_options$ppd==FALSE){
        parameters <- c("thetaBS","psiBS","pEti")

      }
      rst.bugs   <- call.bugs(data, inits, parameters,...);
      rst.bugs
    }

    
    #-----------------BEGINNING OF MODELS----------------------------------#
    
    #
    # write the .bug files into mcmc_options$bugsmodel.dir; could later set it equal to result.folder.
    #
    
    # the one with posterior predictive distribution:
    model_NoReg_BrS_plcm_ppd <- function(){#begin model
      ## 1) accommodates singletons, combos, and NoA;
      ## 2) check-bit to prevent case data informing FPR;
      
      # BrS measurements:
      for (k in 1:Nd){
        for (j in 1:JBrS){
          ind[k,j] <- equals(1,template[Icat[k],j])
          MBS[k,j] ~ dbern(mu_bs[k,j])
          mu_bs[k,j]<-ind[k,j]*thetaBS[j]+(1-ind[k,j])*psiBS.cut[j]
          
          #predictive checking:
          ind.new[k,j] <- equals(1,template[Icat.new[k],j])
          MBS.new[k,j] ~ dbern(mu_bs.new[k,j])
          mu_bs.new[k,j]<-ind.new[k,j]*thetaBS[j]+(1-ind.new[k,j])*psiBS.cut[j]
        }
      }
      
      for (k in (Nd+1):(Nd+Nu)){
        for (j in 1:JBrS){
          MBS[k,j] ~ dbern(mu_bs[k,j])
          mu_bs[k,j]<- psiBS[j]
          
          #predictive checking:
          MBS.new[k,j] ~ dbern(mu_bs[k,j])
        }
      }
      
      ## SS measurements on pathogens that also have BrS:
      #for (k in 1:Nd){
      #  for (j in 1:JSS){
      #	  MSS[k,j]   ~ dbern(mu_ss[k,j])
      #	  mu_ss[k,j] <- ind[k,j]*thetaSS[j]+(1-ind[k,j])*psiSS[j]
      #	}
      #for (j in 1:JSSonly){
      #  ind[k,j+JBrS]<- equals(1,template[Icat[k],j+JBrS])
      #  MSS.only[k,j] ~ dbern(mu_ss.only[k,j])
      #  mu_ss.only[k,j]<-ind[k,j+JBrS]*thetaSS.only[j]
      #}
      #}
      
      # priors
      for (k in 1:Nd){
        Icat[k] ~ dcat(pEti[1:Jcause])
        
        #predictive checking:
        Icat.new[k]~dcat(pEti[1:Jcause])
      }
      pEti[1:Jcause]~ddirch(alpha[])
      
      #for (k in (Nd+1):(Nd+Nu)){
      #  Icat[k] <- Jcause+1
      #}
      
      # bronze-standard measurement characteristics:
      for (j in 1:JBrS){
        thetaBS[j]~dbeta(alphaB[j],betaB[j])
        psiBS[j]~dbeta(1,1)
        psiBS.cut[j]<-cut(psiBS[j])
      }
      
      # silver-standard measurement characteristics:
      #for (j in 1:JSS){
      #	thetaSS[j]~dbeta(alphaS[j],betaS[j])
      #	psiSS[j]<-0
      #}
      
      #for (j in 1:JSSonly){
      #	thetaSS.only[j]~dbeta(alphaS.only[j],betaS.only[j])
      #}
      
    }#end of model
    # the one without ppd:
    model_NoReg_BrS_plcm <- function(){#begin model
        ## 1) accommodates singletons, combos, and NoA;
        ## 2) check-bit to prevent case data informing FPR;
        
        # BrS measurements:
        for (k in 1:Nd){
          for (j in 1:JBrS){
            ind[k,j] <- equals(1,template[Icat[k],j])
            MBS[k,j] ~ dbern(mu_bs[k,j])
            mu_bs[k,j]<-ind[k,j]*thetaBS[j]+(1-ind[k,j])*psiBS.cut[j]
            
            
          }
        }
        
        for (k in (Nd+1):(Nd+Nu)){
          for (j in 1:JBrS){
            MBS[k,j] ~ dbern(mu_bs[k,j])
            mu_bs[k,j]<- psiBS[j]
            
          }
        }
        
        ## SS measurements on pathogens that also have BrS:
        #for (k in 1:Nd){
        #  for (j in 1:JSS){
        #	  MSS[k,j]   ~ dbern(mu_ss[k,j])
        #	  mu_ss[k,j] <- ind[k,j]*thetaSS[j]+(1-ind[k,j])*psiSS[j]
        #	}
        #for (j in 1:JSSonly){
        #  ind[k,j+JBrS]<- equals(1,template[Icat[k],j+JBrS])
        #  MSS.only[k,j] ~ dbern(mu_ss.only[k,j])
        #  mu_ss.only[k,j]<-ind[k,j+JBrS]*thetaSS.only[j]
        #}
        #}
        
        # priors
        for (k in 1:Nd){
          Icat[k] ~ dcat(pEti[1:Jcause])
          
        }
        pEti[1:Jcause]~ddirch(alpha[])
        
        #for (k in (Nd+1):(Nd+Nu)){
        #  Icat[k] <- Jcause+1
        #}
        
        # bronze-standard measurement characteristics:
        for (j in 1:JBrS){
          thetaBS[j]~dbeta(alphaB[j],betaB[j])
          psiBS[j]~dbeta(1,1)
          psiBS.cut[j]<-cut(psiBS[j])
        }
        
        # silver-standard measurement characteristics:
        #for (j in 1:JSS){
        #	thetaSS[j]~dbeta(alphaS[j],betaS[j])
        #	psiSS[j]<-0
        #}
        
        #for (j in 1:JSSonly){
        #	thetaSS.only[j]~dbeta(alphaS.only[j],betaS.only[j])
        #}
        
      }#end of model
    
    #-----------------END OF MODELS----------------------------------#
    
    #
    # run the model:
    #
    if (mcmc_options$ppd==TRUE){
      model_func         <- model_NoReg_BrS_plcm_ppd
      model_bugfile_name <- "model_NoReg_BrS_plcm_ppd.bug"
      #file.show(filename)
     # gs <- mybugs("model_NoReg_BrS_plcm_ppd.bug")
    } else {
      model_func         <- model_NoReg_BrS_plcm
      model_bugfile_name <- "model_NoReg_BrS_plcm.bug"
      #gs <- mybugs("model_NoReg_BrS_plcm.bug")
    }
    
    filename <- file.path(mcmc_options$bugsmodel.dir, model_bugfile_name)
    R2WinBUGS::write.model(model_func, filename)
    gs <- mybugs(model_bugfile_name)
  }
