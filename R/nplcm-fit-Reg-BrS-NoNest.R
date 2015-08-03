if(getRversion() >= "2.15.1") utils::globalVariables(c("set_prior_tpr","set_prior_eti"))

#' Fit nested partially-latent class model with regression (low-level)
#'
#' \code{nplcm_fit_Reg_BrS_NoNest} use design matrix to perform etiology 
#' and/or false positive regressions with following features:
#' \itemize{
#' \item General regression;
#' \item bronze- (BrS) measurements;
#' \item conditional independence;
#' \item all pathogens have BrS measurements.
#' }
#'
#' @inheritParams nplcm
#' @return WinBUGS fit results.
#'
#' @export

nplcm_fit_Reg_BrS_NoNest <- function(data_nplcm,model_options,mcmc_options){
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
    # parsing <- assign_model(data_nplcm,model_options)
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
    
    
    # get design matrices:
    dm_FPR <- try(model.matrix(model_options$X_reg_FPR,
                            data.frame(data_nplcm$X,Y = data_nplcm$Y)))
    dm_Eti <- try(model.matrix(model_options$X_reg_Eti,
                            data.frame(data_nplcm$X, Y = data_nplcm$Y)))

    num.knots.FPR <- length(grep("random",colnames(dm_FPR)))
    num.knots.Eti <- length(grep("dm_Rdate_Eti",colnames(dm_Eti)))
    ncol_dm_FPR   <- ncol(dm_FPR)
    ncol_dm_Eti   <- ncol(dm_Eti)
    
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
    
    mybugs <- function(...){
      inits      <- function(){list(beta = matrix(0,nrow=JBrS,ncol=ncol(dm_FPR)-num.knots.FPR),
                                    b = matrix(0,nrow=JBrS,ncol=num.knots.FPR),
                                    taub = rep(0.01,JBrS),
                                    betaEti = rbind(matrix(0,nrow=JBrS-1,ncol=ncol(dm_Eti)),
                                                    rep(NA,ncol(dm_Eti)))
                                    )};
      data       <- c("Nd","Nu","JBrS","Jcause","ncol_dm_FPR","ncol_dm_Eti",
                      "template",
                      "MBS",
                      "dm_Eti",#"num.knots.Eti",
                      "dm_FPR","num.knots.FPR",
                      "alphaB","betaB");
      
      if (mcmc_options$individual.pred==FALSE & mcmc_options$ppd==TRUE){
        stop("==Not done.==")
        
      } else if(mcmc_options$individual.pred==TRUE & mcmc_options$ppd==TRUE){
        stop("==Not done.==")
        
      } else if (mcmc_options$individual.pred==TRUE & mcmc_options$ppd==FALSE){
        stop("==Not done.==")
        
      } else if (mcmc_options$individual.pred==FALSE & mcmc_options$ppd==FALSE){
        parameters <- c("taub","sigmab","beta","b","betaEti","Icat","pEti",
                        "tpr","mu","fpr.case","fpr.ctrl")
        
      }
      rst.bugs   <- call.bugs(data, inits, parameters,...);
      rst.bugs
    }
    
    
    #-----------------BEGINNING OF MODELS----------------------------------#
    
    #
    # write the .bug files into mcmc_options$bugsmodel.dir; 
    # could later set it equal to result.folder.
    #
    
    # the one with posterior predictive distribution:
    model_Reg_BrS_plcm_ppd <- function(){
    
    }
    # the one without ppd:
    model_Reg_BrS_plcm <- function(){
             chunk1 <-  paste0("for (j in 1:JBrS){
                #cases
                for (k in 1:Nd){
                  ind[k,j] <- equals(1,template[Icat[k],j])
                  MBS[k,j] ~ dbern(mu[k,j])
                  mu[k,j]<-ind[k,j]*tpr[j]+(1-ind[k,j])*fpr.case[k,j] 		
                  logit(fpr.case[k,j])<-m[k,j]
                }
                
                #controls
                for (k in (Nd+1):(Nd+Nu)){
                  ind[k,j]<-0
                  MBS[k,j]~dbern(fpr.ctrl[k-Nd,j]) 
                  logit(fpr.ctrl[k-Nd,j])<-m[k,j]
                }
                
                for (k in 1:(Nd+Nu)){
                  m[k,j]<-mfe[k,j]+mre[k,j]
                  mfe[k,j]<- ",create_bugs_regressor_Eti(ncol(dm_FPR)-num.knots.FPR,"dm_FPR","beta"),"\n 
                  mre[k,j] <- ",create_bugs_regressor_FPR(num.knots.FPR),"\n")
             
             chunk2 <-paste0("}
              }
              
              #############################################
              ##  priors
              #############################################
              for (j in 1:JBrS){
                for (o in 1:num.knots.FPR){
                  b[j,o]~dnorm(0,taub[j])
                }
                
                for (l in 1:(ncol_dm_FPR-num.knots.FPR)){
                  beta[j,l]~dnorm(0,1.0E-2)
                }
                
                taub[j]~dgamma(1.0E-2,1.0E-2)
                
                sigmab[j]<-1/sqrt(taub[j])
              }
              
              
              ## etiologic process prior
              for (k in 1:Nd){
                Icat[k] ~ dcat(pEti[k,1:JBrS])
                for (j in 1:Jcause){     
                  pEti[k,j]     <- phi[k,j]/sum(phi[k,])
                  log(phi[k,j]) <- mEti[k,j]
                  mEti[k,j] <- ",create_bugs_regressor_Eti(ncol(dm_Eti)),"\n")
             chunk3 <-paste0("}
              }
              
              for (j in 1:(Jcause-1)){
                for (l in 1:(ncol_dm_Eti)){
                  betaEti[j,l]~dnorm(0,0.1) #
                  #betaEti[j,l]~dnorm(0,tauEti[l]) # to make betaEti's to be similar
                }
              }
              
              #	for (l in 1:(num.knots.Eti+1)){
              #		tauEti[l]~dgamma(0.1,0.1)
              #		sigmaEti[l]<-1/sqrt(tauEti[l])
              #	}
              
              for (l in 1:(ncol_dm_Eti)){
                betaEti[JBrS,l]<-0
              }
              
              ## priors for tpr and fpr.case
              for (j in 1:JBrS){
                # Delta[j]~dnorm(0,0.1)
                tpr[j]~dbeta(alphaB[j],betaB[j])
              }")
             
             
            paste0("model{#BEGIN OF MODEL:\n",
                   chunk1,"\n",chunk2,"\n",chunk3,"\n",
                   "}#END OF MODEL.")
    }
    
    #-----------------END OF MODELS----------------------------------#
    
    #
    # run the model:
    #
    if (mcmc_options$ppd==TRUE){
      stop("==Posterior predictive checking not implemented for this model. Please check back later or contact the maintainer. Thanks.==")
      model_func         <- model_Reg_BrS_plcm_ppd
      model_bugfile_name <- "model_Reg_BrS_plcm_ppd.bug"
      #file.show(filename)
    } else {
      model_func         <- model_Reg_BrS_plcm()
      model_bugfile_name <- "model_Reg_BrS_plcm.bug"
    }
    
    filename <- file.path(mcmc_options$bugsmodel.dir, model_bugfile_name)
    writeLines(model_func,filename)
    gs <- mybugs(model_bugfile_name)
  }
