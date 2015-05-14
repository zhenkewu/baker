#' Fit nested partially-latent class model (low-level)
#'
#' Features:
#' \itemize{
#' \item no regression;
#' \item bronze- (BrS)and silver-standard (SS) measurements;
#' \item conditional dependence;
#' \item all pathogens have BrS measurements.
#' }
#'
#' @inheritParams nplcm
#' @return WinBUGS fit results.
#'
#' @export

nplcm_fit_NoReg_BrSandSS_Nest <-
  function(data_nplcm,model_options,mcmc_options){
    
    Mobs <- data_nplcm$Mobs
    Y    <- data_nplcm$Y
    X    <- data_nplcm$X
    
    
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

      m.file <- paste(bugsmodel.dir, m.file, sep="");
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

    MSS.case <- Mobs$MSS[Y==1,1:JBrS]
    MSS.case <- as.matrix(MSS.case)

    SS_index <- which(colMeans(is.na(MSS.case))<0.9)#.9 is arbitrary; any number <1 will work.
    JSS      <- length(SS_index)
    MSS      <- MSS.case[,SS_index]

    # set priors:
    alpha          <- eti_prior_set(model_options)

    TPR_prior_list <- TPR_prior_set(model_options,data_nplcm)
    alphaB      <- TPR_prior_list$alphaB
    betaB       <- TPR_prior_list$betaB
    alphaS      <- TPR_prior_list$alphaS
    betaS       <- TPR_prior_list$betaS
    if (parsing$measurement$SSonly){
      MSS.only.case <- Mobs$MSS[Y==1,(1:JSSonly)+JBrS]
      MSS.only <- as.matrix(MSS.only.case)
      alphaS.only <- TPR_prior_list$alphaS.only
      betaS.only  <- TPR_prior_list$betaS.only
    }


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
                      "JSS","MSS",
                      #"MSS.only","JSSonly","alphaS.only","betaS.only",
                      "MBS","alphaB","betaB","alphaS","betaS");

      if (mcmc_options$individual.pred==FALSE & mcmc_options$ppd==TRUE){
        parameters <- c("pEti","Lambda","Eta","alphadp0","MBS.new",
                        "ThetaBS.marg","PsiBS.marg","PsiBS.case",
                        "ThetaBS","PsiBS","thetaSS")

      } else if(mcmc_options$individual.pred==TRUE & mcmc_options$ppd==TRUE){
        parameters <- c("pEti","Lambda","Eta","alphadp0","Icat","MBS.new",
                        "ThetaBS.marg","PsiBS.marg","PsiBS.case",
                        "ThetaBS","PsiBS","thetaSS")

      } else if (mcmc_options$individual.pred==TRUE & mcmc_options$ppd==FALSE){
        parameters <- c("pEti","Lambda","Eta","alphadp0","Icat",
                        "ThetaBS.marg","PsiBS.marg","PsiBS.case",
                        "ThetaBS","PsiBS","thetaSS")

      } else if (mcmc_options$individual.pred==FALSE & mcmc_options$ppd==FALSE){
        parameters <- c("pEti","Lambda","Eta","alphadp0",
                        "ThetaBS.marg","PsiBS.marg","PsiBS.case",
                        "ThetaBS","PsiBS","thetaSS")

      }
      rst.bugs   <- call.bugs(data, inits, parameters,...);
      rst.bugs
    }

    if (mcmc_options$ppd==TRUE){
      gs <- mybugs("model_NoReg_BrSandSS_nplcm_ppd.bug")
    } else {
      gs <- mybugs("model_NoReg_BrSandSS_nplcm.bug")
    }

  }
