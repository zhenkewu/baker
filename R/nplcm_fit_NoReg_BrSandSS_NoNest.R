#' Fit nested partially-latent class model (low-level)
#'
#' Features:
#' \itemize{
#' \item no regression;
#' \item bronze- (BrS) and silver-standard (SS) measurements;
#' \item conditional independence;
#' \item all pathogens have BrS measure.
#' }
#'
#' @inheritParams nplcm
#' @importFrom R2WinBUGS bugs
#' @return WinBUGS fit results.
#'
#' @export

nplcm_fit_NoReg_BrSandSS_NoNest<-
  function(Mobs,Y,X,model_options,mcmc_options){
    # define generic function to call WinBUGS:
    call.bugs <- function(data, inits, parameters,m.file,
                          bugsmodel.dir = mcmc_options$bugsmodel.dir,
                          winbugs.dir   = mcmc_options$winbugs.dir,
                          nitermcmc     = mcmc_options$n.itermcmc,
                          nburnin       = mcmc_options$n.burnin,
                          nthin         = mcmc_options$n.thin,
                          nchains       = mcmc_options$n.chains,
                          dic      = FALSE,
                          is.debug = mcmc_options$debugstatus,
                          workd    = mcmc_options$result.folder,...) {

      m.file <- paste(bugsmodel.dir, m.file, sep="");
      f.tmp <- function() {
        ##winbugs
        gs <- bugs(data, inits, parameters,
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
    parsing <- assign_model(Mobs,Y,X,model_options)
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

    TPR_prior_list <- TPR_prior_set(model_options,Mobs,Y)
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

    # Tailor for each model:
    mybugs <- function(...){
      inits      <- function(){list(thetaBS = rbeta(JBrS,1,1),
                                    psiBS   = rbeta(JBrS,1,1))};
      data       <- c("Nd","Nu","JBrS","Jcause","alpha",
                      "template",
                      "MBS","JSS","MSS",
                      "alphaB","betaB","alphaS","betaS");

      if (mcmc_options$individual.pred==FALSE & mcmc_options$ppd==TRUE){
        parameters <- c("thetaBS","psiBS","pEti","thetaSS","MBS.new");

      } else if(mcmc_options$individual.pred==TRUE & mcmc_options$ppd==TRUE){
        parameters <- c("thetaBS","psiBS","pEti","thetaSS","Icat","MBS.new")

      } else if (mcmc_options$individual.pred==TRUE & mcmc_options$ppd==FALSE){
        parameters <- c("thetaBS","psiBS","pEti","thetaSS","Icat")

      } else if (mcmc_options$individual.pred==FALSE & mcmc_options$ppd==FALSE){
        parameters <- c("thetaBS","psiBS","pEti","thetaSS")

      }
      rst.bugs   <- call.bugs(data, inits, parameters,...);
      rst.bugs
    }

    if (mcmc_options$ppd==TRUE){
      gs <- mybugs("model_NoReg_BrSandSS_plcm_ppd.bug")
    } else {
      gs <- mybugs("model_NoReg_BrSandSS_plcm.bug")
    }

  }
