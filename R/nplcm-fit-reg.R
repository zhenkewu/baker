#' (Low-level function) Fitting nested partially-latent class models with
#' stratification/regression. [needs revision!!]
#'
#' @inheritParams nplcm
#' @return A WinBUGS result, fitted by function \code{bugs()} from the R2WinBUGS package.
#' @export
#'
nplcm_fit_reg<-function(data_nplcm,model_options,mcmc_options){#BEGIN function
  Mobs <- data_nplcm$Mobs
  Y    <- data_nplcm$Y
  X    <- data_nplcm$X
  #define generic function to call WinBUGS:
  call.bugs <- function(data, inits, parameters,m.file,
                        bugsmodel.dir = mcmc_options$bugsmodel.dir,
                        winbugs.dir   = mcmc_options$winbugs.dir,
                        nitermcmc     = mcmc_options$n.itermcmc,
                        nburnin       = mcmc_options$n.burnin,
                        nthin         = mcmc_options$n.thin,
                        nchains       = mcmc_options$n.chains,
                        dic = FALSE, is.debug = mcmc_options$debugstatus,
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

  ## for discrete stratified/regression analysis
  X_group        <- X[model_options$X_reg]
  # dichotomize age variable:
  if ("AGECAT" %in% model_options$X_reg){
   X_group$AGECAT <- as.numeric(X_group$AGECAT > 1)+1
  }
  X_group$group_names <- apply(X_group,1,paste,collapse="&")
  X_group$ID     <- 1:nrow(X_group)

  form_agg <- as.formula(paste0("cbind(group_names,ID)~",
                                paste(model_options$X_reg,collapse="+")))

  grouping <- aggregate(form_agg,X_group,identity)

  ## temporary code to get the count of observations in each group:
  #form_agg2 <- as.formula(paste0("cbind(group_names,ID)~",
  #                              paste(c("Y",model_options$X_reg),collapse="+")))
  #aggregate(form_agg2,X_group,length)

  group_nm <- lapply(grouping$group_names,function(v) unique(as.character(v)))
  names(group_nm) <- 1:length(group_nm)
  X_group$grp <- rep(NA,nrow(X_group))

  for (l in seq_along(group_nm)){
    X_group$grp[unfactor(grouping$ID[[l]])] <- l
  }

  group <- X_group$grp
  N_vec <- table(group)
  N_grp <- length(N_vec)

  if (is.null(model_options$SSonly)|| model_options$SSonly==FALSE){
    ##if all pathogens have NPPCR data:--------------------------
    #check sources of measurement data:
    data_source        <- c("yes","no")[is.na(Mobs)+1]
    names(data_source) <- names(Mobs)
    cat("Available measurements: \n")
    print(data_source)

    #compatibility checking:
    if (length(model_options$M_use)!=length(model_options$TPR_prior)){
      stop("The number of measurement source(s) is different from
           the number of TPR prior option!
           Make them equal, and match with order!")
    }

    #some data preparation:
    Nd <- sum(Y==1)
    Nu <- sum(Y==0)

    model_data_source <- rep(NA,3)
    names(model_data_source) <- c("MBS","MSS","MGS")
    model_data_source[1] <- c("no","yes")["BrS"%in%model_options$M_use+1]
    model_data_source[2] <- c("no","yes")["SS"%in%model_options$M_use+1]
    model_data_source[3] <- c("no","yes")["GS"%in%model_options$M_use+1]

    cat("Actual measurements used in the model: \n")
    print(model_data_source)

    cat("True positive rate (TPR) prior(s) for ",
        c("MBS","MSS","MGS")[model_data_source=="yes"],"\n",
        " is(are respectively) ", model_options$TPR_prior,"\n")

    allowed_list <- model_options$allowed_list
    pathogen_list <-model_options$pathogen_list

    Jfull         <- length(pathogen_list)
    Jallowed      <- length(allowed_list)
    template <-  as.matrix(rbind(symb2I(allowed_list,
                                        pathogen_list),rep(0,Jfull)))

    # BEGIN fit model -------------------------------------------------------------
    if (model_options$k_subclass ==1){#BEGIN plcm:
      cat("Number of subclasses: ", model_options$k_subclass,"\n")
      #plcm: conditional independence model
      #plcm-BrS only:
      if (model_data_source[1]=="yes" &
            model_data_source[2]=="no" &
            model_data_source[3]=="no"){

        MBS.case <- Mobs$MBS[Y==1,]
        MBS.ctrl <- Mobs$MBS[Y==0,]
        MBS      <- as.matrix(rbind(MBS.case,MBS.ctrl))


        alpha <- set_prior_eti(model_options)

        TPR_prior_list <- set_prior_tpr(model_options,data_nplcm)

        alphaB <- TPR_prior_list$alphaB
        betaB <- TPR_prior_list$betaB

        mybugs <- function(...){
          inits      <- function(){list(thetaBS = rbeta(Jfull,1,1),
                                        psiBS   = rbeta(Jfull,1,1))};
          data       <- c("Nd","Nu","Jfull","Jallowed","alpha",
                          "template","MBS","alphaB","betaB");

          if (mcmc_options$individual.pred==FALSE & mcmc_options$ppd==TRUE){
            parameters <- c("thetaBS","psiBS","pEti","MBS.new");

          } else if(mcmc_options$individual.pred==TRUE & mcmc_options$ppd==TRUE){
            parameters <- c("thetaBS","psiBS","pEti","Icat","MBS.new");

          } else if (mcmc_options$individual.pred==TRUE & mcmc_options$ppd==FALSE){
            parameters <- c("thetaBS","psiBS","pEti","Icat");

          } else if (mcmc_options$individual.pred==FALSE & mcmc_options$ppd==FALSE){
            parameters <- c("thetaBS","psiBS","pEti");

          }
          rst.bugs   <- call.bugs(data, inits, parameters,...);
          rst.bugs
        }

        if (mcmc_options$ppd==TRUE){
          gs <- mybugs("model_plcm_brsonly_ppd.bug")
        } else {
          gs <- mybugs("model_plcm_brsonly.bug")
        }
      }

      ##plcm-BrS + SS:
      if (model_data_source[1]=="yes" &
            model_data_source[2]=="yes" &
              model_data_source[3]=="no"){

        MBS.case <- Mobs$MBS[Y==1,]
        MBS.ctrl <- Mobs$MBS[Y==0,]
        MBS      <- as.matrix(rbind(MBS.case,MBS.ctrl))

        MSS.case <- Mobs$MSS[Y==1,1:Jfull]
        MSS.case <- as.matrix(MSS.case)

        SS_index <- which(colMeans(is.na(MSS.case))<0.9)#.9 is arbitrary; any number <1 will work.
        JSS      <- length(SS_index)
        MSS      <- MSS.case[,SS_index]

        alpha <- set_prior_eti(model_options)

        TPR_prior_list <- set_prior_tpr(model_options,data_nplcm)

        alphaB <- TPR_prior_list$alphaB
        betaB <- TPR_prior_list$betaB
        alphaS <- TPR_prior_list$alphaS
        betaS <- TPR_prior_list$betaS

        mybugs <- function(...){
          inits      <- function(){list(thetaBS = rbeta(Jfull,1,1),
                                        psiBS   = matrix(rbeta(Jfull*N_grp,1,1),
                                                     nrow=N_grp,ncol=Jfull))};
          data       <- c("Nd","Nu","Jfull","Jallowed","alpha","template",
                          "MBS","JSS","MSS",
                          "N_grp","group",
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
          stop("== Not yet implemented. Please contact maintainer. Thanks.")
          #gs <- mybugs("model_plcm_ppd.bug")
        } else {
          gs <- mybugs("model_Reg_BrSandSS_plcm.bug")
        }
      }

    }else{#END plcm, BEGIN nplcm:
      #nplcm: conditional dependence model
      cat("Number of subclasses: ", model_options$k_subclass,"\n")
      # nplcm: BrS only:
      if (model_data_source[1]=="yes" &
            model_data_source[2]=="no" &
              model_data_source[3]=="no"){

        MBS.case <- Mobs$MBS[Y==1,]
        MBS.ctrl <- Mobs$MBS[Y==0,]
        MBS      <- as.matrix(rbind(MBS.case,MBS.ctrl))

        K        <- model_options$k_subclass

        alpha <- set_prior_eti(model_options)

        TPR_prior_list <- set_prior_tpr(model_options,data_nplcm)

        alphaB <- TPR_prior_list$alphaB
        betaB <- TPR_prior_list$betaB

        mybugs <- function(...){
          inits      <- function(){list(pEti = rep(1/Jallowed,Jallowed),
                                        r0 = c(rep(.5,K-1),NA),
                                        r1 = cbind(matrix(rep(.5,Jfull*(K-1)),
                                                          nrow=Jfull,ncol=K-1),
                                                   rep(NA,Jfull)),
                                        alphadp0 = 1)};
          data       <- c("Nd","Nu","Jfull","Jallowed",
                          "alpha","template","K",
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

        if (mcmc_options$ppd==TRUE){
          gs <- mybugs("model_nplcm_brsonly_ppd.bug")
        } else {
          gs <- mybugs("model_nplcm_brsonly.bug")
        }
      }


      ## nplcm: BrS + SS:
      if (model_data_source[1]=="yes" &
            model_data_source[2]=="yes" &
              model_data_source[3]=="no"){

        MBS.case <- Mobs$MBS[Y==1,]
        MBS.ctrl <- Mobs$MBS[Y==0,]
        MBS      <- as.matrix(rbind(MBS.case,MBS.ctrl))

        MSS.case <- Mobs$MSS[Y==1,1:Jfull]
        MSS.case <- as.matrix(MSS.case)

        SS_index <- which(colMeans(is.na(MSS.case))<0.9)#.9 is arbitrary; any number <1 will work.
        JSS      <- length(SS_index)
        MSS      <- MSS.case[,SS_index]

        K        <- model_options$k_subclass

        alpha <- set_prior_eti(model_options)

        TPR_prior_list <- set_prior_tpr(model_options,data_nplcm)

        alphaB <- TPR_prior_list$alphaB
        betaB <- TPR_prior_list$betaB
        alphaS <- TPR_prior_list$alphaS
        betaS <- TPR_prior_list$betaS

        mybugs <- function(...){
          inits      <- function(){list(pEti = rep(1/Jallowed,Jallowed),
                                        r0 = c(rep(.5,K-1),NA),
                                        r1 = cbind(matrix(rep(.5,Jfull*(K-1)),
                                                          nrow=Jfull,ncol=K-1),
                                                   rep(NA,Jfull)),
                                        alphadp0 = 1)};
          data       <- c("Nd","Nu","Jfull","Jallowed",
                          "alpha","template","K",
                          "JSS","MSS",
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
          gs <- mybugs("model_nplcm_ppd.bug")
        } else {
          gs <- mybugs("model_nplcm.bug")
        }
      }
    }
    #END fitting model---------------------------------------------------------
    }else{#BEGIN SS only model:
      ##-------------- if SS only pathogens are available (some pathogens have no NPPCR data):
      #check sources of measurement data:
      data_source        <- c("yes","no")[is.na(Mobs)+1]
      names(data_source) <- names(Mobs)
      cat("Available measurements: \n")
      print(data_source)

      #compatibility checking:
      if (length(model_options$M_use)!=length(model_options$TPR_prior)){
        stop("The number of measurement source(s) is different from
             the number of TPR prior option!
             Make them equal, and match them in the same order!")
      }

      #some data preparation:
      Nd <- sum(Y==1)
      Nu <- sum(Y==0)

      model_data_source <- rep(NA,3)
      names(model_data_source) <- c("MBS","MSS","MGS")
      model_data_source[1] <- c("no","yes")["BrS"%in%model_options$M_use+1]
      model_data_source[2] <- c("no","yes")["SS"%in%model_options$M_use+1]
      model_data_source[3] <- c("no","yes")["GS"%in%model_options$M_use+1]

      cat("Actual measurements used in the model: \n")
      print(model_data_source)

      cat("True positive rate (TPR) prior(s) for ",
          c("MBS","MSS","MGS")[model_data_source=="yes"],"\n",
          " is(are respectively) ", model_options$TPR_prior,"\n")

      allowed_list <- model_options$allowed_list
      pathogen_list <-model_options$pathogen_list

      Jfull         <- length(pathogen_list)
      Jallowed      <- length(allowed_list)


      #SS only pathogen parameters:
      pathogen_SSonly_list <- model_options$pathogen_SSonly_list
      JSS_only     <- length(pathogen_SSonly_list)

      template <-  as.matrix(rbind(diag(1,Jfull+JSS_only),rep(0,Jfull+JSS_only)))

      # BEGIN fit model -------------------------------------------------------------
      if (model_options$k_subclass ==1){#BEGIN plcm model:
                #plcm: conditional independence model
                cat("Number of subclasses: ", model_options$k_subclass,"\n")

                #plcm-BrS only (not allowed, will output error.):
                if (model_data_source[1]=="yes" &
                      model_data_source[2]=="no" &
                        model_data_source[3]=="no"){
                  stop("In model_options, you specified pathogens with SS only.
                       Please include SS measurement data.")
                }

                # plcm- SS data only model fit: to be written if necessary.

                ##plcm-BrS + SS:
                if (model_data_source[1]=="yes" &
                      model_data_source[2]=="yes" &
                        model_data_source[3]=="no"){

                      MBS.case <- Mobs$MBS[Y==1,]
                      MBS.ctrl <- Mobs$MBS[Y==0,]
                      MBS      <- as.matrix(rbind(MBS.case,MBS.ctrl))

                      MSS.case <- Mobs$MSS[Y==1,1:Jfull]
                      MSS.case <- as.matrix(MSS.case)

                      MSS.only.case <- Mobs$MSS[Y==1,(1:JSS_only)+Jfull]

                      SS_index <- which(colMeans(is.na(MSS.case))<0.9)#.9 is arbitrary; any number <1 will work.
                      JSS      <- length(SS_index)
                      MSS      <- MSS.case[,SS_index]
                      MSS.only <- as.matrix(MSS.only.case)

                      alpha <- set_prior_eti(model_options)

                      TPR_prior_list <- set_prior_tpr(model_options,data_nplcm)

                      alphaB <- TPR_prior_list$alphaB
                      betaB <- TPR_prior_list$betaB
                      alphaS <- TPR_prior_list$alphaS
                      betaS <- TPR_prior_list$betaS
                      alphaS.only <- TPR_prior_list$alphaS.only
                      betaS.only  <- TPR_prior_list$betaS.only


                      mybugs <- function(...){
                        inits      <- function(){list(thetaBS = rbeta(Jfull,1,1),
                                                      psiBS   = matrix(rbeta(Jfull*N_grp,1,1),
                                                                       nrow=N_grp,ncol=Jfull))};
                        data       <- c("Nd","Nu","Jfull","Jallowed","alpha","template",
                                        "MBS","JSS","MSS",
                                        "MSS.only","JSS_only","alphaS.only","betaS.only",
                                        "N_grp","group",
                                        "alphaB","betaB","alphaS","betaS");

                        if (mcmc_options$individual.pred==FALSE & mcmc_options$ppd==TRUE){
                          parameters <- c("thetaBS","psiBS","pEti","thetaSS","thetaSS.only","MBS.new");

                        } else if(mcmc_options$individual.pred==TRUE & mcmc_options$ppd==TRUE){
                          parameters <- c("thetaBS","psiBS","pEti","thetaSS","thetaSS.only","Icat","MBS.new")

                        } else if (mcmc_options$individual.pred==TRUE & mcmc_options$ppd==FALSE){
                          parameters <- c("thetaBS","psiBS","pEti","thetaSS","thetaSS.only","Icat")

                        } else if (mcmc_options$individual.pred==FALSE & mcmc_options$ppd==FALSE){
                          parameters <- c("thetaBS","psiBS","pEti","thetaSS","thetaSS.only")

                        }
                        rst.bugs   <- call.bugs(data, inits, parameters,...);
                        rst.bugs
                      }

                      if (mcmc_options$ppd==TRUE){
                        stop("== Not yet implemented. Please contact maintainer. Thanks.==")
                        #gs <- mybugs("model_plcm_SSonly_ppd.bug")
                      } else {
                        gs <- mybugs("model_Reg_BrSandSS_SSonly_plcm.bug")
                      }
                }
       }#END plcm model.
   }#END silver only model
  # return the WinBUGS fitted results:
  return(gs)
}#END function






















