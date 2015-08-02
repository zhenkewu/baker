if(getRversion() >= "2.15.1") utils::globalVariables(c("equals","Icat","thetaBS",
                                                       "Icat.new","psiBS"))
#' Fit nested partially-latent class model (low-level)
#'
#' Features:
#' \itemize{
#' \item no regression;
#' \item no nested
#' }
#'
#' @inheritParams nplcm
#' @return WinBUGS fit results.
#' 
#' @seealso \link{write_model_NoReg_plcm} for writing .bug model file; 
#' 
#' @export
nplcm_fit_NoReg_NoNest<-
  function(data_nplcm,model_options,mcmc_options){
    # Record the settings of current analysis:
    cat("==Results stored in: ==","\n",mcmc_options$result.folder)
    #model_options:
    dput(model_options,file.path(mcmc_options$result.folder,"model_options.txt"))
    #mcmc_options:
    dput(mcmc_options,file.path(mcmc_options$result.folder,"mcmc_options.txt"))
    
    # read in data:
    Mobs <- data_nplcm$Mobs
    Y    <- data_nplcm$Y
    
    # read in options:
    likelihood       <- model_options$likelihood
    use_measurements <- model_options$use_measurements
    prior            <- model_options
    
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
    
    # get sample sizes:
    Nd      <- sum(Y==1)
    Nu      <- sum(Y==0)

    # get lengths of vectors:
    cause_list  <- likelihood$cause_list
    Jcause      <- length(cause_list)

    in_data = in_init = out_parameter <- NULL
    #
    # 1. BrS measurement data: 
    #
    if ("BrS" %in% use_measurements){
        JBrS_list  <- lapply(Mobs$MBS,ncol)
        
        # mapping template (by `make_template` function):
        patho_BrS_list <- lapply(Mobs$MBS,colnames)
        template_BrS_list <- lapply(patho_BrS_list,make_template,cause_list)
        
        MBS.case_list <- lapply(Mobs$MBS,"[",which(Y==1),TRUE,drop=FALSE)
        MBS.ctrl_list <- lapply(Mobs$MBS,"[",which(Y==0),TRUE,drop=FALSE)
        MBS_list <- list()
        for (i in seq_along(MBS.case_list)){
          MBS_list[[i]]      <- rbind(MBS.case_list[[i]],MBS.ctrl_list[[i]])
        }
        names(MBS_list) <- names(MBS.case_list)
        
        single_column_MBS <- which(lapply(MBS_list,ncol)==1)
        for(i in seq_along(JBrS_list)){
          assign(paste("JBrS", i, sep = "_"), JBrS_list[[i]])    
          assign(paste("MBS", i, sep = "_"), as.matrix_or_vec(MBS_list[[i]])) 
          assign(paste("templateBS", i, sep = "_"), as.matrix_or_vec(template_BrS_list[[i]]))   
        }
      
        # set BrS measurement priors: (need to add informative specification!)
        # hyperparameter for sensitivity (can add for specificity if necessary): 
        for(i in seq_along(JBrS_list)){
          assign(paste("alphaB", i, sep = "_"), rep(1,JBrS_list[[i]]))    
          assign(paste("betaB", i, sep = "_"),  rep(1,JBrS_list[[i]]))    
        }
        
        if (length(single_column_MBS)==0){
            # summarize into one name (for all measurements):
            in_data       <- c(in_data,"Nd","Nu","Jcause","alpha",
                               paste("JBrS",1:length(JBrS_list),sep="_"),
                               paste("MBS",1:length(JBrS_list),sep="_"),
                               paste("templateBS",1:length(JBrS_list),sep="_"),
                               paste("alphaB",1:length(JBrS_list),sep="_"),
                               paste("betaB",1:length(JBrS_list),sep="_")
                               )
        } else {
            # summarize into one name (for all measurements):
            in_data       <- c(in_data,"Nd","Nu","Jcause","alpha",
                               paste("JBrS",1:length(JBrS_list),sep="_")[-single_column_MBS],
                               paste("MBS",1:length(JBrS_list),sep="_"),
                               paste("templateBS",1:length(JBrS_list),sep="_"),
                               paste("alphaB",1:length(JBrS_list),sep="_"),
                               paste("betaB",1:length(JBrS_list),sep="_")
            )
        }
        out_parameter <- c(out_parameter, paste("thetaBS", seq_along(JBrS_list), sep = "_"),
                           paste("psiBS", seq_along(JBrS_list), sep = "_"),
                           "pEti")
    }
    
    #
    # 2. SS measurement data: 
    #
    if ("SS" %in% use_measurements){
        JSS_list   <- lapply(Mobs$MSS,ncol)
        
        # mapping template (by `make_template` function):
        patho_SS_list <- lapply(Mobs$MSS,colnames)
        template_SS_list <- lapply(patho_SS_list,make_template,cause_list)
        
        MSS_list <- lapply(Mobs$MSS,"[",which(Y==1),TRUE,drop=FALSE)
        
        single_column_MSS <- which(lapply(MSS_list,ncol)==1)
        
        for(i in seq_along(JSS_list)){
          assign(paste("JSS", i, sep = "_"), JSS_list[[i]])    
          assign(paste("MSS", i, sep = "_"), as.matrix_or_vec(MSS_list[[i]])) 
          assign(paste("templateSS", i, sep = "_"), as.matrix_or_vec(template_SS_list[[i]]))   
        }
        
        # set SS measurement priors: (need to add informative specification!)
        # hyper parameters for sensitivity:
        for(i in seq_along(JSS_list)){
          assign(paste("alphaS", i, sep = "_"), rep(1,JSS_list[[i]]))    
          assign(paste("betaS", i, sep = "_"),  rep(1,JSS_list[[i]]))    
        }
        
        if (length(single_column_MSS)==0){
          # summarize into one name (for all measurements):
          in_data       <- unique(c(in_data,"Nd","Jcause","alpha",
                             paste("JSS",1:length(JSS_list),sep="_"),
                             paste("MSS",1:length(JSS_list),sep="_"),
                             paste("templateSS",1:length(JSS_list),sep="_"),
                             paste("alphaS",1:length(JSS_list),sep="_"),
                             paste("betaS",1:length(JSS_list),sep="_")))
        } else{
          in_data       <- unique(c(in_data,"Nd","Nu","Jcause","alpha",
                                  paste("JSS",1:length(JSS_list),sep="_")[-single_column_MSS],
                                  paste("MSS",1:length(JSS_list),sep="_"),
                                  paste("templateSS",1:length(JSS_list),sep="_"),
                                  paste("alphaS",1:length(JSS_list),sep="_"),
                                  paste("betaS",1:length(JSS_list),sep="_")))
        }
        
        out_parameter <- unique(c(out_parameter, paste("thetaSS", seq_along(JSS_list), sep = "_"),
                           "pEti"))
    }
    
    # set up initializing function:
    if ("BrS" %in% use_measurements & !("SS" %in% use_measurements)){
        in_init       <-   function(){
          tmp_thetaBS <- list()
          tmp_psiBS <- list()
          for(i in seq_along(JBrS_list)){
            tmp_thetaBS[[i]] <- rbeta(JBrS_list[[i]],1,1)
          }
          for(i in seq_along(JBrS_list)){
            tmp_psiBS[[i]] <- rbeta(JBrS_list[[i]],1,1)       
          }
          res <- c(tmp_thetaBS,tmp_psiBS)
          names(res) <- c(paste("thetaBS", seq_along(JBrS_list), sep = "_"),
                          paste("psiBS", seq_along(JBrS_list), sep = "_"))
          res
        }
    }
    
    if (!("BrS" %in% use_measurements) & "SS" %in% use_measurements){
      in_init       <-   function(){
        tmp_thetaSS <- list()
        tmp_psiSS <- list()
        for(i in seq_along(JSS_list)){
          tmp_thetaSS[[i]] <- rbeta(JSS_list[[i]],1,1)
        }
        res <- c(tmp_thetaSS,tmp_psiSS)
        names(res) <- c(paste("thetaSS", seq_along(JSS_list), sep = "_"))
        res
      }
    } 
    if ("BrS" %in% use_measurements & "SS" %in% use_measurements){    
      in_init       <-   function(){
        tmp_thetaBS <- list()
        tmp_psiBS <- list()
        for(i in seq_along(JBrS_list)){
          tmp_thetaBS[[i]] <- rbeta(JBrS_list[[i]],1,1)
        }
        for(i in seq_along(JBrS_list)){
          tmp_psiBS[[i]] <- rbeta(JBrS_list[[i]],1,1)       
        }
        res <- c(tmp_thetaBS,tmp_psiBS)
        names(res) <- c(paste("thetaBS", seq_along(JBrS_list), sep = "_"),
                        paste("psiBS", seq_along(JBrS_list), sep = "_"))
        
        tmp_thetaSS <- list()
        tmp_psiSS <- list()
        for(i in seq_along(JSS_list)){
          tmp_thetaSS[[i]] <- rbeta(JSS_list[[i]],1,1)
        }
        res2 <- c(tmp_thetaSS,tmp_psiSS)
        names(res2) <- c(paste("thetaSS", seq_along(JSS_list), sep = "_"))
        
        c(res,res2)
      }
    }

    
    # etiology (measurement independent)
    alpha          <- rep(1,Jcause)
    
    #
    # fit model :
    #
    
    mybugs <- function(...){
      inits      <- in_init;
      data       <- in_data;
      parameters <- out_parameter;
      rst.bugs   <- call.bugs(data, inits, parameters,...);
      rst.bugs
    }
    
    #
    # write the .bug files into mcmc_options$bugsmodel.dir; 
    # could later set it equal to result.folder.
    #
    
    if (mcmc_options$ppd==TRUE){
      stop("not done.")
    } else {
      model_func         <- write_model_NoReg_plcm(data_nplcm$Mobs,model_options$likelihood$cause_list,model_options$use_measurements)
      model_bugfile_name <- "model_NoReg_plcm.bug"
    }
    
    filename <- file.path(mcmc_options$bugsmodel.dir, model_bugfile_name)
    writeLines(model_func, filename)


    #
    # run the model:
    #
    gs <- mybugs(model_bugfile_name)
}
