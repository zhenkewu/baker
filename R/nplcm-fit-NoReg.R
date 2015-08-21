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
#' @seealso \link{write_model_NoReg} for writing .bug model file; 
#' 
#' @export
#' 
#' 
nplcm_fit_NoReg<-
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
    prior            <- model_options$prior
    
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
    
    
    #####################################################################
    # 1. prepare data (including hyper-parameters):
    #####################################################################
    
    # get sample sizes:
    Nd      <- sum(Y==1)
    Nu      <- sum(Y==0)
    
    # get lengths of vectors:
    cause_list  <- likelihood$cause_list
    Jcause      <- length(cause_list)
    
    in_data <- in_init <- out_parameter <- NULL
    

    if ("BrS" %in% use_measurements){
      
      #
      # BrS measurement data: 
      #
      JBrS_list  <- lapply(Mobs$MBS,ncol)
      # mapping template (by `make_template` function):
      patho_BrS_list    <- lapply(Mobs$MBS,colnames)
      template_BrS_list <- lapply(patho_BrS_list,make_template,cause_list)
      for (s in seq_along(template_BrS_list)){
        if (sum(template_BrS_list[[s]])==0){
          warning(paste0("== Bronze-standard slice ", names(data_nplcm$Mobs$MBS)[s], " has no measurements informative of the causes! Please check if measurements' columns correspond to causes.=="))  
        }
      }
      
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
      
      if (length(single_column_MBS)==0){
        # summarize into one name (for all measurements):
        in_data       <- c(in_data,"Nd","Nu","Jcause","alpha",
                           paste("JBrS",1:length(JBrS_list),sep="_"),
                           paste("MBS",1:length(JBrS_list),sep="_"),
                           paste("templateBS",1:length(JBrS_list),sep="_")
                          # paste("alphaB",1:length(JBrS_list),sep="_"),
                          # paste("betaB",1:length(JBrS_list),sep="_")
        )
      } else {
        # summarize into one name (for all measurements):
        in_data       <- c(in_data,"Nd","Nu","Jcause","alpha",
                           paste("JBrS",1:length(JBrS_list),sep="_")[-single_column_MBS],
                           paste("MBS",1:length(JBrS_list),sep="_"),
                           paste("templateBS",1:length(JBrS_list),sep="_")
                          # paste("alphaB",1:length(JBrS_list),sep="_"),
                          # paste("betaB",1:length(JBrS_list),sep="_")
        )
      }
      
      #
      # hyper-parameters:
      #
      
      # set BrS measurement priors:
      # hyperparameter for sensitivity (can add for specificity if necessary): 
      
      for (s in seq_along(Mobs$MBS)){
        if (likelihood$k_subclass[s]==1){BrS_tpr_prior <- set_prior_tpr_BrS_NoNest(s,model_options,data_nplcm)}
        if (likelihood$k_subclass[s] > 1){BrS_tpr_prior <- set_prior_tpr_BrS_NoNest(s,model_options,data_nplcm)}
        
        assign(paste("alphaB", s, sep = "_"), BrS_tpr_prior[[1]]$alpha)     # <---- input BrS TPR prior here.
        assign(paste("betaB", s, sep = "_"),  BrS_tpr_prior[[1]]$beta)    
        
        if (likelihood$k_subclass[s]==1){
         
         out_parameter <- c(out_parameter,paste(c("thetaBS","psiBS"), s, sep="_"))  
        }else{
         assign(paste("K", s, sep = "_"), likelihood$k_subclass[s])
         in_data       <- c(in_data,paste0("K_",s)) # <---- not prior, but data about the subclasses for this slice.
         out_parameter <- c(out_parameter,
                            paste(c("ThetaBS","PsiBS","Lambda","Eta","alphadp0"),s,sep="_")
                            )
        }
      }
      
      #
      # collect in_data, out_parameter together:
      #
      in_data       <- c(in_data,
                         paste("alphaB",1:length(JBrS_list),sep="_"),
                         paste("betaB",1:length(JBrS_list),sep="_")
      )
      out_parameter <- c(out_parameter,"pEti")
    }
    
  
    if ("SS" %in% use_measurements){
      #
      # 2. SS measurement data: 
      #
      
      JSS_list   <- lapply(Mobs$MSS,ncol)
      
      # mapping template (by `make_template` function):
      patho_SS_list <- lapply(Mobs$MSS,colnames)
      template_SS_list <- lapply(patho_SS_list,make_template,cause_list)
      
      for (s in seq_along(template_SS_list)){
        if (sum(template_SS_list[[s]])==0){
          warning(paste0("== Silver-standard slice ", names(data_nplcm$Mobs$MSS)[s], " has no measurements informative of the causes! Please check if measurements' columns correspond to causes.=="))  
        }
      }
      
      MSS_list <- lapply(Mobs$MSS,"[",which(Y==1),TRUE,drop=FALSE)
      
      single_column_MSS <- which(lapply(MSS_list,ncol)==1)
      
      for(i in seq_along(JSS_list)){
        assign(paste("JSS", i, sep = "_"), JSS_list[[i]])    
        assign(paste("MSS", i, sep = "_"), as.matrix_or_vec(MSS_list[[i]])) 
        assign(paste("templateSS", i, sep = "_"), as.matrix_or_vec(template_SS_list[[i]]))   
      }
      
      # setup groupwise TPR for SS:
      SS_TPR_strat <- FALSE
      prior_SS <- model_options$prior$TPR_prior$SS
      parsed_model <- assign_model(model_options,data_nplcm)
      if (parsed_model$SS_grp){
        SS_TPR_strat <- TRUE
        for(i in seq_along(JSS_list)){
          assign(paste("GSS_TPR", i, sep = "_"),  length(unique(prior_SS$grp)))    
          assign(paste("SS_TPR_grp", i, sep = "_"),  prior_SS$grp)    
        }
      }
      
      # add GSS_TPR_1, or 2 if we want to index by slices:
      for (i in seq_along(JSS_list)){
        if (!is.null(prior_SS$grp)){ # <--- need to change to list if we have multiple slices.
          assign(paste("GSS_TPR", i, sep = "_"), length(unique(prior_SS$grp)))
        }
        if (is.null(prior_SS$grp)){ # <--- need to change to list if we have multiple slices.
          assign(paste("GSS_TPR", i, sep = "_"), 1)
        }
      }
      
      
      SS_tpr_prior <- set_prior_tpr_SS(model_options,data_nplcm)
      
      # set SS measurement priors: 
      # hyper-parameters for sensitivity:
      for(i in seq_along(JSS_list)){
        
        GSS_TPR_curr <- eval(parse(text = paste0("GSS_TPR_",i)))
        alpha_mat <- list() # dimension for slices.
        beta_mat  <- list()
        alpha_mat[[i]] <- matrix(NA, nrow=GSS_TPR_curr,ncol=JSS_list[[i]])
        beta_mat[[i]]  <- matrix(NA, nrow=GSS_TPR_curr,ncol=JSS_list[[i]])
        
        colnames(alpha_mat[[i]]) <- patho_SS_list[[i]]
        colnames(beta_mat[[i]]) <- patho_SS_list[[i]]
        
        for (g in 1:GSS_TPR_curr){
          alpha_mat[[i]][g,] <- unlist(SS_tpr_prior[[i]][[g]]$alpha)
          beta_mat[[i]][g,]  <- unlist(SS_tpr_prior[[i]][[g]]$beta)
        }
        names(alpha_mat) <- names(beta_mat) <- names(Mobs$MSS)
        
        if (GSS_TPR_curr>1){
          assign(paste("alphaS", i, sep = "_"), alpha_mat[[i]])    # <---- input SS TPR prior here.
          assign(paste("betaS", i, sep = "_"),  beta_mat[[i]])    
        }else{
          assign(paste("alphaS", i, sep = "_"), c(alpha_mat[[i]]))   # <---- input SS TPR prior here.
          assign(paste("betaS", i, sep = "_"),  c(beta_mat[[i]]))    
        }
      }
      
      if (!SS_TPR_strat){
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
      }else {
        if (length(single_column_MSS)==0){
          # summarize into one name (for all measurements):
          in_data       <- unique(c(in_data,"Nd","Jcause","alpha",
                                    paste("JSS",1:length(JSS_list),sep="_"),
                                    paste("GSS_TPR",1:length(JSS_list),sep="_"),
                                    paste("SS_TPR_grp",1:length(JSS_list),sep="_"),
                                    paste("MSS",1:length(JSS_list),sep="_"),
                                    paste("templateSS",1:length(JSS_list),sep="_"),
                                    paste("alphaS",1:length(JSS_list),sep="_"),   
                                    paste("betaS",1:length(JSS_list),sep="_")))
        } else{
          in_data       <- unique(c(in_data,"Nd","Nu","Jcause","alpha",
                                    paste("JSS",1:length(JSS_list),sep="_")[-single_column_MSS],
                                    paste("GSS_TPR",1:length(JSS_list),sep="_"),
                                    paste("SS_TPR_grp",1:length(JSS_list),sep="_"),
                                    paste("MSS",1:length(JSS_list),sep="_"),
                                    paste("templateSS",1:length(JSS_list),sep="_"),
                                    paste("alphaS",1:length(JSS_list),sep="_"),
                                    paste("betaS",1:length(JSS_list),sep="_")))
        }
        
      }
      out_parameter <- unique(c(out_parameter, paste("thetaSS", seq_along(JSS_list), sep = "_"),
                                "pEti"))
    }
    
    #####################################
    # set up initialization function:
    #####################################
    if ("BrS" %in% use_measurements & !("SS" %in% use_measurements)){
      in_init       <-   function(){
        res <- list()
        for (s in seq_along(Mobs$MBS)){
          res_curr <- list()
          if (likelihood$k_subclass[s]==1){
            res_curr[[1]] <- rbeta(JBrS_list[[s]],1,1)
            res_curr[[2]] <- rbeta(JBrS_list[[s]],1,1)
            names(res_curr) <- paste(c("thetaBS","psiBS"),s,sep="_")
            res <- c(res,res_curr)
          }
          if (likelihood$k_subclass[s] > 1){
            K_curr <- likelihood$k_subclass[s]
            res_curr[[1]] <- c(rep(.5,K_curr-1),NA)
            res_curr[[2]] <- c(rep(.5,K_curr-1),NA)
            # for different Eta's:
            #             res_curr[[2]] <- cbind(matrix(rep(.5,Jcause*(K_curr-1)),
            #                                           nrow=Jcause,ncol=K_curr-1),
            #                                    rep(NA,Jcause))
            res_curr[[3]] <- 1
            
            names(res_curr) <- paste(c("r0","r1","alphadp0"),s,sep="_")
            res <- c(res,res_curr)
          }
        }
        res
      }
    }
    
    if (!("BrS" %in% use_measurements) & "SS" %in% use_measurements){
      in_init       <-   function(){
        tmp_thetaSS <- list()
        tmp_psiSS <- list()
        for(i in seq_along(JSS_list)){
          GSS_TPR_curr <- eval(parse(text = paste0("GSS_TPR_",i)))
          if (GSS_TPR_curr==1){
            tmp_thetaSS[[i]] <- rbeta(JSS_list[[i]],1,1)
          } else{
            tmp_thetaSS[[i]] <- matrix(rbeta(GSS_TPR_curr*JSS_list[[i]],1,1),nrow=GSS_TPR_curr,ncol=JSS_list[[i]])
          }
        }
        res <- c(tmp_thetaSS,tmp_psiSS)
        names(res) <- c(paste("thetaSS", seq_along(JSS_list), sep = "_"))
        res
      }
    } 
    
    if ("BrS" %in% use_measurements & "SS" %in% use_measurements){    
      in_init       <-   function(){
        res <- list()
        for (s in seq_along(Mobs$MBS)){
          res_curr <- list()
          if (likelihood$k_subclass[s]==1){
            res_curr[[1]] <- rbeta(JBrS_list[[s]],1,1)
            res_curr[[2]] <- rbeta(JBrS_list[[s]],1,1)
            names(res_curr) <- paste(c("thetaBS","psiBS"),s,sep="_")
            res <- c(res,res_curr)
          }
          if (likelihood$k_subclass[s] > 1){
            K_curr <- likelihood$k_subclass[s]
            res_curr[[1]] <- c(rep(.5,K_curr-1),NA)
            res_curr[[2]] <- c(rep(.5,K_curr-1),NA)
            # for different Eta's:
            # res_curr[[2]] <- cbind(matrix(rep(.5,Jcause*(K_curr-1)),
            #                               nrow=Jcause,ncol=K_curr-1),
            #                              rep(NA,Jcause))
            res_curr[[3]] <- 1
            
            names(res_curr) <- paste(c("r0","r1","alphadp0"),s,sep="_")
            res <- c(res,res_curr)
          }
        }
        
        tmp_thetaSS <- list()
        tmp_psiSS <- list()
        for(i in seq_along(JSS_list)){
          GSS_TPR_curr <- eval(parse(text = paste0("GSS_TPR_",i)))
          if (GSS_TPR_curr==1){
            tmp_thetaSS[[i]] <- rbeta(JSS_list[[i]],1,1)
          } else{
            tmp_thetaSS[[i]] <- matrix(rbeta(GSS_TPR_curr*JSS_list[[i]],1,1),nrow=GSS_TPR_curr,ncol=JSS_list[[i]])
          }
        }
        res2 <- c(tmp_thetaSS,tmp_psiSS)
        names(res2) <- c(paste("thetaSS", seq_along(JSS_list), sep = "_"))
        #print(res2)
        c(res,res2)
      }
    }
    
    # etiology (measurement independent)
    alpha          <- prior$Eti_prior    # <-------- input etiology prior here.
    
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
    }  
    if (mcmc_options$ppd==FALSE){
      model_func         <- write_model_NoReg(model_options$likelihood$k_subclass,
                                              data_nplcm$Mobs,
                                              model_options$prior,
                                              model_options$likelihood$cause_list,
                                              model_options$use_measurements)
      model_bugfile_name <- "model_NoReg.bug"
    }
    
    filename <- file.path(mcmc_options$bugsmodel.dir, model_bugfile_name)
    writeLines(model_func, filename)
    
    
    #
    # run the model:
    #
    gs <- mybugs(model_bugfile_name)
}
