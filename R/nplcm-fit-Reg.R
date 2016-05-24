if(getRversion() >= "2.15.1") utils::globalVariables(c("set_prior_tpr","set_prior_eti"))

#' Fit nested partially-latent class model with regression (low-level)
#'
#' @details This function prepares data, specifies hyperparameters in priors 
#' (true positive rates and etiology fractions), initializes the posterior
#' sampling chain, writes the model file (for JAGS or WinBUGS with slight
#' differences in syntax), and fits the model. Features:
#' \itemize{
#' \item regression;
#' \item no nested subclasses, i.e. conditional indepdnence of 
#' multivariate measurements given disease class and covariates;
#' \item multiple BrS + multiple SS.
#' }
#' If running JAGS on windows, please go to control panel to add the directory to
#' jags into ENVIRONMENTAL VARIABLE!
#'
#' @inheritParams nplcm
#' @return BUGS fit results.
#' 
#' @seealso \link{write_model_NoReg} for constructing .bug model file; This function
#' then put it in the folder \code{mcmc_options$bugsmodel.dir}.
#' 
#' @family model fitting functions 
#' 
#' @export
nplcm_fit_Reg_NoNest <- 
  function(data_nplcm,model_options,mcmc_options){
    # Record the settings of current analysis:
    cat("==[baker] Results stored in: ==","\n",mcmc_options$result.folder,"\n")
    #model_options:
    dput(model_options,file.path(mcmc_options$result.folder,"model_options.txt"))
    #mcmc_options:
    dput(mcmc_options,file.path(mcmc_options$result.folder,"mcmc_options.txt"))
    
    # read in data:
    Mobs <- data_nplcm$Mobs
    Y    <- data_nplcm$Y
    X    <- data_nplcm$X
    # read in options:
    likelihood       <- model_options$likelihood
    use_measurements <- model_options$use_measurements
    prior            <- model_options$prior
    
    #####################################################################
    # 1. prepare data (including hyper-parameters):
    #####################################################################
    # get sample size:
    Nd      <- sum(Y==1)
    Nu      <- sum(Y==0)
    
    # get lengths of vectors:
    cause_list  <- likelihood$cause_list
    Jcause      <- length(cause_list)
    
    in_data <- in_init <- out_parameter <- NULL
    
    # generate design matrix for etiology regression:
    Eti_formula <- likelihood$Eti_formula
    FPR_formula <- likelihood$FPR_formula
    
    # design matrix for etiology regression:
    Z_Eti       <- stats::model.matrix(Eti_formula,data.frame(X,Y)[Y==1,,drop=FALSE])
    ncol_dm_Eti <- ncol(Z_Eti)
    attributes(Z_Eti)[names(attributes(Z_Eti))!="dim"] <- NULL
    
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
          warning(paste0("==[baker] Bronze-standard slice ", names(data_nplcm$Mobs$MBS)[s], " has no measurements informative of the causes! Please check if measurements' columns correspond to causes.=="))  
        }
      }
      
      MBS.case_list <- lapply(Mobs$MBS,"[",which(Y==1),TRUE,drop=FALSE)
      MBS.ctrl_list <- lapply(Mobs$MBS,"[",which(Y==0),TRUE,drop=FALSE)
      MBS_list <- list()
      for (i in seq_along(MBS.case_list)){
        MBS_list[[i]]      <- rbind(MBS.case_list[[i]],MBS.ctrl_list[[i]])
      }
      names(MBS_list)   <- names(MBS.case_list)
      
      single_column_MBS <- which(lapply(MBS_list,ncol)==1)
      
      
      # input design matrix for FPR regressions:
      int_Y <- as.integer(Y)
      if (any(int_Y[1:sum(int_Y)]==0) | any(int_Y[-(1:sum(int_Y))])==1){
        stop("==[baker] Please put cases at the top of controls in `data_nplcm`.==\n")
      }
      
      Z_FPR_list <- lapply(FPR_formula,function(form){stats::model.matrix(form,data.frame(X,Y))}) # <-- make sure that the row orders are the same.
      
      #intercept_only_MBS <- which(lapply(Z_FPR_list,function(x) (ncol(x)==1 & all(x==1)))==TRUE)
      
      has_basis_list <- lapply(Z_FPR_list, function(Z) {
        if (length(grep("^s_",dimnames(Z)[[2]]))==0){
          FALSE
        } else{
          TRUE
        }
      })
      basis_id_list <- lapply(Z_FPR_list, function(Z) {
        if (length(grep("^s_",dimnames(Z)[[2]]))==0){
          NULL
        } else{
          0+grep("^s_",dimnames(Z)[[2]])
        }
      })
      n_basis_list  <- lapply(basis_id_list,length) 
      non_basis_id_list <- lapply(Z_FPR_list,function(Z){
        if (length(grep("^s_",dimnames(Z)[[2]]))==0){
          0+(1:ncol(Z))
        } else{
          (1:ncol(Z))[-grep("^s_",dimnames(Z)[[2]])]
        }
      })
      
      for(i in seq_along(JBrS_list)){
        attributes(Z_FPR_list[[i]])[names(attributes(Z_FPR_list[[i]]))!="dim"] <- NULL
        
        assign(paste("JBrS", i, sep = "_"), JBrS_list[[i]])    
        assign(paste("MBS", i, sep = "_"), as.matrix_or_vec(MBS_list[[i]])) 
        assign(paste("templateBS", i, sep = "_"), as.matrix_or_vec(template_BrS_list[[i]]))   
        assign(paste("Z_FPR",i,sep="_"),Z_FPR_list[[i]])
        assign(paste("basis_id",i,sep="_"),basis_id_list[[i]])
        assign(paste("n_basis",i,sep="_"),n_basis_list[[i]])
        assign(paste("non_basis_id",i,sep="_"),non_basis_id_list[[i]])
      }
      
      # summarize into one name (for all measurements):
      if (length(single_column_MBS)==0){
        # if all slices have >2 columns:
        in_data       <- c(in_data,"Nd","Nu","Jcause","Z_Eti","ncol_dm_Eti",
                           paste("JBrS",1:length(JBrS_list),sep="_"),
                           paste("MBS",1:length(JBrS_list),sep="_"),
                           paste("templateBS",1:length(JBrS_list),sep="_"),
                           paste("Z_FPR",1:length(JBrS_list),sep="_"),
                           paste("basis_id",1:length(JBrS_list),sep="_")[unlist(has_basis_list)],
                           paste("n_basis",1:length(JBrS_list),sep="_")[unlist(has_basis_list)],
                           paste("non_basis_id",1:length(JBrS_list),sep="_")
                           # paste("alphaB",1:length(JBrS_list),sep="_"),
                           # paste("betaB",1:length(JBrS_list),sep="_")
        )
      } else {
        # if there exist slices with 1 column:
        in_data       <- c(in_data,"Nd","Nu","Jcause","Z_Eti","ncol_dm_Eti",
                           paste("JBrS",1:length(JBrS_list),sep="_")[-single_column_MBS], # <---- no need to iterate in .bug file for a slice with one column.
                           paste("MBS",1:length(JBrS_list),sep="_"),
                           paste("templateBS",1:length(JBrS_list),sep="_"),
                           paste("Z_FPR",1:length(JBrS_list),sep="_"),
                           paste("basis_id",1:length(JBrS_list),sep="_")[unlist(has_basis_list)],
                           paste("n_basis",1:length(JBrS_list),sep="_")[unlist(has_basis_list)],
                           paste("non_basis_id",1:length(JBrS_list),sep="_")
                           # paste("alphaB",1:length(JBrS_list),sep="_"),
                           # paste("betaB",1:length(JBrS_list),sep="_")
        )
      }
      
      #
      # hyper-parameters:
      #
      
      # Set BrS measurement priors:
      # hyperparameter for sensitivity (can add for specificity if necessary): 
      for (s in seq_along(Mobs$MBS)){
        if (likelihood$k_subclass[s] == 1){BrS_tpr_prior <- set_prior_tpr_BrS_NoNest(s,model_options,data_nplcm)}
        if (likelihood$k_subclass[s] > 1){BrS_tpr_prior <- set_prior_tpr_BrS_NoNest(s,model_options,data_nplcm)}
        
        assign(paste("alphaB", s, sep = "_"), BrS_tpr_prior[[1]]$alpha)     # <---- input BrS TPR prior here.
        assign(paste("betaB", s, sep = "_"),  BrS_tpr_prior[[1]]$beta)    
        
        if (likelihood$k_subclass[s]==1){
          out_parameter <- c(out_parameter,paste(c("thetaBS","betaFPR"), s, sep="_"))
          out_parameter <- c(out_parameter,paste(c("taubeta"), s, sep="_")[unlist(has_basis_list)[s]])
        }else{
          stop("==[baker] Regression with nested subclasses coming soon.==\n")
          #assign(paste("K", s, sep = "_"), likelihood$k_subclass[s])
          #in_data       <- c(in_data,paste0("K_",s)) # <---- not prior, but data about the subclasses for this slice.
          #out_parameter <- c(out_parameter,
          #                   paste(c("ThetaBS","PsiBS","Lambda","Eta","alphadp0","alphadp0_case"),s,sep="_")
          #)
        }
      }
      
      #
      # collect in_data, out_parameter together:
      #
      in_data       <- c(in_data,
                         paste("alphaB",1:length(JBrS_list),sep="_"),
                         paste("betaB",1:length(JBrS_list),sep="_")
      )
      out_parameter <- c(out_parameter,"betaEti")
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
          warning(paste0("==[baker] Silver-standard slice ", names(data_nplcm$Mobs$MSS)[s], 
                         " has no measurements informative of the causes! Please check if measurements' columns correspond to causes.==\n"))  
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
      prior_SS     <- model_options$prior$TPR_prior$SS
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
          assign(paste("GSS_TPR", i, sep = "_"), length(unique(prior_SS$grp))) # <--- need to change to depending on i if grp change wrt specimen.
        }
        if (is.null(prior_SS$grp)){ # <--- need to change to list if we have multiple slices.
          assign(paste("GSS_TPR", i, sep = "_"), 1)
        }
      }
      
      SS_tpr_prior <- set_prior_tpr_SS(model_options,data_nplcm)
      
      # set SS measurement priors: 
      # hyper-parameters for sensitivity:
      alpha_mat <- list() # dimension for slices.
      beta_mat  <- list()
      
      for(i in seq_along(JSS_list)){
        
        GSS_TPR_curr <- eval(parse(text = paste0("GSS_TPR_",i)))
        alpha_mat[[i]] <- matrix(NA, nrow=GSS_TPR_curr,ncol=JSS_list[[i]])
        beta_mat[[i]]  <- matrix(NA, nrow=GSS_TPR_curr,ncol=JSS_list[[i]])
        
        colnames(alpha_mat[[i]]) <- patho_SS_list[[i]]
        colnames(beta_mat[[i]]) <- patho_SS_list[[i]]
        
        for (g in 1:GSS_TPR_curr){
          alpha_mat[[i]][g,] <- unlist(SS_tpr_prior[[i]][[g]]$alpha)
          beta_mat[[i]][g,]  <- unlist(SS_tpr_prior[[i]][[g]]$beta)
        }
        
        if (GSS_TPR_curr>1){
          assign(paste("alphaS", i, sep = "_"), alpha_mat[[i]])    # <---- input SS TPR prior here.
          assign(paste("betaS", i, sep = "_"),  beta_mat[[i]])    
        }else{
          assign(paste("alphaS", i, sep = "_"), c(alpha_mat[[i]]))   # <---- input SS TPR prior here.
          assign(paste("betaS", i, sep = "_"),  c(beta_mat[[i]]))    
        }
      }
      names(alpha_mat) <- names(beta_mat)<- names(Mobs$MSS)
      
      if (!SS_TPR_strat){
        if (length(single_column_MSS)==0){
          # summarize into one name (for all measurements):
          in_data       <- unique(c(in_data,"Nd","Jcause","ncol_dm_Eti",
                                    paste("JSS",1:length(JSS_list),sep="_"),
                                    paste("MSS",1:length(JSS_list),sep="_"),
                                    paste("templateSS",1:length(JSS_list),sep="_"),
                                    paste("alphaS",1:length(JSS_list),sep="_"),   
                                    paste("betaS",1:length(JSS_list),sep="_")))
        } else{
          in_data       <- unique(c(in_data,"Nd","Nu","Jcause","ncol_dm_Eti",
                                    paste("JSS",1:length(JSS_list),sep="_")[-single_column_MSS],
                                    paste("MSS",1:length(JSS_list),sep="_"),
                                    paste("templateSS",1:length(JSS_list),sep="_"),
                                    paste("alphaS",1:length(JSS_list),sep="_"),
                                    paste("betaS",1:length(JSS_list),sep="_")))
        }
      }else {
        if (length(single_column_MSS)==0){
          # summarize into one name (for all measurements):
          in_data       <- unique(c(in_data,"Nd","Jcause","ncol_dm_Eti",
                                    paste("JSS",1:length(JSS_list),sep="_"),
                                    paste("GSS_TPR",1:length(JSS_list),sep="_"),
                                    paste("SS_TPR_grp",1:length(JSS_list),sep="_"),
                                    paste("MSS",1:length(JSS_list),sep="_"),
                                    paste("templateSS",1:length(JSS_list),sep="_"),
                                    paste("alphaS",1:length(JSS_list),sep="_"),   
                                    paste("betaS",1:length(JSS_list),sep="_")))
        } else{
          in_data       <- unique(c(in_data,"Nd","Nu","Jcause","ncol_dm_Eti",
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
                                "betaEti"))
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
            res_curr[[1]] <- stats::rbeta(JBrS_list[[s]],1,1)
            res_curr[[2]] <- matrix(0,nrow=ncol(Z_FPR_list[[s]]),ncol=JBrS_list[[s]])
            names(res_curr) <- paste(c("thetaBS","betaFPR"),s,sep="_")
            
            res <- c(res,res_curr)
          }
          if (likelihood$k_subclass[s] > 1){
            stop("==[baker] Nested subclasses for regression coming soon.==\n")
            # K_curr <- likelihood$k_subclass[s]
            # res_curr[[1]] <- c(rep(.5,K_curr-1),NA)
            # res_curr[[2]] <- c(rep(.5,K_curr-1),NA)
            # # for different Eta's:
            # #             res_curr[[2]] <- cbind(matrix(rep(.5,Jcause*(K_curr-1)),
            # #                                           nrow=Jcause,ncol=K_curr-1),
            # #                                    rep(NA,Jcause))
            # res_curr[[3]] <- 1
            # res_curr[[4]] <- 1 #<---- added together with 'alphadp0_case' below.
            # 
            # names(res_curr) <- paste(c("r0","r1","alphadp0","alphadp0_case"),s,sep="_")
            # res <- c(res,res_curr)
          }
        }
        
        res[[length(res)+1]] <- cbind(matrix(0,nrow=ncol(Z_Eti),ncol=Jcause-1),rep(NA,ncol(Z_Eti)))
        names(res)[length(res)] <- "betaEti"
        res
      }
    }
    
    if (!("BrS" %in% use_measurements) & "SS" %in% use_measurements){
      in_init       <-   function(){
        tmp_thetaSS <- list()
        #tmp_psiSS <- list()
        tmp_Icat_case <- list()
        for(i in seq_along(JSS_list)){
          GSS_TPR_curr <- eval(parse(text = paste0("GSS_TPR_",i)))
          if (GSS_TPR_curr==1){
            tmp_thetaSS[[i]] <- stats::rbeta(JSS_list[[i]],1,1)
          } else{
            tmp_thetaSS[[i]] <- matrix(stats::rbeta(GSS_TPR_curr*JSS_list[[i]],1,1),
                                       nrow=GSS_TPR_curr,ncol=JSS_list[[i]])
          }
          if (i==1){
            #if (length(JSS_list)>1){
            #  warning("==[baker] Only the first slice of silver-standard data is used to 
            #          initialize 'Icat' in JAGS fitting. Choose wisely!\n ==")
            #}  
            tmp_Icat_case[[i]] <- init_latent_jags_multipleSS(MSS_list,likelihood$cause_list)
          }
        }
        
        res <- c(tmp_thetaSS,tmp_Icat_case)
        names(res) <- c(paste("thetaSS", seq_along(JSS_list), sep = "_"),"Icat")
        res
      }
    } 
    
    if ("BrS" %in% use_measurements & "SS" %in% use_measurements){    
      in_init       <-   function(){
        res <- list()
        for (s in seq_along(Mobs$MBS)){
          res_curr <- list()
          if (likelihood$k_subclass[s]==1){
            res_curr[[1]] <- stats::rbeta(JBrS_list[[s]],1,1)
            res_curr[[2]] <- matrix(0,nrow=ncol(Z_FPR_list[[s]]),ncol=JBrS_list[[s]])
            names(res_curr) <- paste(c("thetaBS","betaFPR"),s,sep="_")
            res <- c(res,res_curr)
          }
          if (likelihood$k_subclass[s] > 1){
            stop("==[baker] Nested subclasses for regression coming soon.==\n")
            # K_curr <- likelihood$k_subclass[s]
            # res_curr[[1]] <- c(rep(.5,K_curr-1),NA)
            # res_curr[[2]] <- c(rep(.5,K_curr-1),NA)
            # # for different Eta's:
            # # res_curr[[2]] <- cbind(matrix(rep(.5,Jcause*(K_curr-1)),
            # #                               nrow=Jcause,ncol=K_curr-1),
            # #                              rep(NA,Jcause))
            # res_curr[[3]] <- 1
            # res_curr[[4]] <- 1
            # 
            # names(res_curr) <- paste(c("r0","r1","alphadp0","alphadp0_case"),s,sep="_")
            # res <- c(res,res_curr)
          }
        }
        
        tmp_thetaSS <- list()
        #tmp_psiSS <- list()
        tmp_Icat_case <- list()
        for(i in seq_along(JSS_list)){
          GSS_TPR_curr <- eval(parse(text = paste0("GSS_TPR_",i)))
          if (GSS_TPR_curr==1){
            tmp_thetaSS[[i]] <- stats::rbeta(JSS_list[[i]],1,1)
          } else{
            tmp_thetaSS[[i]] <- matrix(stats::rbeta(GSS_TPR_curr*JSS_list[[i]],1,1),nrow=GSS_TPR_curr,ncol=JSS_list[[i]])
          }
          if (i==1){
            #if (length(JSS_list)>1){
            # warning("==[baker] Only the first slice of silver-standard data is
            #         used to initialize 'Icat' in JAGS fitting. Choose wisely!==\n ")
            #}
            tmp_Icat_case[[i]] <- init_latent_jags_multipleSS(MSS_list,likelihood$cause_list)
          }
        }
        res2 <- c(tmp_thetaSS,tmp_Icat_case)
        names(res2) <- c(paste("thetaSS", seq_along(JSS_list), sep = "_"),"Icat")
        #print(res2)
        res[[length(res)+1]] <- cbind(matrix(0,nrow=ncol(Z_Eti),ncol=Jcause-1),rep(NA,ncol(Z_Eti)))
        names(res)[length(res)] <- "betaEti"
        
        c(res,res2)
      }
    }
    
    #
    # fit model :
    #
    
    # special operations:
    # get individual prediction:
    if (!is.null(mcmc_options$individual.pred) && mcmc_options$individual.pred){out_parameter <- c(out_parameter,"Icat")}
    # get posterior predictive distribution of BrS measurments:
    if (!is.null(mcmc_options$ppd) && mcmc_options$ppd){out_parameter <- c(out_parameter,paste("MBS.new",seq_along(Mobs$MBS),sep = "_"))}
    # get pEti: (for regression models; besides the coefficients of etiology regressions are always recorded)
    if (!is.null(mcmc_options$get.pEti) && mcmc_options$get.pEti){out_parameter <- c(out_parameter,"pEti")}
    
    #
    # write the .bug files into mcmc_options$bugsmodel.dir; 
    # could later set it equal to result.folder.
    #
    use_jags <- (!is.null(mcmc_options$use_jags) && mcmc_options$use_jags)
    model_func         <- write_model_Reg_NoNest(data_nplcm$Mobs,
                                                 model_options$prior,
                                                 model_options$likelihood$cause_list,
                                                 model_options$likelihood$Eti_formula,
                                                 model_options$likelihood$FPR_formula,
                                                 model_options$use_measurements,
                                                 mcmc_options$ppd,
                                                 use_jags)
    model_bugfile_name <- "model_Reg_NoNest.bug"
    
    filename <- file.path(mcmc_options$bugsmodel.dir, model_bugfile_name)
    writeLines(model_func, filename)
    
    in_data <- unique(in_data)
    
    if (Jcause > 2){ # use dmnorm to speed up JAGS calculations, so we need mean 
      # and covariance vector.
      zero_Jcause_1 <- rep(0,Jcause-1)
      I_Jcause_1    <- diag(1,Jcause-1)
      in_data <- c(in_data,"zero_Jcause_1","I_Jcause_1")
    } 
    
    for (s in seq_along(JBrS_list)){
      if (JBrS_list[[s]]>1){
          assign(paste("zero_JBrS", s, sep = "_"), rep(0,JBrS_list[[s]]))    
          assign(paste("I_JBrS", s, sep = "_"), diag(1,JBrS_list[[s]]))    
          in_data <- c(in_data,
                       paste("zero_JBrS", s, sep = "_"),
                       paste("I_JBrS", s, sep = "_"))
      }
    }
    
    #
    # run the model:
    #
    if (!use_jags){
      stop("==[baker] Regression model in WinBUGS coming soon. Please use JAGS for now.==\n")
      # ##winbugs is the only current option:
      # gs <- R2WinBUGS::bugs(data     = in_data,
      #                       inits    = in_init, 
      #                       parameters.to.save = out_parameter,
      #                       model.file = filename,
      #                       working.directory=mcmc_options$result.folder,
      #                       bugs.directory  = mcmc_options$winbugs.dir,  #<- special to WinBUGS.
      #                       n.iter         = mcmc_options$n.itermcmc,
      #                       n.burnin       = mcmc_options$n.burnin,
      #                       n.thin         = mcmc_options$n.thin,
      #                       n.chains       = mcmc_options$n.chains,
      #                       DIC      = FALSE,
      #                       debug    = mcmc_options$debugstatus);
      # return(gs)
    }
    
    here <- environment()
    
    if (use_jags){
      ##JAGS
      in_data.list <- lapply(as.list(in_data), get, envir= here)
      names(in_data.list) <- in_data
      #lapply(names(in_data.list), dump, append = TRUE, envir = here,
      #       file = file.path(mcmc_options$result.folder,"jagsdata.txt"))
      dump(names(in_data.list), append = FALSE, envir = here,
             file = file.path(mcmc_options$result.folder,"jagsdata.txt"))
      gs <- R2jags::jags2(data   = in_data,
                          inits  = in_init,
                          parameters.to.save = out_parameter,
                          model.file = filename,
                          working.directory = mcmc_options$result.folder,
                          n.iter         = as.integer(mcmc_options$n.itermcmc),
                          n.burnin       = as.integer(mcmc_options$n.burnin),
                          n.thin         = as.integer(mcmc_options$n.thin),
                          n.chains       = as.integer(mcmc_options$n.chains),
                          DIC            = FALSE,
                          clearWD        = FALSE,              #<--- special to JAGS.
                          jags.path      = mcmc_options$jags.dir# <- special to JAGS.
      );
      return(gs)
    }
    
  }


#nplcm_fit_Reg_NoNest(data_nplcm,model_options,mcmc_options)
