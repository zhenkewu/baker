if(getRversion() >= "2.15.1") utils::globalVariables(c("set_prior_tpr","set_prior_eti"))


#' Fit nested partially-latent class model with regression (low-level)
#'
#' @details This function prepares data, specifies hyperparameters in priors 
#' (true positive rates and etiology fractions), initializes the posterior
#' sampling chain, writes the model file (for JAGS or WinBUGS with slight
#' differences in syntax), and fits the model. Features:
#' \itemize{
#' \item regression;
#' \item no nested subclasses, i.e. conditional independence of 
#' multivariate measurements given disease class and covariates;
#' \item multiple BrS + multiple SS.
#' }
#' If running JAGS on windows, please go to control panel to add the directory to
#' jags into ENVIRONMENTAL VARIABLE!
#'
#' @inheritParams nplcm
#' @return BUGS fit results.
#' 
#' @seealso \link{write_model_NoReg} for automatically generate \code{.bug} model
#'file; This present function store it in location: \code{mcmc_options$bugsmodel.dir}.
#' 
#' @family model fitting functions 
#' 
#' @export
nplcm_fit_Reg_discrete_predictor_NoNest <- 
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
    
    is_discrete_Eti     <- is_discrete(data.frame(X,Y)[Y==1,,drop=FALSE], Eti_formula)
    is_discrete_FPR_vec <- rep(NA,length(FPR_formula))
    names(is_discrete_FPR_vec) <- names(FPR_formula)
    for (s in seq_along(FPR_formula)){
      is_discrete_FPR_vec[s] <- is_discrete(X, FPR_formula[[s]])
    }
    
    # test discrete predictors (why do this when we already assigned the model?):
    if (!is_discrete_Eti){
      stop("==[baker]Etiology regression has non-discrete covariates! ==")
    }
    if (any(!is_discrete_FPR_vec)){
      stop("==[baker]FPR regression has non-discrete covariates! ==")
    }
    
    # design matrix for etiology regression (because this function is 
    # for all discrete predictors, the usual model.matrix is sufficient):
    Z_Eti       <- stats::model.matrix(Eti_formula,data.frame(X,Y)[Y==1,,drop=FALSE])
    # Z_Eti0       <- stats::model.matrix(Eti_formula,data.frame(X,Y)[Y==1,,drop=FALSE])
    # a.eig        <- eigen(t(Z_Eti0)%*%Z_Eti0)
    # sqrt_Z_Eti0  <- a.eig$vectors %*% diag(sqrt(a.eig$values)) %*% solve(a.eig$vectors)
    # 
    # Z_Eti  <- Z_Eti0%*%solve(sqrt_Z_Eti0)
    
    ncol_dm_Eti        <- ncol(Z_Eti)
    Eti_colname_design_mat <- attributes(Z_Eti)$dimnames[[2]]
    attributes(Z_Eti)[names(attributes(Z_Eti))!="dim"] <- NULL 
    # this to prevent issues when JAGS reads in data.
    
    unique_Eti_level   <- unique(Z_Eti)
    n_unique_Eti_level <- nrow(unique_Eti_level)
    Eti_stratum_id     <- apply(Z_Eti,1,function(v) 
      which(rowSums(abs(unique_Eti_level-t(replicate(n_unique_Eti_level,v))))==0))
    rownames(unique_Eti_level) <- 1:n_unique_Eti_level
    colnames(unique_Eti_level) <- Eti_colname_design_mat
    
    #stratum names for etiology regression:
    dput(unique_Eti_level,file.path(mcmc_options$result.folder,"unique_Eti_level.txt"))
    
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
          warning(paste0("==[baker] Bronze-standard slice ", names(data_nplcm$Mobs$MBS)[s], " has no measurements informative of the causes specified in 'cause_list', except 'NoA'! Please check if you need this measurement slice columns correspond to causes other than 'NoA'.=="))  
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
        stop("==[baker] Please put cases on top of controls in `data_nplcm`.==\n")
      }
      
      Z_FPR_list <- lapply(FPR_formula,function(form){stats::model.matrix(form,data.frame(X,Y))}) # <-- make sure that the row orders are the same.
      
      #intercept_only_MBS <- which(lapply(Z_FPR_list,function(x) (ncol(x)==1 & all(x==1)))==TRUE)
      
      unique_FPR_level_list <- lapply(Z_FPR_list, function(Z){
        unique(Z)
      })
      
      n_unique_FPR_level_list <- lapply(unique_FPR_level_list,nrow)
      
      FPR_stratum_id_list <- list()
      
      #Z_FPR_list_orth <- list()
      for(i in seq_along(JBrS_list)){
        FPR_stratum_id_list[[i]] <- apply(Z_FPR_list[[i]],1,function(v) 
          which(rowSums(abs(unique_FPR_level_list[[i]]-t(replicate(n_unique_FPR_level_list[[i]],v))))==0))
      }
      names(FPR_stratum_id_list) <- names(MBS_list)
      
      FPR_colname_design_mat_list        <- vector("list",length(Mobs$MBS))
      names(FPR_colname_design_mat_list) <-  names(Mobs$MBS)
      for (i in seq_along(JBrS_list)){
        assign(paste("unique_FPR_level",i,sep="_"),unique_FPR_level_list[[i]])
        assign(paste("n_unique_FPR_level",i,sep="_"),n_unique_FPR_level_list[[i]])
        assign(paste("FPR_stratum_id",i,sep="_"),unname(FPR_stratum_id_list[[i]])) # <-- caused problem when this is the only variable in the data_nplcm$X. So unname solves the problem.
        
        assign(paste("FPR_colname_design_mat", i, sep = "_"), attributes(Z_FPR_list[[i]])$dimnames[[2]]) 
        FPR_colname_design_mat_list[[i]] <- eval(parse(text = paste0("FPR_colname_design_mat_",i)))
        unique_FPR_level_list[[i]]       <- eval(parse(text = paste0("unique_FPR_level_",i)))
        
        rownames(unique_FPR_level_list[[i]]) <- 1:eval(parse(text = paste0("n_unique_FPR_level_",i)))
        colnames(unique_FPR_level_list[[i]]) <- FPR_colname_design_mat_list[[i]]
        
        attributes(Z_FPR_list[[i]])[names(attributes(Z_FPR_list[[i]]))!="dim"] <- NULL
        
        assign(paste("JBrS", i, sep = "_"), JBrS_list[[i]])    
        assign(paste("MBS", i, sep = "_"), as.matrix_or_vec(MBS_list[[i]])) 
        assign(paste("templateBS", i, sep = "_"), as.matrix_or_vec(template_BrS_list[[i]]))  
        
        # Z_FPR0       <- Z_FPR_list[[i]]
        # a.eig        <- eigen(t(Z_FPR0)%*%Z_FPR0)
        # sqrt_Z_FPR0  <- a.eig$vectors %*% diag(sqrt(a.eig$values)) %*% solve(a.eig$vectors)
        # 
        # Z_FPR_list_orth[[i]]  <- Z_FPR0%*%solve(sqrt_Z_FPR0)
        # assign(paste("Z_FPR",i,sep="_"),Z_FPR_list_orth[[i]])
        
      }
      
      #stratum names for FPR regression:
      dput(unique_FPR_level_list,file.path(mcmc_options$result.folder,"unique_FPR_level_list.txt"))
      
      
      # setup groupwise TPR for BrS:
      BrS_TPR_strat <- FALSE
      prior_BrS     <- model_options$prior$TPR_prior$BrS
      parsed_model <- assign_model(model_options,data_nplcm)
      if (parsed_model$BrS_grp){
        BrS_TPR_strat <- TRUE
        for(i in seq_along(JBrS_list)){
          assign(paste("GBrS_TPR", i, sep = "_"),  length(unique(prior_BrS$grp)))    
          assign(paste("BrS_TPR_grp", i, sep = "_"),  prior_BrS$grp)    
        }
      }
      
      # add GBrS_TPR_1, or 2 if we want to index by slices:
      for (i in seq_along(JBrS_list)){
        if (!is.null(prior_BrS$grp)){ # <--- need to change to list if we have multiple slices.
          assign(paste("GBrS_TPR", i, sep = "_"), length(unique(prior_BrS$grp))) # <--- need to change to depending on i if grp changes wrt specimen.
        }
        if (is.null(prior_BrS$grp)){ # <--- need to change to list if we have multiple slices.
          assign(paste("GBrS_TPR", i, sep = "_"), 1)
        }
      }
      
      # set BrS measurement priors: 
      # hyper-parameters for sensitivity:
      
      alpha_mat <- list() # dimension for slices.
      beta_mat  <- list()
      
      for(i in seq_along(JBrS_list)){
        
        if (likelihood$k_subclass[i] == 1){BrS_tpr_prior <- set_prior_tpr_BrS_NoNest(i,model_options,data_nplcm)}
        if (likelihood$k_subclass[i] > 1){BrS_tpr_prior <- set_prior_tpr_BrS_NoNest(i,model_options,data_nplcm)}
        
        GBrS_TPR_curr <- eval(parse(text = paste0("GBrS_TPR_",i)))
        alpha_mat[[i]] <- matrix(NA, nrow=GBrS_TPR_curr,ncol=JBrS_list[[i]])
        beta_mat[[i]]  <- matrix(NA, nrow=GBrS_TPR_curr,ncol=JBrS_list[[i]])
        
        colnames(alpha_mat[[i]]) <- patho_BrS_list[[i]]
        colnames(beta_mat[[i]]) <- patho_BrS_list[[i]]
        
        for (g in 1:GBrS_TPR_curr){
          alpha_mat[[i]][g,] <- unlist(BrS_tpr_prior[[1]][[g]]$alpha)
          beta_mat[[i]][g,]  <- unlist(BrS_tpr_prior[[1]][[g]]$beta)
        }
        
        if (GBrS_TPR_curr>1){
          assign(paste("alphaB", i, sep = "_"), alpha_mat[[i]])      # <---- input BrS TPR prior here.
          assign(paste("betaB", i, sep = "_"),  beta_mat[[i]])    
        }else{
          assign(paste("alphaB", i, sep = "_"), c(alpha_mat[[i]]))   # <---- input BrS TPR prior here.
          assign(paste("betaB", i, sep = "_"),  c(beta_mat[[i]]))    
        }
      }
      names(alpha_mat) <- names(beta_mat)<- names(Mobs$MBS)
      
      
      if (!BrS_TPR_strat){# no stratification.
        # summarize into one name (for all measurements):
        if (length(single_column_MBS)==0){
          # if all slices have >2 columns:
          in_data       <- c(in_data,"Nd","Nu","Jcause","n_unique_Eti_level","alphaEti",
                             paste("JBrS",1:length(JBrS_list),sep="_"),
                             paste("MBS",1:length(JBrS_list),sep="_"),
                             paste("templateBS",1:length(JBrS_list),sep="_"),
                             "Eti_stratum_id",
                             paste("n_unique_FPR_level",1:length(JBrS_list),sep="_"),
                             paste("FPR_stratum_id",1:length(JBrS_list),sep="_")
                             # paste("alphaB",1:length(JBrS_list),sep="_"),
                             # paste("betaB",1:length(JBrS_list),sep="_")
          )
        } else {
          # if there exist slices with 1 column:
          in_data       <- c(in_data,"Nd","Nu","Jcause","n_unique_Eti_level","alphaEti",
                             paste("JBrS",1:length(JBrS_list),sep="_")[-single_column_MBS], # <---- no need to iterate in .bug file for a slice with one column.
                             paste("MBS",1:length(JBrS_list),sep="_"),
                             paste("templateBS",1:length(JBrS_list),sep="_"),
                             "Eti_stratum_id",
                             paste("n_unique_FPR_level",1:length(JBrS_list),sep="_"),
                             paste("FPR_stratum_id",1:length(JBrS_list),sep="_")
                             # paste("alphaB",1:length(JBrS_list),sep="_"),
                             # paste("betaB",1:length(JBrS_list),sep="_")
          )
        }
      } else{ # with stratification.
        # summarize into one name (for all measurements):
        if (length(single_column_MBS)==0){
          # if all slices have >2 columns:
          in_data       <- c(in_data,"Nd","Nu","Jcause","n_unique_Eti_level","alphaEti",
                             paste("JBrS",1:length(JBrS_list),sep="_"),
                             paste("GBrS_TPR",1:length(JBrS_list),sep="_"),   # <-- added for TPR strata.
                             paste("BrS_TPR_grp",1:length(JBrS_list),sep="_"),# <-- added for TPR strata.
                             paste("MBS",1:length(JBrS_list),sep="_"),
                             paste("templateBS",1:length(JBrS_list),sep="_"),
                             "Eti_stratum_id",
                             paste("n_unique_FPR_level",1:length(JBrS_list),sep="_"),
                             paste("FPR_stratum_id",1:length(JBrS_list),sep="_")
                             # paste("alphaB",1:length(JBrS_list),sep="_"),
                             # paste("betaB",1:length(JBrS_list),sep="_")
          )
        } else {
          # if there exist slices with 1 column:
          in_data       <- c(in_data,"Nd","Nu","Jcause","n_unique_Eti_level","alphaEti",
                             paste("JBrS",1:length(JBrS_list),sep="_")[-single_column_MBS], # <---- no need to iterate in .bug file for a slice with one column.
                             paste("GBrS_TPR",1:length(JBrS_list),sep="_"),   # <-- added for TPR strata.
                             paste("BrS_TPR_grp",1:length(JBrS_list),sep="_"),# <-- added for TPR strata.
                             paste("MBS",1:length(JBrS_list),sep="_"),
                             paste("templateBS",1:length(JBrS_list),sep="_"),
                             "Eti_stratum_id",
                             paste("n_unique_FPR_level",1:length(JBrS_list),sep="_"),
                             paste("FPR_stratum_id",1:length(JBrS_list),sep="_")
                             # paste("alphaB",1:length(JBrS_list),sep="_"),
                             # paste("betaB",1:length(JBrS_list),sep="_")
          )
        }
      }
      
      #
      # hyper-parameters:
      #
      
      # Set BrS measurement priors:
      # hyperparameter for sensitivity (can add for specificity if necessary): 
      for (s in seq_along(Mobs$MBS)){
        #if (likelihood$k_subclass[s] == 1){BrS_tpr_prior <- set_prior_tpr_BrS_NoNest(s,model_options,data_nplcm)}
        #if (likelihood$k_subclass[s] > 1){BrS_tpr_prior <- set_prior_tpr_BrS_NoNest(s,model_options,data_nplcm)}
        
        #assign(paste("alphaB", s, sep = "_"), BrS_tpr_prior[[1]]$alpha)     # <---- input BrS TPR prior here.
        #assign(paste("betaB", s, sep = "_"),  BrS_tpr_prior[[1]]$beta)    
        
        if (likelihood$k_subclass[s]==1){
          out_parameter <- c(out_parameter,paste(c("thetaBS","psiBS"), s, sep="_"))
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
          in_data       <- unique(c(in_data,"Nd","Jcause",
                                    paste("JSS",1:length(JSS_list),sep="_"),
                                    paste("MSS",1:length(JSS_list),sep="_"),
                                    paste("templateSS",1:length(JSS_list),sep="_"),
                                    paste("alphaS",1:length(JSS_list),sep="_"),   
                                    paste("betaS",1:length(JSS_list),sep="_")))
        } else{
          in_data       <- unique(c(in_data,"Nd","Nu","Jcause",
                                    paste("JSS",1:length(JSS_list),sep="_")[-single_column_MSS],
                                    paste("MSS",1:length(JSS_list),sep="_"),
                                    paste("templateSS",1:length(JSS_list),sep="_"),
                                    paste("alphaS",1:length(JSS_list),sep="_"),
                                    paste("betaS",1:length(JSS_list),sep="_")))
        }
      }else {
        if (length(single_column_MSS)==0){
          # summarize into one name (for all measurements):
          in_data       <- unique(c(in_data,"Nd","Jcause",
                                    paste("JSS",1:length(JSS_list),sep="_"),
                                    paste("GSS_TPR",1:length(JSS_list),sep="_"),
                                    paste("SS_TPR_grp",1:length(JSS_list),sep="_"),
                                    paste("MSS",1:length(JSS_list),sep="_"),
                                    paste("templateSS",1:length(JSS_list),sep="_"),
                                    paste("alphaS",1:length(JSS_list),sep="_"),   
                                    paste("betaS",1:length(JSS_list),sep="_")))
        } else{
          in_data       <- unique(c(in_data,"Nd","Nu","Jcause",
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
            #res_curr[[1]] <- stats::rbeta(JBrS_list[[s]],1,1)
            GBrS_TPR_curr <- eval(parse(text = paste0("GBrS_TPR_",s)))
            if (GBrS_TPR_curr==1){
              res_curr[[1]] <- stats::rbeta(JBrS_list[[s]],1,1)
            } else{
              res_curr[[1]] <- matrix(stats::rbeta(GBrS_TPR_curr*JBrS_list[[s]],1,1),
                                      nrow=GBrS_TPR_curr,ncol=JBrS_list[[s]])
              if (JBrS_list[[s]]==1){
                res_curr[[1]] <- c(res_curr[[1]])
              }
            }
            res_curr[[2]] <- matrix(stats::rbeta(n_unique_FPR_level_list[[s]]*JBrS_list[[s]],1,1),
                                    nrow=n_unique_FPR_level_list[[s]],ncol=JBrS_list[[s]])
            names(res_curr) <- paste(c("thetaBS","psiBS"),s,sep="_")
            
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
        
        tmp <- matrix(stats::rgamma(n_unique_Eti_level*Jcause, 1,1),nrow=n_unique_Eti_level,ncol=Jcause)
        res[[length(res)+1]] <- diag(1/rowSums(tmp),nrow=nrow(tmp),ncol=nrow(tmp))%*%tmp
        names(res)[length(res)] <- "pEti"
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
            if (JSS_list[[i]]==1){
              tmp_thetaSS[[i]] <- c(tmp_thetaSS[[i]])
            }
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
            #res_curr[[1]] <- stats::rbeta(JBrS_list[[s]],1,1)
            GBrS_TPR_curr <- eval(parse(text = paste0("GBrS_TPR_",s)))
            if (GBrS_TPR_curr==1){
              res_curr[[1]] <- stats::rbeta(JBrS_list[[s]],1,1)
            } else{
              res_curr[[1]] <- matrix(stats::rbeta(GBrS_TPR_curr*JBrS_list[[s]],1,1),
                                      nrow=GBrS_TPR_curr,ncol=JBrS_list[[s]])
              if (JBrS_list[[s]]==1){
                res_curr[[1]] <- c(res_curr[[1]])
              }
            }
            res_curr[[2]] <- matrix(stats::rbeta(n_unique_FPR_level_list[[s]]*JBrS_list[[s]],1,1),
                                    nrow=n_unique_FPR_level_list[[s]],ncol=JBrS_list[[s]])
            names(res_curr) <- paste(c("thetaBS","psiBS"),s,sep="_")
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
            if (JSS_list[[i]]==1){
              tmp_thetaSS[[i]] <- c(tmp_thetaSS[[i]])
            }
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
        tmp <- matrix(stats::rgamma(n_unique_Eti_level*Jcause, 1,1),nrow=n_unique_Eti_level,ncol=Jcause)
        res[[length(res)+1]] <- diag(1/rowSums(tmp),nrow=nrow(tmp), ncol=nrow(tmp))%*%tmp
        names(res)[length(res)] <- "pEti"
        
        c(res,res2)
      }
    }
    
    # etiology (measurement independent)
    alphaEti          <- prior$Eti_prior    # <-------- input etiology prior here.
    attributes(alphaEti)[names(attributes(alphaEti))!="dim"]<- NULL # the dimnames caused issues.
    
    if(length(c(alphaEti)) < Jcause*n_unique_Eti_level){stop("==[baker] You've specified stratified etiologies, 
          one per stratum defined by 'Eti_formula' in 'model_options$likelihood'. Please specify a matrix of 'Eti_prior' with 
          (#rows, #columns) = (#Etiology strata, #causes). Also be aware that etiology prior for each stratum corresponds to pseudo-observations,
          and they altogether comprise the total pseudo-observations. For example, if you specify a matrix of 1s for 10 strata and 30 causes,
          it says you observed in a priori in each stratum 30 individuals, each caused by one of the causes!===")}
    
    #
    # fit model :
    #
    
    # special operations:
    # get individual prediction:
    if (!is.null(mcmc_options$individual.pred) && mcmc_options$individual.pred){out_parameter <- c(out_parameter,"Icat")}
    # get posterior predictive distribution of BrS measurments:
    if (!is.null(mcmc_options$ppd) && mcmc_options$ppd){out_parameter <- c(out_parameter,paste("MBS.new",seq_along(Mobs$MBS),sep = "_"))}
    # get pEti: (for regression models; besides the coefficients of etiology regressions are always recorded)
    if (!is.null(mcmc_options$get.pEti) && mcmc_options$get.pEti){out_parameter <- unique(c(out_parameter,"pEti"))}
    
    #
    # write the .bug files into mcmc_options$bugsmodel.dir; 
    # could later set it equal to result.folder.
    #
    use_jags <- (!is.null(mcmc_options$use_jags) && mcmc_options$use_jags)
    model_func         <- write_model_Reg_discrete_predictor_NoNest(data_nplcm$Mobs,
                                                                    model_options$prior,
                                                                    model_options$likelihood$cause_list,
                                                                    model_options$use_measurements,
                                                                    mcmc_options$ppd,
                                                                    use_jags)
    model_bugfile_name <- "model_Reg_discrete_predictor_NoNest.bug"
    
    filename <- file.path(mcmc_options$bugsmodel.dir, model_bugfile_name)
    writeLines(model_func, filename)
    
    in_data <- unique(c(in_data))
    
    # #
    # # # uncomment below if using dmnorm for betaEti and betaFPR (also need to edit plug-and-play.R):
    # #
    # # if (Jcause > 2){ # use dmnorm to speed up JAGS calculations, so we need mean 
    # #   # and covariance vector.
    # #   zero_Jcause_1 <- rep(0,Jcause-1)
    # #   I_Jcause_1    <- diag(1,Jcause-1)
    # #   in_data <- c(in_data,"zero_Jcause_1","I_Jcause_1")
    # # } 
    # 
    # # for (s in seq_along(JBrS_list)){
    # #   if (JBrS_list[[s]]>1){
    # #       assign(paste("zero_JBrS", s, sep = "_"), rep(0,JBrS_list[[s]]))    
    # #       assign(paste("I_JBrS", s, sep = "_"), diag(1,JBrS_list[[s]]))    
    # #       in_data <- c(in_data,
    # #                    paste("zero_JBrS", s, sep = "_"),
    # #                    paste("I_JBrS", s, sep = "_"))
    # #   }
    # # }
    
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
      curr_data_txt_file <- file.path(mcmc_options$result.folder,"jagsdata.txt")
      if(file.exists(curr_data_txt_file)){file.remove(curr_data_txt_file)}
      dump(names(in_data.list), append = FALSE, envir = here,
           file = curr_data_txt_file)
      # fix dimension problem.... convert say .Dmi=7:6 to c(7,6) (an issue for templateBS_1):
      bad_jagsdata_txt <- readLines(curr_data_txt_file)
      good_jagsdata_txt <- gsub( ".Dim = ([0-9]+):([0-9]+)", ".Dim = c(\\1,\\2)", bad_jagsdata_txt,fixed = FALSE)
      writeLines(good_jagsdata_txt, curr_data_txt_file)
      
      
      
      gs <- jags2_baker(data   = curr_data_txt_file,
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





#' Fit nested partially-latent class model with regression (low-level)
#'
#' @details This function prepares data, specifies hyperparameters in priors 
#' (true positive rates and etiology fractions), initializes the posterior
#' sampling chain, writes the model file (for JAGS or WinBUGS with slight
#' differences in syntax), and fits the model. Features:
#' \itemize{
#' \item regression (not all discrete covariates);
#' \item no nested subclasses, i.e. conditional independence of 
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
    # Z_Eti0       <- stats::model.matrix(Eti_formula,data.frame(X,Y)[Y==1,,drop=FALSE])
    # a.eig        <- eigen(t(Z_Eti0)%*%Z_Eti0)
    # sqrt_Z_Eti0  <- a.eig$vectors %*% diag(sqrt(a.eig$values)) %*% solve(a.eig$vectors)
    # 
    # Z_Eti  <- Z_Eti0%*%solve(sqrt_Z_Eti0)
    
    ncol_dm_Eti <- ncol(Z_Eti)
    # if need to do PS, need to insert many basis operations here: #<-----!
    ER_has_basis <- ifelse(length(grep("^s_",dimnames(Z_Eti)[[2]]))==0,FALSE,TRUE)
    ER_has_non_basis <- has_non_basis(Eti_formula)
    ER_basis_id <- c(sapply(1, function(s) {
      if (length(grep("^s_",dimnames(Z_Eti)[[2]]))==0){
        NULL
      } else{
        0+grep("^s_",dimnames(Z_Eti)[[2]])
      }
    }))
    ER_n_basis  <- length(ER_basis_id)
    ER_non_basis_id <- c(sapply(1, function(s) {
      if (length(grep("^s_",dimnames(Z_Eti)[[2]]))==0){
        0+(1:ncol(Z_Eti))
      } else{
        (1:ncol(Z_Eti))[-grep("^s_",dimnames(Z_Eti)[[2]])]
      }
    }))
    #END <-----!
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
      has_non_basis_list <- lapply(Z_FPR_list, function(Z) {
        if (length(grep("^s_",dimnames(Z)[[2]]))==ncol(Z)){
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
      
      #Z_FPR_list_orth <- list()
      for(i in seq_along(JBrS_list)){
        attributes(Z_FPR_list[[i]])[names(attributes(Z_FPR_list[[i]]))!="dim"] <- NULL
        
        assign(paste("JBrS", i, sep = "_"), JBrS_list[[i]])    
        assign(paste("MBS", i, sep = "_"), as.matrix_or_vec(MBS_list[[i]])) 
        assign(paste("templateBS", i, sep = "_"), as.matrix_or_vec(template_BrS_list[[i]]))  
        
        # Z_FPR0       <- Z_FPR_list[[i]]
        # a.eig        <- eigen(t(Z_FPR0)%*%Z_FPR0)
        # sqrt_Z_FPR0  <- a.eig$vectors %*% diag(sqrt(a.eig$values)) %*% solve(a.eig$vectors)
        # 
        # Z_FPR_list_orth[[i]]  <- Z_FPR0%*%solve(sqrt_Z_FPR0)
        # assign(paste("Z_FPR",i,sep="_"),Z_FPR_list_orth[[i]])
        assign(paste("Z_FPR",i,sep="_"),Z_FPR_list[[i]])
        assign(paste("basis_id",i,sep="_"),basis_id_list[[i]])
        assign(paste("n_basis",i,sep="_"),n_basis_list[[i]])
        assign(paste("non_basis_id",i,sep="_"),non_basis_id_list[[i]])
      }
      
      # setup groupwise TPR for BrS:
      BrS_TPR_strat <- FALSE
      prior_BrS     <- model_options$prior$TPR_prior$BrS
      parsed_model <- assign_model(model_options,data_nplcm)
      if (parsed_model$BrS_grp){
        BrS_TPR_strat <- TRUE
        for(i in seq_along(JBrS_list)){
          assign(paste("GBrS_TPR", i, sep = "_"),  length(unique(prior_BrS$grp)))    
          assign(paste("BrS_TPR_grp", i, sep = "_"),  prior_BrS$grp)    
        }
      }
      
      # add GBrS_TPR_1, or 2 if we want to index by slices:
      for (i in seq_along(JBrS_list)){
        if (!is.null(prior_BrS$grp)){ # <--- need to change to list if we have multiple slices.
          assign(paste("GBrS_TPR", i, sep = "_"), length(unique(prior_BrS$grp))) # <--- need to change to depending on i if grp change wrt specimen.
        }
        if (is.null(prior_BrS$grp)){ # <--- need to change to list if we have multiple slices.
          assign(paste("GBrS_TPR", i, sep = "_"), 1)
        }
      }
      
      # set BrS measurement priors: 
      # hyper-parameters for sensitivity:
      
      alpha_mat <- list() # dimension for slices.
      beta_mat  <- list()
      
      for(i in seq_along(JBrS_list)){
        
        if (likelihood$k_subclass[i] == 1){BrS_tpr_prior <- set_prior_tpr_BrS_NoNest(i,model_options,data_nplcm)}
        if (likelihood$k_subclass[i] > 1){BrS_tpr_prior <- set_prior_tpr_BrS_NoNest(i,model_options,data_nplcm)}
        
        GBrS_TPR_curr <- eval(parse(text = paste0("GBrS_TPR_",i)))
        alpha_mat[[i]] <- matrix(NA, nrow=GBrS_TPR_curr,ncol=JBrS_list[[i]])
        beta_mat[[i]]  <- matrix(NA, nrow=GBrS_TPR_curr,ncol=JBrS_list[[i]])
        
        colnames(alpha_mat[[i]]) <- patho_BrS_list[[i]]
        colnames(beta_mat[[i]]) <- patho_BrS_list[[i]]
        
        for (g in 1:GBrS_TPR_curr){
          alpha_mat[[i]][g,] <- unlist(BrS_tpr_prior[[1]][[g]]$alpha)
          beta_mat[[i]][g,]  <- unlist(BrS_tpr_prior[[1]][[g]]$beta)
        }
        
        if (GBrS_TPR_curr>1){
          assign(paste("alphaB", i, sep = "_"), alpha_mat[[i]])      # <---- input BrS TPR prior here.
          assign(paste("betaB", i, sep = "_"),  beta_mat[[i]])    
        }else{
          assign(paste("alphaB", i, sep = "_"), c(alpha_mat[[i]]))   # <---- input BrS TPR prior here.
          assign(paste("betaB", i, sep = "_"),  c(beta_mat[[i]]))    
        }
      }
      names(alpha_mat) <- names(beta_mat)<- names(Mobs$MBS)
      
      
      if (!BrS_TPR_strat){# no stratification.
        # summarize into one name (for all measurements):
        if (length(single_column_MBS)==0){
          # if all slices have >2 columns:
          in_data       <- c(in_data,"Nd","Nu","Jcause","Z_Eti",
                             paste("JBrS",1:length(JBrS_list),sep="_"),
                             paste("MBS",1:length(JBrS_list),sep="_"),
                             paste("templateBS",1:length(JBrS_list),sep="_"),
                             paste("Z_FPR",1:length(JBrS_list),sep="_"),
                             paste("basis_id",1:length(JBrS_list),sep="_")[unlist(has_basis_list)],
                             paste("n_basis",1:length(JBrS_list),sep="_")[unlist(has_basis_list)],
                             paste("non_basis_id",1:length(JBrS_list),sep="_")[unlist(has_non_basis_list)]
                             # paste("alphaB",1:length(JBrS_list),sep="_"),
                             # paste("betaB",1:length(JBrS_list),sep="_")
          )
        } else {
          # if there exist slices with 1 column:
          in_data       <- c(in_data,"Nd","Nu","Jcause","Z_Eti",
                             paste("JBrS",1:length(JBrS_list),sep="_")[-single_column_MBS], # <---- no need to iterate in .bug file for a slice with one column.
                             paste("MBS",1:length(JBrS_list),sep="_"),
                             paste("templateBS",1:length(JBrS_list),sep="_"),
                             paste("Z_FPR",1:length(JBrS_list),sep="_"),
                             paste("basis_id",1:length(JBrS_list),sep="_")[unlist(has_basis_list)],
                             paste("n_basis",1:length(JBrS_list),sep="_")[unlist(has_basis_list)],
                             paste("non_basis_id",1:length(JBrS_list),sep="_")[unlist(has_non_basis_list)]
                             # paste("alphaB",1:length(JBrS_list),sep="_"),
                             # paste("betaB",1:length(JBrS_list),sep="_")
          )
        }
      } else{ # with stratification.
        # summarize into one name (for all measurements):
        if (length(single_column_MBS)==0){
          # if all slices have >2 columns:
          in_data       <- c(in_data,"Nd","Nu","Jcause","Z_Eti",
                             paste("JBrS",1:length(JBrS_list),sep="_"),
                             paste("GBrS_TPR",1:length(JBrS_list),sep="_"),   # <-- added for TPR strata.
                             paste("BrS_TPR_grp",1:length(JBrS_list),sep="_"),# <-- added for TPR strata.
                             paste("MBS",1:length(JBrS_list),sep="_"),
                             paste("templateBS",1:length(JBrS_list),sep="_"),
                             paste("Z_FPR",1:length(JBrS_list),sep="_"),
                             paste("basis_id",1:length(JBrS_list),sep="_")[unlist(has_basis_list)],
                             paste("n_basis",1:length(JBrS_list),sep="_")[unlist(has_basis_list)],
                             paste("non_basis_id",1:length(JBrS_list),sep="_")[unlist(has_non_basis_list)]
                             # paste("alphaB",1:length(JBrS_list),sep="_"),
                             # paste("betaB",1:length(JBrS_list),sep="_")
          )
        } else {
          # if there exist slices with 1 column:
          in_data       <- c(in_data,"Nd","Nu","Jcause","Z_Eti",
                             paste("JBrS",1:length(JBrS_list),sep="_")[-single_column_MBS], # <---- no need to iterate in .bug file for a slice with one column.
                             paste("GBrS_TPR",1:length(JBrS_list),sep="_"),   # <-- added for TPR strata.
                             paste("BrS_TPR_grp",1:length(JBrS_list),sep="_"),# <-- added for TPR strata.
                             paste("MBS",1:length(JBrS_list),sep="_"),
                             paste("templateBS",1:length(JBrS_list),sep="_"),
                             paste("Z_FPR",1:length(JBrS_list),sep="_"),
                             paste("basis_id",1:length(JBrS_list),sep="_")[unlist(has_basis_list)],
                             paste("n_basis",1:length(JBrS_list),sep="_")[unlist(has_basis_list)],
                             paste("non_basis_id",1:length(JBrS_list),sep="_")[unlist(has_non_basis_list)]
                             # paste("alphaB",1:length(JBrS_list),sep="_"),
                             # paste("betaB",1:length(JBrS_list),sep="_")
          )
        }
      }
      
      #
      # hyper-parameters:
      #
      
      # Set BrS measurement priors:
      # hyperparameter for sensitivity (can add for specificity if necessary): 
      for (s in seq_along(Mobs$MBS)){
        #if (likelihood$k_subclass[s] == 1){BrS_tpr_prior <- set_prior_tpr_BrS_NoNest(s,model_options,data_nplcm)}
        #if (likelihood$k_subclass[s] > 1){BrS_tpr_prior <- set_prior_tpr_BrS_NoNest(s,model_options,data_nplcm)}
        
        #assign(paste("alphaB", s, sep = "_"), BrS_tpr_prior[[1]]$alpha)     # <---- input BrS TPR prior here.
        #assign(paste("betaB", s, sep = "_"),  BrS_tpr_prior[[1]]$beta)    
        
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
      if (ER_has_basis){
        out_parameter <- c(out_parameter,"betaEti","ER_taubeta") #<-- may not have this parameter in ncs.
      } else{
        out_parameter <- c(out_parameter,"betaEti")
      }
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
          in_data       <- unique(c(in_data,"Nd","Jcause",
                                    paste("JSS",1:length(JSS_list),sep="_"),
                                    paste("MSS",1:length(JSS_list),sep="_"),
                                    paste("templateSS",1:length(JSS_list),sep="_"),
                                    paste("alphaS",1:length(JSS_list),sep="_"),   
                                    paste("betaS",1:length(JSS_list),sep="_")))
        } else{
          in_data       <- unique(c(in_data,"Nd","Nu","Jcause",
                                    paste("JSS",1:length(JSS_list),sep="_")[-single_column_MSS],
                                    paste("MSS",1:length(JSS_list),sep="_"),
                                    paste("templateSS",1:length(JSS_list),sep="_"),
                                    paste("alphaS",1:length(JSS_list),sep="_"),
                                    paste("betaS",1:length(JSS_list),sep="_")))
        }
      }else {
        if (length(single_column_MSS)==0){
          # summarize into one name (for all measurements):
          in_data       <- unique(c(in_data,"Nd","Jcause",
                                    paste("JSS",1:length(JSS_list),sep="_"),
                                    paste("GSS_TPR",1:length(JSS_list),sep="_"),
                                    paste("SS_TPR_grp",1:length(JSS_list),sep="_"),
                                    paste("MSS",1:length(JSS_list),sep="_"),
                                    paste("templateSS",1:length(JSS_list),sep="_"),
                                    paste("alphaS",1:length(JSS_list),sep="_"),   
                                    paste("betaS",1:length(JSS_list),sep="_")))
        } else{
          in_data       <- unique(c(in_data,"Nd","Nu","Jcause",
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
    # 2. set up initialization function:
    #####################################
    if ("BrS" %in% use_measurements & !("SS" %in% use_measurements)){
      in_init       <-   function(){
        res <- list()
        for (s in seq_along(Mobs$MBS)){
          res_curr <- list()
          if (likelihood$k_subclass[s]==1){
            #res_curr[[1]] <- stats::rbeta(JBrS_list[[s]],1,1)
            GBrS_TPR_curr <- eval(parse(text = paste0("GBrS_TPR_",s)))
            if (GBrS_TPR_curr==1){
              res_curr[[1]] <- stats::rbeta(JBrS_list[[s]],1,1)
            } else{
              res_curr[[1]] <- matrix(stats::rbeta(GBrS_TPR_curr*JBrS_list[[s]],1,1),
                                      nrow=GBrS_TPR_curr,ncol=JBrS_list[[s]])
              if (JBrS_list[[s]]==1){
                res_curr[[1]] <- c(res_curr[[1]])
              }
            }
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
            if (JSS_list[[i]]==1){
              tmp_thetaSS[[i]] <- c(tmp_thetaSS[[i]])
            }
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
            #res_curr[[1]] <- stats::rbeta(JBrS_list[[s]],1,1)
            GBrS_TPR_curr <- eval(parse(text = paste0("GBrS_TPR_",s)))
            if (GBrS_TPR_curr==1){
              res_curr[[1]] <- stats::rbeta(JBrS_list[[s]],1,1)
            } else{
              res_curr[[1]] <- matrix(stats::rbeta(GBrS_TPR_curr*JBrS_list[[s]],1,1),
                                      nrow=GBrS_TPR_curr,ncol=JBrS_list[[s]])
              if (JBrS_list[[s]]==1){
                res_curr[[1]] <- c(res_curr[[1]])
              }
            }
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
            if (JSS_list[[i]]==1){
              tmp_thetaSS[[i]] <- c(tmp_thetaSS[[i]])
            }
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
    if (!is.null(mcmc_options$get.pEti) && mcmc_options$get.pEti){out_parameter <- unique(c(out_parameter,"pEti"))}
    
    # write the .bug files into mcmc_options$bugsmodel.dir; 
    # could later set it equal to result.folder.
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
    
    # if (length(prior$Eti_prior)>1){
    #   stop("== [baker] Regression model used. Please change `model_options$prior$Eti_prior` 
    #        into a single, positive real number representing
    #        the stanadrd deviation of beta coefficients in etiology regression! ==")
    # }
    sd_betaEti_basis    <- prior$Eti_prior[1]
    sd_betaEti_nonbasis <- prior$Eti_prior[2]   #<--place to adjust (currently identical for basis vs non-basis. check "insert_bugfile_chunk_reg_etiology")
    sd_betaFPR_basis   <- prior$FPR_coef_prior[1]
    sd_betaFPR_nonbasis <- prior$FPR_coef_prior[2]
    ER_alpha    <- prior$Eti_hyper_pflex[1]
    ER_beta     <- prior$Eti_hyper_pflex[2]
    
    in_data <- unique(c(in_data,
                        "ER_basis_id"[ER_has_basis],
                        "ER_n_basis"[ER_has_basis],
                        "ER_non_basis_id"[ER_has_non_basis],
                        "sd_betaEti_basis"[ER_has_basis],
                        "ER_alpha"[ER_has_basis],
                        "ER_beta"[ER_has_basis],
                        "sd_betaEti_nonbasis"[ER_has_non_basis],
                        "sd_betaFPR_basis"[any(unlist(has_basis_list))],
                        "sd_betaFPR_nonbasis"[any(unlist(has_non_basis_list))]
    ))
    
    # #
    # # # uncomment below if using dmnorm for betaEti and betaFPR (also need to edit plug-and-play.R):
    # #
    # # if (Jcause > 2){ # use dmnorm to speed up JAGS calculations, so we need mean 
    # #   # and covariance vector.
    # #   zero_Jcause_1 <- rep(0,Jcause-1)
    # #   I_Jcause_1    <- diag(1,Jcause-1)
    # #   in_data <- c(in_data,"zero_Jcause_1","I_Jcause_1")
    # # } 
    # 
    # # for (s in seq_along(JBrS_list)){
    # #   if (JBrS_list[[s]]>1){
    # #       assign(paste("zero_JBrS", s, sep = "_"), rep(0,JBrS_list[[s]]))    
    # #       assign(paste("I_JBrS", s, sep = "_"), diag(1,JBrS_list[[s]]))    
    # #       in_data <- c(in_data,
    # #                    paste("zero_JBrS", s, sep = "_"),
    # #                    paste("I_JBrS", s, sep = "_"))
    # #   }
    # # }
    
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
      curr_data_txt_file <- file.path(mcmc_options$result.folder,"jagsdata.txt")
      if(file.exists(curr_data_txt_file)){file.remove(curr_data_txt_file)}
      dump(names(in_data.list), append = FALSE, envir = here,
           file = curr_data_txt_file)
      ## fix dimension problem.... convert say .Dmi=7:6 to c(7,6) (an issue for templateBS_1):
      bad_jagsdata_txt <- readLines(curr_data_txt_file)
      good_jagsdata_txt <- gsub( "([0-9]+):([0-9]+)", "c(\\1,\\2)", bad_jagsdata_txt,fixed = FALSE)
      writeLines(good_jagsdata_txt, curr_data_txt_file)
      
      # # fix dimension problem.... convert say 7:6 to c(7,6) (an issue for a dumped matrix):
      # inits_fnames <- list.files(mcmc_options$result.folder,pattern = "^jagsinits[0-9]+.txt",
      #                            full.names = TRUE)
      # for (fiter in seq_along(inits_fnames)){
      #   curr_inits_txt_file <- inits_fnames[fiter]
      #   bad_jagsinits_txt <- readLines(curr_inits_txt_file)
      #   good_jagsinits_txt <- gsub( "([0-9]+):([0-9]+)", "c(\\1,\\2)", bad_jagsinits_txt,fixed = FALSE)
      #   writeLines(good_jagsinits_txt, curr_inits_txt_file)
      # }
      
      gs <- jags2_baker(data   = curr_data_txt_file,
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




#' Fit nested partially-latent class model with regression (low-level)
#'
#' Called by \link{nplcm}() upon being assigned to this nested regression by 
#' \link{assign_model}
#'
#' @details This function prepares data, specifies hyperparameters in priors 
#' (true positive rates and etiology fractions), initializes the posterior
#' sampling chain, writes the model file (for JAGS or WinBUGS with slight
#' differences in syntax), and fits the model. Features:
#' \itemize{
#' \item regression (not all discrete covariates);
#' \item nested subclasses, i.e. conditional dependence of 
#' multivariate measurements given disease class and covariates;
#' \item multiple BrS + multiple SS.
#' }
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

nplcm_fit_Reg_Nest <- function(data_nplcm,model_options,mcmc_options){
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
  
  for (s in seq_along(Mobs$MBS)){
    if (likelihood$k_subclass[s]>1 && ncol(Mobs$MBS[[s]])==1){
      stop(paste0("==[baker] slice ", s, " bronze-standard data: Cannot do nested modeling for BrS measurements with only one column! =="))
    }
  }
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
  # Z_Eti0       <- stats::model.matrix(Eti_formula,data.frame(X,Y)[Y==1,,drop=FALSE])
  # a.eig        <- eigen(t(Z_Eti0)%*%Z_Eti0)
  # sqrt_Z_Eti0  <- a.eig$vectors %*% diag(sqrt(a.eig$values)) %*% solve(a.eig$vectors)
  # 
  # Z_Eti  <- Z_Eti0%*%solve(sqrt_Z_Eti0)
  
  ncol_dm_Eti <- ncol(Z_Eti)
  # if need to do PS, need to insert many basis operations here: #<-----!
  ER_has_basis <- ifelse(length(grep("^s_",dimnames(Z_Eti)[[2]]))==0,FALSE,TRUE)
  ER_has_non_basis <- has_non_basis(Eti_formula)
  ER_basis_id <- c(sapply(1, function(s) {
    if (length(grep("^s_",dimnames(Z_Eti)[[2]]))==0){
      NULL
    } else{
      0+grep("^s_",dimnames(Z_Eti)[[2]])
    }
  }))
  ER_n_basis  <- length(ER_basis_id)
  ER_non_basis_id <- c(sapply(1, function(s) {
    if (length(grep("^s_",dimnames(Z_Eti)[[2]]))==0){
      0+(1:ncol(Z_Eti))
    } else{
      (1:ncol(Z_Eti))[-grep("^s_",dimnames(Z_Eti)[[2]])]
    }
  }))
  #END <-----!
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
    has_non_basis_list <- lapply(Z_FPR_list, function(Z) {
      if (length(grep("^s_",dimnames(Z)[[2]]))==ncol(Z)){
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
    
    #Z_FPR_list_orth <- list()
    for(i in seq_along(JBrS_list)){
      attributes(Z_FPR_list[[i]])[names(attributes(Z_FPR_list[[i]]))!="dim"] <- NULL
      
      assign(paste("JBrS", i, sep = "_"), JBrS_list[[i]])    
      assign(paste("MBS", i, sep = "_"), as.matrix_or_vec(MBS_list[[i]])) 
      assign(paste("templateBS", i, sep = "_"), as.matrix_or_vec(template_BrS_list[[i]]))  
      
      # Z_FPR0       <- Z_FPR_list[[i]]
      # a.eig        <- eigen(t(Z_FPR0)%*%Z_FPR0)
      # sqrt_Z_FPR0  <- a.eig$vectors %*% diag(sqrt(a.eig$values)) %*% solve(a.eig$vectors)
      # 
      # Z_FPR_list_orth[[i]]  <- Z_FPR0%*%solve(sqrt_Z_FPR0)
      # assign(paste("Z_FPR",i,sep="_"),Z_FPR_list_orth[[i]])
      assign(paste("Z_FPR",i,sep="_"),Z_FPR_list[[i]])
      assign(paste("basis_id",i,sep="_"),basis_id_list[[i]])
      assign(paste("n_basis",i,sep="_"),n_basis_list[[i]])
      assign(paste("non_basis_id",i,sep="_"),non_basis_id_list[[i]])
    }
    
    # setup groupwise TPR for BrS:
    BrS_TPR_strat <- FALSE
    prior_BrS     <- model_options$prior$TPR_prior$BrS
    parsed_model <- assign_model(model_options,data_nplcm)
    if (parsed_model$BrS_grp){
      BrS_TPR_strat <- TRUE
      for(i in seq_along(JBrS_list)){
        assign(paste("GBrS_TPR", i, sep = "_"),  length(unique(prior_BrS$grp)))    
        assign(paste("BrS_TPR_grp", i, sep = "_"),  prior_BrS$grp)    
      }
    }
    
    # add GBrS_TPR_1, or 2 if we want to index by slices:
    for (i in seq_along(JBrS_list)){
      if (!is.null(prior_BrS$grp)){ # <--- need to change to list if we have multiple slices.
        assign(paste("GBrS_TPR", i, sep = "_"), length(unique(prior_BrS$grp))) # <--- need to change to depending on i if grp change wrt specimen.
      }
      if (is.null(prior_BrS$grp)){ # <--- need to change to list if we have multiple slices.
        assign(paste("GBrS_TPR", i, sep = "_"), 1)
      }
    }
    
    # set BrS measurement priors: 
    # hyper-parameters for sensitivity:
    
    alpha_mat <- list() # dimension for slices.
    beta_mat  <- list()
    
    for(i in seq_along(JBrS_list)){
      
      if (likelihood$k_subclass[i] == 1){BrS_tpr_prior <- set_prior_tpr_BrS_NoNest(i,model_options,data_nplcm)}
      if (likelihood$k_subclass[i] > 1){BrS_tpr_prior <- set_prior_tpr_BrS_NoNest(i,model_options,data_nplcm)}
      
      GBrS_TPR_curr <- eval(parse(text = paste0("GBrS_TPR_",i)))
      alpha_mat[[i]] <- matrix(NA, nrow=GBrS_TPR_curr,ncol=JBrS_list[[i]])
      beta_mat[[i]]  <- matrix(NA, nrow=GBrS_TPR_curr,ncol=JBrS_list[[i]])
      
      colnames(alpha_mat[[i]]) <- patho_BrS_list[[i]]
      colnames(beta_mat[[i]]) <- patho_BrS_list[[i]]
      
      for (g in 1:GBrS_TPR_curr){
        alpha_mat[[i]][g,] <- unlist(BrS_tpr_prior[[1]][[g]]$alpha)
        beta_mat[[i]][g,]  <- unlist(BrS_tpr_prior[[1]][[g]]$beta)
      }
      
      if (GBrS_TPR_curr>1){
        assign(paste("alphaB", i, sep = "_"), alpha_mat[[i]])      # <---- input BrS TPR prior here.
        assign(paste("betaB", i, sep = "_"),  beta_mat[[i]])    
      }else{
        assign(paste("alphaB", i, sep = "_"), c(alpha_mat[[i]]))   # <---- input BrS TPR prior here.
        assign(paste("betaB", i, sep = "_"),  c(beta_mat[[i]]))    
      }
    }
    names(alpha_mat) <- names(beta_mat)<- names(Mobs$MBS)
    
    
    if (!BrS_TPR_strat){# no stratification.
      # summarize into one name (for all measurements):
      if (length(single_column_MBS)==0){
        # if all slices have >2 columns:
        in_data       <- c(in_data,"Nd","Nu","Jcause","Z_Eti",
                           paste("JBrS",1:length(JBrS_list),sep="_"),
                           paste("MBS",1:length(JBrS_list),sep="_"),
                           paste("templateBS",1:length(JBrS_list),sep="_"),
                           paste("Z_FPR",1:length(JBrS_list),sep="_"),
                           paste("basis_id",1:length(JBrS_list),sep="_")[unlist(has_basis_list)],
                           paste("n_basis",1:length(JBrS_list),sep="_")[unlist(has_basis_list)],
                           paste("non_basis_id",1:length(JBrS_list),sep="_")[unlist(has_non_basis_list)]
                           # paste("alphaB",1:length(JBrS_list),sep="_"),
                           # paste("betaB",1:length(JBrS_list),sep="_")
        )
      } else {
        # if there exist slices with 1 column:
        in_data       <- c(in_data,"Nd","Nu","Jcause","Z_Eti",
                           paste("JBrS",1:length(JBrS_list),sep="_")[-single_column_MBS], # <---- no need to iterate in .bug file for a slice with one column.
                           paste("MBS",1:length(JBrS_list),sep="_"),
                           paste("templateBS",1:length(JBrS_list),sep="_"),
                           paste("Z_FPR",1:length(JBrS_list),sep="_"),
                           paste("basis_id",1:length(JBrS_list),sep="_")[unlist(has_basis_list)],
                           paste("n_basis",1:length(JBrS_list),sep="_")[unlist(has_basis_list)],
                           paste("non_basis_id",1:length(JBrS_list),sep="_")[unlist(has_non_basis_list)]
                           # paste("alphaB",1:length(JBrS_list),sep="_"),
                           # paste("betaB",1:length(JBrS_list),sep="_")
        )
      }
    } else{ # with stratification.
      # summarize into one name (for all measurements):
      if (length(single_column_MBS)==0){
        # if all slices have >2 columns:
        in_data       <- c(in_data,"Nd","Nu","Jcause","Z_Eti",
                           paste("JBrS",1:length(JBrS_list),sep="_"),
                           paste("GBrS_TPR",1:length(JBrS_list),sep="_"),   # <-- added for TPR strata.
                           paste("BrS_TPR_grp",1:length(JBrS_list),sep="_"),# <-- added for TPR strata.
                           paste("MBS",1:length(JBrS_list),sep="_"),
                           paste("templateBS",1:length(JBrS_list),sep="_"),
                           paste("Z_FPR",1:length(JBrS_list),sep="_"),
                           paste("basis_id",1:length(JBrS_list),sep="_")[unlist(has_basis_list)],
                           paste("n_basis",1:length(JBrS_list),sep="_")[unlist(has_basis_list)],
                           paste("non_basis_id",1:length(JBrS_list),sep="_")[unlist(has_non_basis_list)]
                           # paste("alphaB",1:length(JBrS_list),sep="_"),
                           # paste("betaB",1:length(JBrS_list),sep="_")
        )
      } else {
        # if there exist slices with 1 column:
        in_data       <- c(in_data,"Nd","Nu","Jcause","Z_Eti",
                           paste("JBrS",1:length(JBrS_list),sep="_")[-single_column_MBS], # <---- no need to iterate in .bug file for a slice with one column.
                           paste("GBrS_TPR",1:length(JBrS_list),sep="_"),   # <-- added for TPR strata.
                           paste("BrS_TPR_grp",1:length(JBrS_list),sep="_"),# <-- added for TPR strata.
                           paste("MBS",1:length(JBrS_list),sep="_"),
                           paste("templateBS",1:length(JBrS_list),sep="_"),
                           paste("Z_FPR",1:length(JBrS_list),sep="_"),
                           paste("basis_id",1:length(JBrS_list),sep="_")[unlist(has_basis_list)],
                           paste("n_basis",1:length(JBrS_list),sep="_")[unlist(has_basis_list)],
                           paste("non_basis_id",1:length(JBrS_list),sep="_")[unlist(has_non_basis_list)]
                           # paste("alphaB",1:length(JBrS_list),sep="_"),
                           # paste("betaB",1:length(JBrS_list),sep="_")
        )
      }
    }
    
    #
    # hyper-parameters:
    #
    
    # Set BrS measurement priors:
    # hyperparameter for sensitivity (can add for specificity if necessary): 
    for (s in seq_along(Mobs$MBS)){
      #if (likelihood$k_subclass[s] == 1){BrS_tpr_prior <- set_prior_tpr_BrS_NoNest(s,model_options,data_nplcm)}
      #if (likelihood$k_subclass[s] > 1){BrS_tpr_prior <- set_prior_tpr_BrS_NoNest(s,model_options,data_nplcm)}
      
      #assign(paste("alphaB", s, sep = "_"), BrS_tpr_prior[[1]]$alpha)     # <---- input BrS TPR prior here.
      #assign(paste("betaB", s, sep = "_"),  BrS_tpr_prior[[1]]$beta)    
      
      if (likelihood$k_subclass[s]==1){
        out_parameter <- c(out_parameter,paste(c("thetaBS","betaFPR"), s, sep="_"))
        out_parameter <- c(out_parameter,paste(c("taubeta"), s, sep="_")[unlist(has_basis_list)[s]])
      }else{
        cat("==[baker] Yay! Regression with nested subclasses is ready to run.==\n")
        assign(paste("K", s, sep = "_"), likelihood$k_subclass[s])
        in_data       <- unique(c(in_data,paste0("K_",s))) # <---- No. of subclasses for this slice.
        out_parameter <- unique(c(out_parameter,
                                  add_meas_BrS_param_Nest_reg_Slice_jags(s,Mobs,prior,cause_list,FPR_formula)$parameter,
                                  paste(c("ThetaBS","PsiBS"), s, sep="_"))
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
    if (ER_has_basis){
      out_parameter <- c(out_parameter,"betaEti","ER_taubeta") #<-- may not have this parameter in ncs.
    } else{
      out_parameter <- c(out_parameter,"betaEti")
    }
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
        in_data       <- unique(c(in_data,"Nd","Jcause",
                                  paste("JSS",1:length(JSS_list),sep="_"),
                                  paste("MSS",1:length(JSS_list),sep="_"),
                                  paste("templateSS",1:length(JSS_list),sep="_"),
                                  paste("alphaS",1:length(JSS_list),sep="_"),   
                                  paste("betaS",1:length(JSS_list),sep="_")))
      } else{
        in_data       <- unique(c(in_data,"Nd","Nu","Jcause",
                                  paste("JSS",1:length(JSS_list),sep="_")[-single_column_MSS],
                                  paste("MSS",1:length(JSS_list),sep="_"),
                                  paste("templateSS",1:length(JSS_list),sep="_"),
                                  paste("alphaS",1:length(JSS_list),sep="_"),
                                  paste("betaS",1:length(JSS_list),sep="_")))
      }
    }else {
      if (length(single_column_MSS)==0){
        # summarize into one name (for all measurements):
        in_data       <- unique(c(in_data,"Nd","Jcause",
                                  paste("JSS",1:length(JSS_list),sep="_"),
                                  paste("GSS_TPR",1:length(JSS_list),sep="_"),
                                  paste("SS_TPR_grp",1:length(JSS_list),sep="_"),
                                  paste("MSS",1:length(JSS_list),sep="_"),
                                  paste("templateSS",1:length(JSS_list),sep="_"),
                                  paste("alphaS",1:length(JSS_list),sep="_"),   
                                  paste("betaS",1:length(JSS_list),sep="_")))
      } else{
        in_data       <- unique(c(in_data,"Nd","Nu","Jcause",
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
  # 2. set up initialization function:
  #####################################
  if ("BrS" %in% use_measurements & !("SS" %in% use_measurements)){
    in_init       <-   function(){
      res <- list()
      for (s in seq_along(Mobs$MBS)){
        res_curr <- list()
        if (likelihood$k_subclass[s]==1){
          #res_curr[[1]] <- stats::rbeta(JBrS_list[[s]],1,1)
          GBrS_TPR_curr <- eval(parse(text = paste0("GBrS_TPR_",s)))
          if (GBrS_TPR_curr==1){
            res_curr[[1]] <- stats::rbeta(JBrS_list[[s]],1,1)
          } else{
            res_curr[[1]] <- matrix(stats::rbeta(GBrS_TPR_curr*JBrS_list[[s]],1,1),
                                    nrow=GBrS_TPR_curr,ncol=JBrS_list[[s]])
            if (JBrS_list[[s]]==1){
              res_curr[[1]] <- c(res_curr[[1]])
            }
          }
          res_curr[[2]] <- matrix(0,nrow=ncol(Z_FPR_list[[s]]),ncol=JBrS_list[[s]])
          names(res_curr) <- paste(c("thetaBS","betaFPR"),s,sep="_")
          
          res <- c(res,res_curr)
        }
        if (likelihood$k_subclass[s] > 1){
          GBrS_TPR_curr <- eval(parse(text = paste0("GBrS_TPR_",s)))
          if (GBrS_TPR_curr>1){
            stop("==[baker] Cannot deal with multiple TPRs for nested subclasses for regression. Please ask maintainer.==\n")
          }
          K_curr <- likelihood$k_subclass[s]
          res_curr[[1]] <- matrix(0.5,nrow=JBrS_list[[s]],ncol=K_curr)
          res_curr[[2]] <- cbind(matrix(0,nrow=ncol(Z_FPR_list[[s]]),ncol=K_curr-1),rep(NA,ncol(Z_FPR_list[[s]])))
          res_curr[[3]] <- sample(1:K_curr,Nd+Nu,replace=TRUE)
          res_curr[[4]] <- cbind(matrix(0,nrow=ncol(Z_FPR_list[[s]]),ncol=K_curr-1),rep(NA,ncol(Z_FPR_list[[s]])))
          res_curr[[5]] <- matrix(0.5,nrow=JBrS_list[[s]],ncol=K_curr)
          
          names(res_curr) <- paste(c("PsiBS","betaFPR","Z","case_betaFPR","ThetaBS"),s,sep="_")
          res <- c(res,res_curr)
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
          if (JSS_list[[i]]==1){
            tmp_thetaSS[[i]] <- c(tmp_thetaSS[[i]])
          }
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
          #res_curr[[1]] <- stats::rbeta(JBrS_list[[s]],1,1)
          GBrS_TPR_curr <- eval(parse(text = paste0("GBrS_TPR_",s)))
          if (GBrS_TPR_curr==1){
            res_curr[[1]] <- stats::rbeta(JBrS_list[[s]],1,1)
          } else{
            res_curr[[1]] <- matrix(stats::rbeta(GBrS_TPR_curr*JBrS_list[[s]],1,1),
                                    nrow=GBrS_TPR_curr,ncol=JBrS_list[[s]])
            if (JBrS_list[[s]]==1){
              res_curr[[1]] <- c(res_curr[[1]])
            }
          }
          res_curr[[2]] <- matrix(0,nrow=ncol(Z_FPR_list[[s]]),ncol=JBrS_list[[s]])
          names(res_curr) <- paste(c("thetaBS","betaFPR"),s,sep="_")
          res <- c(res,res_curr)
        }
        if (likelihood$k_subclass[s] > 1){
          GBrS_TPR_curr <- eval(parse(text = paste0("GBrS_TPR_",s)))
          if (GBrS_TPR_curr>1){
            stop("==[baker] Cannot deal with multiple TPRs for nested subclasses for regression. Please ask maintainer.==\n")
          }
          K_curr <- likelihood$k_subclass[s]
          res_curr[[1]] <- matrix(0.5,nrow=JBrS_list[[s]],ncol=K_curr)
          res_curr[[2]] <- cbind(matrix(0,nrow=ncol(Z_FPR_list[[s]]),ncol=K_curr-1),rep(NA,ncol(Z_FPR_list[[s]])))
          res_curr[[3]] <- sample(1:K_curr,Nd+Nu,replace=TRUE)
          res_curr[[4]] <- cbind(matrix(0,nrow=ncol(Z_FPR_list[[s]]),ncol=K_curr-1),rep(NA,ncol(Z_FPR_list[[s]])))
          res_curr[[5]] <- matrix(0.5,nrow=JBrS_list[[s]],ncol=K_curr)
          
          names(res_curr) <- paste(c("PsiBS","betaFPR","Z","case_betaFPR","ThetaBS"),s,sep="_")
          res <- c(res,res_curr)
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
          if (JSS_list[[i]]==1){
            tmp_thetaSS[[i]] <- c(tmp_thetaSS[[i]])
          }
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
  
  ######################################
  # 3. fit model :
  ######################################
  # special operations:
  # get individual prediction:
  if (!is.null(mcmc_options$individual.pred) && mcmc_options$individual.pred){out_parameter <- c(out_parameter,"Icat")}
  # get posterior predictive distribution of BrS measurments:
  if (!is.null(mcmc_options$ppd) && mcmc_options$ppd){out_parameter <- c(out_parameter,paste("MBS.new",seq_along(Mobs$MBS),sep = "_"))}
  # get pEti: (for regression models; besides the coefficients of etiology regressions are always recorded)
  if (!is.null(mcmc_options$get.pEti) && mcmc_options$get.pEti){out_parameter <- unique(c(out_parameter,"pEti"))}
  
  # write the .bug files into mcmc_options$bugsmodel.dir; 
  # could later set it equal to result.folder.
  use_jags <- (!is.null(mcmc_options$use_jags) && mcmc_options$use_jags)
  model_func         <- write_model_Reg_Nest(data_nplcm$Mobs,
                                             model_options$prior,
                                             model_options$likelihood$cause_list,
                                             model_options$likelihood$Eti_formula,
                                             model_options$likelihood$FPR_formula,
                                             model_options$use_measurements,
                                             mcmc_options$ppd,
                                             use_jags)
  model_bugfile_name <- "model_Reg_Nest.bug"
  
  filename <- file.path(mcmc_options$bugsmodel.dir, model_bugfile_name)
  writeLines(model_func, filename)
  
  # if (length(prior$Eti_prior)>1){
  #   stop("== [baker] Regression model used. Please change `model_options$prior$Eti_prior` 
  #        into a single, positive real number representing
  #        the stanadrd deviation of beta coefficients in etiology regression! ==")
  # }
  sd_betaEti_basis    <- prior$Eti_prior[1]
  sd_betaEti_nonbasis <- prior$Eti_prior[2]   #<--place to adjust (currently identical for basis vs non-basis. check "insert_bugfile_chunk_reg_etiology")
  sd_betaFPR_basis   <- prior$FPR_coef_prior[1]
  sd_betaFPR_nonbasis <- prior$FPR_coef_prior[2]
  ER_alpha    <- prior$Eti_hyper_pflex[1]
  ER_beta     <- prior$Eti_hyper_pflex[2]
  
  in_data <- unique(c(in_data,
                      "ER_basis_id"[ER_has_basis],
                      "ER_n_basis"[ER_has_basis],
                      "ER_non_basis_id"[ER_has_non_basis],
                      "sd_betaEti_basis"[ER_has_basis],
                      "ER_alpha"[ER_has_basis],
                      "ER_beta"[ER_has_basis],
                      "sd_betaEti_nonbasis"[ER_has_non_basis],
                      "sd_betaFPR_basis"[any(unlist(has_basis_list))],
                      "sd_betaFPR_nonbasis"[any(unlist(has_non_basis_list))]
                      
                      
  ))
  
  
  for (s in seq_along(Mobs$MBS)){
    assign(paste("d_FPR",s,sep="_"),ncol(Z_FPR_list[[s]]))
      in_data <- c(in_data,paste("d_FPR",s,sep="_"))
  }
  
  # # # uncomment below if using dmnorm for betaEti and betaFPR (also need to edit plug-and-play.R):
  # #
  # # if (Jcause > 2){ # use dmnorm to speed up JAGS calculations, so we need mean 
  # #   # and covariance vector.
  # #   zero_Jcause_1 <- rep(0,Jcause-1)
  # #   I_Jcause_1    <- diag(1,Jcause-1)
  # #   in_data <- c(in_data,"zero_Jcause_1","I_Jcause_1")
  # # } 
  # 
  # # for (s in seq_along(JBrS_list)){
  # #   if (JBrS_list[[s]]>1){
  # #       assign(paste("zero_JBrS", s, sep = "_"), rep(0,JBrS_list[[s]]))    
  # #       assign(paste("I_JBrS", s, sep = "_"), diag(1,JBrS_list[[s]]))    
  # #       in_data <- c(in_data,
  # #                    paste("zero_JBrS", s, sep = "_"),
  # #                    paste("I_JBrS", s, sep = "_"))
  # #   }
  # # }
  
  #
  # run the model:
  #
  if (!use_jags){stop("==[baker] Not supporting WinBUGS any more; Use JAGS for now.==\n")}
  here <- environment()
  if (use_jags){
    ##JAGS
    in_data.list <- lapply(as.list(in_data), get, envir= here)
    names(in_data.list) <- in_data
    #lapply(names(in_data.list), dump, append = TRUE, envir = here,
    #       file = file.path(mcmc_options$result.folder,"jagsdata.txt"))
    curr_data_txt_file <- file.path(mcmc_options$result.folder,"jagsdata.txt")
    if(file.exists(curr_data_txt_file)){file.remove(curr_data_txt_file)}
    dump(names(in_data.list), append = FALSE, envir = here,
         file = curr_data_txt_file)
    ## fix dimension problem.... convert say .Dmi=7:6 to c(7,6) (an issue for templateBS_1):
    bad_jagsdata_txt <- readLines(curr_data_txt_file)
    good_jagsdata_txt <- gsub( "([0-9]+):([0-9]+)", "c(\\1,\\2)", bad_jagsdata_txt,fixed = FALSE)
    writeLines(good_jagsdata_txt, curr_data_txt_file)
    
    gs <- jags2_baker(data   = curr_data_txt_file,
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















