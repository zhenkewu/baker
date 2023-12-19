if(getRversion() >= "2.15.1") utils::globalVariables(c("set_prior_tpr","set_prior_eti"))


#' Fit nested partially-latent class models (highest-level wrapper function)
#'
#' Uses `JAGS` (OSX or Windows) operating system for Bayesian posterior inference
#' (see `README` file for an instruction to install `JAGS`). If running `JAGS` on windows, 
#' please go to control panel to add the directory to `JAGS` into ENVIRONMENTAL VARIABLE.
#' 
#' @param data_nplcm Cases are on top of controls in the rows of diagnostic
#' test results and the covariate matrix. This is assumed by `baker` to automatically 
#' write model files (`.bug`). 
#' \itemize{
#' \item  `Mobs` A list of measurements of distinct qualities (Bronze-, Silver, and Gold-Standard:
#' `MBS`,`MSS`,`MGS`). The elements of the list
#' should include `MBS`, `MSS`, and `MGS`. If any of the component
#' is not available, please specify it as, e.g., `MGS=NULL`
#' (effectively deleting `MGS` from `Mobs`).
#' \itemize{
#' \item `MBS` a list of data frame of bronze-standard (BrS) measurements.
#' For each data frame (referred to as a 'slice'), 
#' rows are subjects, columns are causative agents (e.g., pathogen species). 
#' We use `list` here to accommodate the possibility of multiple sets of BrS data.
#' They have imperfect sensitivity/specificity (e.g. nasopharyngeal polymerase chain
#' reaction - NPPCR). 
#' \item `MSS` a list of data frame of silver-standard (SS) measurements. 
#' Rows are subjects, columns are causative agents measured in specimen (e.g. blood culture).
#' These measurements have perfect specificity but imperfect sensitivity.
#' \item `MGS` a list of data frame of gold-standard (GS) measurements. 
#' Rows are subject, columns are measured causative agents
#' These measurements have perfect sensitivity and specificity.
#' }
#'
#' \item `Y` Vector of disease status: `1` for case, `0` for control.
#' \item `X` Covariate matrix. A subset of columns are primary covariates in cause-specific-
#' case-fraction (CSCF) functions and hence must be available for cases, and another subset
#' are covariates that are available in the cases and the controls. 
#' The two sets of covariates may be identical, overlapping or completely different.
#' In general, this is not the design matrix for regression models,
#' because for enrollment date in a study which may have non-linear effect,
#' basis expansion is often needed for approximation.
#' }
#' 
#' @param model_options A list of model options: likelihood and prior. 
#' \describe{
#' \item{`use_measurements`}{
#'   A vector of characters strings; can be one or more from `"BrS"`, `"SS"`, `"GS"`.
#' }
#' \item{`likelihood`}{
#'     \describe{
#'          \item{cause_list}{ The vector of causes (NB: specify);}
#'          \item{k_subclass}{ The number of nested subclasses in each 
#'          disease class (one of case classes or the control class; the same `k_subclass`
#'          is assumed for each class) and each slice of BrS measurements. 
#'          `1` for conditional independence; larger than `1` for conditional dependence. 
#'          It is only available for BrS measurements. It is a vector of length equal to 
#'          the number of slices of BrS measurements;}
#'          \item{Eti_formula}{ Formula for etiology regressions. You can use  
#'          [s_date_Eti()] to specify the design matrix for `R` format enrollment date; 
#'          it will produce natural cubic spline basis. Specify `~ 1` if no regression is intended.}
#'          \item{FPR_formula}{formula for false positive rates (FPR) regressions; see [formula()]. 
#'          You can use [s_date_FPR()] to specify part of the design matrix for `R` 
#'          format enrollment date; it will produce penalized-spline basis (based on B-splines). 
#'          Specify `~ 1` if no regression is intended. (NB: If `effect="fixed"`, [dm_Rdate_FPR()]
#'          will just specify a design matrix with appropriately standardized dates.)}
#'     }
#' }
#' 
#' \item{`prior`}{
#'     \describe{
#'          \item{Eti_prior}{Description of etiology prior (e.g., `overall_uniform` - 
#'          all hyperparameters are `1`; or `0_1` - all hyperparameters are `0.1`);}
#'          \item{TPR_prior}{Description of priors for the measurements 
#'          (e.g., informative vs non-informative). Its length should be the 
#'          same as `use_measurements` above. Please see examples for how to specify.
#'          The package can also handle multiple slices of BrS, SS data, so separate
#'          specification of the TPR priors are needed. 
#'        }
#'     }
#' }
#' }
#' @param mcmc_options A list of Markov chain Monte Carlo (MCMC) options.
#' \itemize{
#' \item `debugstatus` Logical - whether to pause WinBUGS after it finishes
#' model fitting; (NB: is this obsolete? Test.)
#' \item `n.chains` Number of MCMC chains;
#' \item `n.burnin` Number of burn-in iterations;
#' \item `n.thin` To keep every other `n.thin` samples after burn-in period;
#' \item `individual.pred` `TRUE` to perform individual prediction (`Icat`
#' variables in the `.bug` file); `FALSE` otherwise;
#' \item `ppd` `TRUE` to simulate new data (`XXX.new`
#' variables in the `.bug` file) from the posterior predictive distribution (ppd); 
#' `FALSE` otherwise;
#' \item `get.pEti` `TRUE` for getting posterior samples of individual etiologic fractions; 
#' `FALSE` otherwise. For non-regression, or regression models with all discrete predictors, 
#' by default this is `TRUE`, so no need to specify this entry. It is only relevant for regression models
#' with non-discrete covariates. Because individuals have distinct CSCFs at their specific covariate values, 
#' it's easier to just store the posterior samples of the regression coefficients and reconstruct the pies afterwards,
#' rather than storing them through `JAGS`. 
#' \item `result.folder` Path to folder storing the results;
#' \item `bugsmodel.dir` Path to `.bug` model files;
#' \item `jags.dir` Path to where JAGS is installed; if `NULL`, this will be set 
#' to `jags.dir=""`.
#' }
#' @return A `JAGS` output result, fitted by function [R2jags::jags2()] from `R2jags`. 
#' It is an object of class `nplcm` and `bugs`.
#' Current implemented models follow the hierarchy below:
#' \itemize{
#' \item no regression:  Fitted by at low level by [nplcm_fit_NoReg]
#'        
#' \item regression: 
#'  Given disease class (control or a class of cases with the same subset of causative agents):
#'       \itemize{
#'          \item local independence model for BrS measures:
#'          Fitted at lower level by 
#'          \itemize{
#'          \item [nplcm_fit_Reg_NoNest] deals with the setting with two sets of covariates,
#'          one for CSCF regression and the other for FPR regression. The two sets of
#'          covariates may be identical, overlapping or non-overlapping.
#'          This function is called when there exists one or more than one discrete covariate among
#'          the union of the two covariate sets. The method implemented by this function
#'          directly lets FPR depend upon covariates. 
#'          This is different from Wu and Chen (2021), which let the subclass 
#'          weights depend upon covariates. We implemented this function for methods comparison.
#'           \item [nplcm_fit_Reg_discrete_predictor_NoNest] deals with the setting
#'           with all discrete covariates for FPRs and CSCFs. The strata defined by the two sets of 
#'           covariates need not be identical, e.g., as a result of distinct sets of covariates. Again,
#'           this is directly to let FPR be stratified by covariates, hence different from Wu and Chen (2020+)
#'           We implemented this function for methods comparison.
#'           }
#'          \item local dependence model for BrS measures: 
#'          Fitted at lower level by [nplcm_fit_Reg_Nest]: This is the method introduced in 
#'          Wu and Chen (2021): CSCF regression + case/control subclass weight regression.
#'          It does not provide a specialized function for the setting with all discrete covariates.
#'        }
#' }
#' 
#' @examples 
#' 
#' \donttest{
#' data(data_nplcm_noreg)
#' cause_list <- LETTERS[1:6]
#' J.BrS      <- 6
#' model_options_no_reg <- list(
#'   likelihood   = list(
#'     cause_list = cause_list,
#'     k_subclass = 2,
#'     Eti_formula = ~-1, # no covariate for the etiology regression
#'     FPR_formula = list(
#'       MBS1 =   ~-1)    # no covariate for the subclass weight regression
#'   ),
#'   use_measurements = c("BrS"), 
#'   # use bronze-standard data only for model estimation.
#'   prior= list(
#'     Eti_prior = overall_uniform(1,cause_list), 
#'     # Dirichlet(1,...,1) prior for the etiology.
#'     TPR_prior  = list(BrS = list(
#'       info  = "informative", # informative prior for TPRs
#'       input = "match_range", 
#'       # specify the informative prior for TPRs by specifying a plausible range.
#'       val = list(MBS1 = list(up =  list(rep(0.99,J.BrS)), 
#'                              # upper ranges: matched to 97.5% quantile of a Beta prior
#'                              low = list(rep(0.55,J.BrS))))
#'       # lower ranges: matched to 2.5% quantile of a Beta prior
#'     )
#'     )
#'   )
#' )     
#' 
#' 
#' set.seed(1)
#' # include stratification information in file name:
#' thedir    <- paste0(tempdir(),"_no_reg")
#' 
#' # create folders to store the model results 
#' dir.create(thedir, showWarnings = FALSE)
#' result_folder_no_reg <- file.path(thedir,paste("results",collapse="_"))
#' thedir <- result_folder_no_reg
#' dir.create(thedir, showWarnings = FALSE)
#' 
#' # options for MCMC chains:
#' mcmc_options_no_reg <- list(
#'   debugstatus = TRUE,
#'   n.chains = 1,
#'   n.itermcmc = as.integer(200), 
#'   n.burnin = as.integer(100), 
#'   n.thin = 1,
#'   individual.pred = TRUE, # <- must set to TRUE! <------- NOTE! 
#'   ppd = FALSE,
#'   result.folder = thedir,
#'   bugsmodel.dir = thedir
#' )
#' 
#' BrS_object_1 <- make_meas_object(patho = LETTERS[1:6], 
#'                                  specimen = "MBS", test = "1", 
#'                                  quality = "BrS", cause_list = cause_list)
#' clean_options <- list(BrS_objects = make_list(BrS_object_1))
#' # place the nplcm data and cleaning options into the results folder
#' dput(data_nplcm_noreg,file.path(thedir,"data_nplcm.txt")) 
#' dput(clean_options, file.path(thedir, "data_clean_options.txt"))
#' 
#' rjags::load.module("glm")
#' 
#' nplcm_noreg <- nplcm(data_nplcm_noreg,model_options_no_reg,mcmc_options_no_reg)
#' 
#' 
#' }
#' 
#' 
#' @export
nplcm <- function(data_nplcm,model_options,mcmc_options){
  mcmc_options$use_jags <- TRUE  
  Mobs <- data_nplcm$Mobs
  Y    <- data_nplcm$Y
  if (length(rle(Y)[["values"]])>2 | Y[1]!=1) {
    stop("==[baker] 'data_nplcm$Y' must have cases on top of controls. 
  Use 'baker::subset_data_nplcm_by_index()' to shuffle the rows. Then retry.==\n")}
  X    <- data_nplcm$X
  
  parsed_model <- assign_model(model_options,data_nplcm)
  
  do_reg_vec <- unlist(parsed_model$regression[grep("^do_reg_",names(parsed_model$regression))])
  do_discrete_vec <- unlist(parsed_model$regression[grep("^is_discrete_predictor",names(parsed_model$regression))])
  do_discrete_and_reg <- do_reg_vec&do_discrete_vec
  do_reg <- any(do_reg_vec)
  do_discrete <- sum(do_discrete_and_reg)==sum(do_reg_vec) 
  # only call nplcm_fit_Reg_discrete_predictor_NoNest
  # when ALL regressions use discrete predictors!
  
  do_nested <- parsed_model$nested
  
  fitted_type <- NA
  if(!do_reg){
    res <- nplcm_fit_NoReg(data_nplcm,model_options,mcmc_options)
    fitted_type <- "no_reg"
  } 
  if (do_reg & !any(do_nested)){
    if (do_discrete){
      fitted_type   <- "reg_nonest_strat"
      # when every regression is a regression upon discrete variables;
      # when it is not a regression, the fitting function treats it as a single stratum
      # when specifying the model in the .bug file (in the assign_model function, ~1 
      # is considered not a regression, i.e., covariate independent)
      res <- nplcm_fit_Reg_discrete_predictor_NoNest(data_nplcm,model_options,mcmc_options)
    } else{
      fitted_type <- "reg_nonest"
      # a mix of discrete and continuous regressions:
      res <- nplcm_fit_Reg_NoNest(data_nplcm,model_options,mcmc_options)
    }
  }
  if (do_reg & any(do_nested)){
    if (do_discrete) {fitted_type <- "reg_nest_strat"}else{fitted_type <- "reg_nest"}
    res <- nplcm_fit_Reg_Nest(data_nplcm,model_options,mcmc_options)
  }
  res$DIR_NPLCM <- mcmc_options$result.folder  # <--- add info about where results are stored.
  res$fitted_type <- fitted_type
  class(res) <- c("nplcm","bugs")
  return(res)
}



#' Interpret the specified model structure
#'
#' `assign_model` translates options specified by a user (e.g., in 
#' `model_options`) into information that can be understood by `baker`.
#' 
#' @details `assign_model` will be modified to check if data are conformable
#' to specified model.
#' 
#' @param data_nplcm Data. See [nplcm()] function for data structure.
#' @param model_options See [nplcm()] function.
#' @param silent Default is `TRUE` for no messages; `FALSE` otherwise.
#' @return A list of model specifications:
#' \itemize{
#'    \item `num_slice` A vector counting the No. of measurement slices for each
#'    level of measurement quality (e.g., MBS, MSS, MGS representing
#'    Bronze-Standard Measurements - case-control, 
#'    Silver-Standard Measurements and Gold-Standard
#'    Measurements - case-only);
#'    \item `nested` Local dependence specification for modeling bronze-standard
#'    data. `TRUE` for nested models (conditional dependence given disease class); 
#'    `FALSE` for non-nested models (conditional independence given disease class). 
#'    One for each BrS slice.
#'    \item `regression`
#'        \itemize{
#'            \item `do_reg_Eti` `TRUE` for doing etiology regression.
#'            It means let the etiology fractions vary with explanatory variables. 
#'            `FALSE` otherwise;
#'            \item `do_reg_FPR` A vector whose names represent the slices
#'            of bronze-standard data. For each slice of BrS measurements, 
#'            `TRUE` does false positive rate regression. It means the false
#'            positive rates, estimatable from controls, can vary with 
#'            covariates; `FALSE` otherwise.
#'            \item `is_discrete_predictor` A list of names "Eti", and 
#'            the names for every slice of bronze-standard data. `TRUE`
#'            if all predictors are discrete; `FALSE` otherwise.
#'        }
#' }
#' 
#' @examples
#' cause_list <- c(LETTERS[1:6]) 
#' J.BrS <- 6
#' model_options_no_reg <- list(
#' likelihood   = list(
#'   cause_list = cause_list,
#'   k_subclass = 2,
#'   Eti_formula = ~-1, 
#'   # no covariate for the etiology regression
#'   FPR_formula = list(
#'     MBS1 =   ~-1)    
#'     # no covariate for the subclass weight regression
#' ),
#' use_measurements = c("BrS"), 
#' # use bronze-standard data only for model estimation.
#' prior= list(
#'   Eti_prior = overall_uniform(1,cause_list), 
#'   # Dirichlet(1,...,1) prior for the etiology.
#'   TPR_prior  = list(BrS = list(
#'     info  = "informative", # informative prior for TPRs
#'     input = "match_range", 
#'     # specify the informative prior for TPRs by specifying a plausible range.
#'     val = list(MBS1 = list(up =  list(rep(0.99,J.BrS)), 
#'     # upper ranges: matched to 97.5% quantile of a Beta prior
#'                            low = list(rep(0.55,J.BrS))))
#'                            # lower ranges: matched to 2.5% quantile of a Beta prior
#'   )
#'   )
#' )
#' )     
#' data("data_nplcm_noreg")
#' 
#' assign_model(model_options_no_reg,data_nplcm_noreg)
#' 
#' @family specification checking functions
#' @export
assign_model <- function(model_options,data_nplcm, silent=TRUE){
  # load options:
  likelihood       <- model_options$likelihood
  use_measurements <- model_options$use_measurements
  prior            <- model_options$prior
  
  # load data: 
  Mobs <- data_nplcm$Mobs
  Y    <- data_nplcm$Y
  X    <- data_nplcm$X
  
  nested       <- likelihood$k_subclass > 1
  
  # test the match between actual data and model_options:
  use_data_sources   <- c("MBS","MSS","MGS")[lookup_quality(use_measurements)]
  input_data_sources <-  names(Mobs)
  if (!all(use_data_sources%in%input_data_sources)){
    stop("==[baker] Please supply actual datasets as specified by 'use_measurements' in 'model_options'.==\n")
  }
  
  # get the length of each measurement quality:
  num_slice <- rep(0,3)
  names(num_slice) <- c("MBS","MSS","MGS")
  for (i in seq_along(use_data_sources)){
    num_slice[use_data_sources[i]] <- length(Mobs[[use_data_sources[i]]])
  }
  
  # specify regression for FPR: (only available for bronze-standard data. Silver-standard data automatically have FPR==0.)
  do_reg_FPR <- is_discrete_FPR <- rep(NA,length(likelihood$FPR_formula)) #  <---- a regression for each measurement slice?
  names(do_reg_FPR) <- names(is_discrete_FPR) <- names(likelihood$FPR_formula)
  for (i in seq_along(Mobs$MBS)){
    ind_tmp <-which(names(likelihood$FPR_formula) == names(Mobs$MBS)[i])
    form_tmp <- stats::as.formula(likelihood$FPR_formula[[ind_tmp]])
    if (!length(ind_tmp)) { # don't do regression if no regression formula is found:
      do_reg_FPR[i] <- FALSE
    } else{ # do regression if there is matched regression formula:
      do_reg_FPR[i] <-
        parse_nplcm_reg(form_tmp,data_nplcm,silent=silent)
    }
    
    is_discrete_FPR[i] <- FALSE
    if (!is.null(X)){
      is_discrete_FPR[i] <- (!is_intercept_only(form_tmp) & 
                               !is_intercept_only(form_tmp)& 
                               is_discrete(data.frame(X,Y), form_tmp))
    }
  }
  #
  # specify regression for TPR: (every measurement slice has it.)
  #
  # do_reg_TPR <- list() #  <---- a regression for each measurement slice?
  # for (i in seq_along(Mobs$MBS)){
  #   ind_tmp <-
  #     which(names(likelihood$TPR_formula) == names(Mobs$MBS)[i])
  #   if (!length(ind_tmp)) { # don't do regression if no regression formula is found:
  #     do_reg_TPR[[i]] <- FALSE
  #   } else{ # do regression if there is matched regression formula:
  #     do_reg_TPR[[i]] <-
  #       parse_nplcm_reg(stats::as.formula(likelihood$TPR_formula[[ind_tmp]]),data_nplcm,silent=silent)
  #   }
  # }
  # 
  # names(do_reg_TPR) <- names(Mobs$MBS)
  # 
  # # if using silver-standard data:
  # if ("MSS"%in% use_data_sources){
  #   for (i in length(Mobs$MBS)+seq_along(Mobs$MSS)){
  #     ind_tmp <-
  #       which(names(likelihood$TPR_formula) == names(Mobs$MSS)[i])
  #     if (!length(ind_tmp)) { # don't do regression if no regression formula is found:
  #       do_reg_TPR[[i]] <- FALSE
  #     } else{ # do regression if there is matched regression formula:
  #       do_reg_TPR[[i]] <-
  #         parse_nplcm_reg(stats::as.formula(likelihood$TPR_formula[[ind_tmp]]),data_nplcm,silent=silent)
  #     }
  #   }
  #   names(do_reg_TPR) <- c(names(Mobs$MBS),names(Mobs$MSS))
  # }
  # specify regression for etiology:
  form_tmp   <- stats::as.formula(likelihood$Eti_formula)
  do_reg_Eti <- parse_nplcm_reg(form_tmp,data_nplcm,silent=silent)
  
  is_discrete_Eti <- FALSE
  if (!is.null(X)){ # <--- potential problem if a user input more data than needed. need fixing.
    is_discrete_Eti <- (!stats::is.empty.model(form_tmp) & 
                          !is_intercept_only(form_tmp) & 
                          is_discrete(data.frame(X,Y)[Y==1,,drop=FALSE], form_tmp))
  }
  
  is_discrete_predictor <- list(is_discrete_Eti, is_discrete_FPR)
  names(is_discrete_predictor)[1] <- "Eti"
  names(is_discrete_predictor)[2] <- "FPR"
  regression <- make_list(do_reg_Eti, do_reg_FPR,is_discrete_predictor)#, do_reg_TPR)
  
  # check BrS group:
  BrS_grp   <- FALSE
  prior_BrS <- model_options$prior$TPR_prior$BrS
  GBrS_TPR  <- length(unique(prior_BrS$grp))
  grp_spec  <- (!is.null(prior_BrS$grp) && GBrS_TPR >1 )
  if (grp_spec) {
    for (s in seq_along(prior_BrS$val)){
      if (prior_BrS$input=="match_range" && 
          (length(prior_BrS$val[[s]]$up)!=GBrS_TPR | 
           length(prior_BrS$val[[s]]$low)!=GBrS_TPR) ){
        stop(paste0("==[baker] ",names(prior_BrS$val)[s])," needs ", GBrS_TPR,
             " sets of sensitivity ranges.==\n")
      }
    }
    BrS_grp <- TRUE
  }
  
  
  # check SS group:
  SS_grp   <- FALSE
  prior_SS <- model_options$prior$TPR_prior$SS
  grp_spec <- (!is.null(prior_SS$grp) && length(unique(prior_SS$grp)) >1 )
  if (grp_spec) {SS_grp <- TRUE}
  
  ## <-------- the following are more strict grp specifications (may cause error when running old folders):
  #   val_spec <- (num_slice["MSS"]>0 && any(lapply(prior_SS$val,length)>1))
  #   
  #   if (grp_spec && val_spec){SS_grp <- TRUE}
  #   if (grp_spec && !val_spec){stop("==Specified TPR group in 'grp' of 'model_options$prior$TPR_prior$SS',
  #                                   but either there is no SS data or the length of 'val' does not match the no. of TPR groups. ==")}
  #   if (!grp_spec && val_spec){stop("==No 'grp' specified in 'model_options$prior$TPR_prior$SS',
  #                                   but we have >1 sets of TPRs. ==")}
  
  # return results:
  make_list(num_slice, nested, regression,BrS_grp,SS_grp)
}




#' Fit nested partially-latent class model (low-level)
#'
#' This function prepares data, specifies hyperparameters in priors 
#' (true positive rates and etiology fractions), initializes the posterior
#' sampling chain, writes the model file (for JAGS or WinBUGS with slight differences
#' in syntax), and fits the model. Features:
#' \itemize{
#' \item no regression;
#' \item no nested subclasses
#' }
#'
#' @inheritParams nplcm
#' @return BUGS fit results.
#' 
#' @seealso [write_model_NoReg] for constructing `.bug` model file; This function
#' then put it in the folder `mcmc_options$bugsmodel.dir`.
#' 
#' @family model fitting functions 
#' 
#'    
nplcm_fit_NoReg<-
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
    
    # read in options:
    likelihood       <- model_options$likelihood
    use_measurements <- model_options$use_measurements
    prior            <- model_options$prior
    
    #####################################################################
    # 1. prepare data (including hyper-parameters):
    #####################################################################
    
    # sample sizes:
    Nd      <- sum(Y==1)
    Nu      <- sum(Y==0)
    
    # lengths of vectors:
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
          warning(paste0("==[baker] Bronze-standard slice ", names(data_nplcm$Mobs$MBS)[s], " has no measurements informative of the causes specified in 'cause_list', except 'NoA'! 
                         Please check if you need this measurement slice columns correspond to causes other than 'NoA'.=="))  
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
      
      for(i in seq_along(JBrS_list)){
        assign(paste("JBrS", i, sep = "_"), JBrS_list[[i]])    
        xx <- as.matrix_or_vec(MBS_list[[i]]);  attr(xx,"dimnames") <- NULL # remove dimnames.
        assign(paste("MBS", i, sep = "_"), xx)
        assign(paste("templateBS", i, sep = "_"), as.matrix_or_vec(template_BrS_list[[i]]))   
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
        if (likelihood$k_subclass[i] > 1) {BrS_tpr_prior <- set_prior_tpr_BrS_NoNest(i,model_options,data_nplcm)}
        
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
      
      if (!BrS_TPR_strat){
        # summarize into one name (for all measurements):
        if (length(single_column_MBS)==0){
          # if all slices have >2 columns:
          in_data       <- c(in_data,"Nd","Nu","Jcause","alphaEti",
                             paste("JBrS",1:length(JBrS_list),sep="_"),
                             paste("MBS",1:length(JBrS_list),sep="_"),
                             paste("templateBS",1:length(JBrS_list),sep="_")
                             # paste("alphaB",1:length(JBrS_list),sep="_"),
                             # paste("betaB",1:length(JBrS_list),sep="_")
          )
        } else {
          # if there exist slices with 1 column:
          in_data       <- c(in_data,"Nd","Nu","Jcause","alphaEti",
                             paste("JBrS",1:length(JBrS_list),sep="_")[-single_column_MBS], # <---- no need to iterate in .bug file for a slice with one column.
                             paste("MBS",1:length(JBrS_list),sep="_"),
                             paste("templateBS",1:length(JBrS_list),sep="_")
                             # paste("alphaB",1:length(JBrS_list),sep="_"),
                             # paste("betaB",1:length(JBrS_list),sep="_")
          )
        }
      } else{
        # summarize into one name (for all measurements):
        if (length(single_column_MBS)==0){
          # if all slices have >2 columns:
          in_data       <- c(in_data,"Nd","Nu","Jcause","alphaEti",
                             paste("JBrS",1:length(JBrS_list),sep="_"),
                             paste("GBrS_TPR",1:length(JBrS_list),sep="_"),   # <-- added for TPR strata.
                             paste("BrS_TPR_grp",1:length(JBrS_list),sep="_"),# <-- added for TPR strata.
                             paste("MBS",1:length(JBrS_list),sep="_"),
                             paste("templateBS",1:length(JBrS_list),sep="_")
                             # paste("alphaB",1:length(JBrS_list),sep="_"),
                             # paste("betaB",1:length(JBrS_list),sep="_")
          )
        } else {
          # if there exist slices with 1 column:
          in_data       <- c(in_data,"Nd","Nu","Jcause","alphaEti",
                             paste("JBrS",1:length(JBrS_list),sep="_")[-single_column_MBS], # <---- no need to iterate in .bug file for a slice with one column.
                             paste("GBrS_TPR",1:length(JBrS_list),sep="_"),   # <-- added for TPR strata.
                             paste("BrS_TPR_grp",1:length(JBrS_list),sep="_"),# <-- added for TPR strata.
                             paste("MBS",1:length(JBrS_list),sep="_"),
                             paste("templateBS",1:length(JBrS_list),sep="_")
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
        }else{  # <--- TPR stratification for BrS data not implemented for K>1.
          assign(paste("K", s, sep = "_"), likelihood$k_subclass[s])
          in_data       <- unique(c(in_data,paste0("K_",s))) # <---- No. of subclasses for this slice.
          out_parameter <- unique(c(out_parameter,
                                    paste(c("ThetaBS","PsiBS","Lambda","Eta","alphadp0","alphadp0_case"),s,sep="_")
          ))
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
          in_data       <- unique(c(in_data,"Nd","Jcause","alphaEti",
                                    paste("JSS",1:length(JSS_list),sep="_"),
                                    paste("MSS",1:length(JSS_list),sep="_"),
                                    paste("templateSS",1:length(JSS_list),sep="_"),
                                    paste("alphaS",1:length(JSS_list),sep="_"),   
                                    paste("betaS",1:length(JSS_list),sep="_")))
        } else{
          in_data       <- unique(c(in_data,"Nd","Nu","Jcause","alphaEti",
                                    paste("JSS",1:length(JSS_list),sep="_")[-single_column_MSS],
                                    paste("MSS",1:length(JSS_list),sep="_"),
                                    paste("templateSS",1:length(JSS_list),sep="_"),
                                    paste("alphaS",1:length(JSS_list),sep="_"),
                                    paste("betaS",1:length(JSS_list),sep="_")))
        }
      }else {
        if (length(single_column_MSS)==0){
          # summarize into one name (for all measurements):
          in_data       <- unique(c(in_data,"Nd","Jcause","alphaEti",
                                    paste("JSS",1:length(JSS_list),sep="_"),
                                    paste("GSS_TPR",1:length(JSS_list),sep="_"),    # <-- added for TPR strata.
                                    paste("SS_TPR_grp",1:length(JSS_list),sep="_"), # <-- added for TPR strata.
                                    paste("MSS",1:length(JSS_list),sep="_"),
                                    paste("templateSS",1:length(JSS_list),sep="_"),
                                    paste("alphaS",1:length(JSS_list),sep="_"),   
                                    paste("betaS",1:length(JSS_list),sep="_")))
        } else{
          in_data       <- unique(c(in_data,"Nd","Nu","Jcause","alphaEti",
                                    paste("JSS",1:length(JSS_list),sep="_")[-single_column_MSS],
                                    paste("GSS_TPR",1:length(JSS_list),sep="_"),   # <-- added for TPR strata.
                                    paste("SS_TPR_grp",1:length(JSS_list),sep="_"),# <-- added for TPR strata.
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
            res_curr[[2]] <- stats::rbeta(JBrS_list[[s]],1,1)
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
            res_curr[[4]] <- 1 #<---- added together with 'alphadp0_case' below.
            
            names(res_curr) <- paste(c("r0","r1","alphadp0","alphadp0_case"),s,sep="_")
            res <- c(res,res_curr)
          }
        }
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
            res_curr[[2]] <- stats::rbeta(JBrS_list[[s]],1,1)
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
            res_curr[[4]] <- 1
            
            names(res_curr) <- paste(c("r0","r1","alphadp0","alphadp0_case"),s,sep="_")
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
        c(res,res2)
      }
    }
    
    # etiology (measurement independent)
    alphaEti          <- prior$Eti_prior    # <-------- input etiology prior here.
    
    #
    # fit model :
    #
    
    # special operations:
    # get individual prediction:
    if (!is.null(mcmc_options$individual.pred) && mcmc_options$individual.pred){out_parameter <- c(out_parameter,"Icat")}
    # get posterior predictive distribution of BrS measurments:
    if (!is.null(mcmc_options$ppd) && mcmc_options$ppd){out_parameter <- c(out_parameter,paste("MBS.new",seq_along(Mobs$MBS),sep = "_"))}
    
    #
    # write the .bug files into mcmc_options$bugsmodel.dir; 
    # could later set it equal to result.folder.
    #
    
    use_jags <- (!is.null(mcmc_options$use_jags) && mcmc_options$use_jags)
    model_func         <- write_model_NoReg(model_options$likelihood$k_subclass,
                                            data_nplcm$Mobs,
                                            model_options$prior,
                                            model_options$likelihood$cause_list,
                                            model_options$use_measurements,
                                            mcmc_options$ppd,
                                            use_jags)
    
    model_bugfile_name <- "model_NoReg.bug"
    
    filename <- file.path(mcmc_options$bugsmodel.dir, model_bugfile_name)
    writeLines(model_func, filename)
    #
    # run the model:
    #
    here <- environment()
    ##JAGS
    in_data.list <- lapply(as.list(in_data),get, envir= here)
    names(in_data.list) <- in_data
    #lapply(names(in_data.list), dump, append = TRUE, envir = here,
    #       file = file.path(mcmc_options$result.folder,"jagsdata.txt"))
    #do.call(file.remove, list(list.files(mcmc_options$result.folder, full.names = TRUE)))
    curr_data_txt_file <- file.path(mcmc_options$result.folder,"jagsdata.txt")
    if(file.exists(curr_data_txt_file)){file.remove(curr_data_txt_file)}
    dump(names(in_data.list), append = FALSE, envir = here,
         file = curr_data_txt_file)
    
    # ## fix dimension problem.... convert say .Dmi=7:6 to c(7,6) (an issue for templateBS_1):
    # bad_jagsdata_txt <- readLines(curr_data_txt_file)
    # #good_jagsdata_txt <- gsub( "([0-9]+):([0-9]+)", "c(\\1,\\2)", bad_jagsdata_txt,fixed = FALSE)
    # ## to add an additional complicatoin of dump creates a text file with , dim = but the JAGS only accepts .Dim=
    # good_jagsdata_txt <- gsub( ", dim =", ", .Dim=", 
    #                            gsub( "([0-9]+):([0-9]+)", "c(\\1,\\2)", bad_jagsdata_txt,fixed = FALSE),
    #                            fixed = FALSE)
    # 
    # 
    # writeLines(good_jagsdata_txt, curr_data_txt_file)

    ## fixed some problems of JAGS 4.3.2 not having cut function; and I(a,b) functions weirdly, even though
    ## the two elements are already constants (errors says that are not constant).
    curr_model_txt_file <- file.path(mcmc_options$result.folder,model_bugfile_name)
    bad_model_txt <- readLines(curr_model_txt_file)
    good_model_txt <- gsub( "cut\\(", "(", bad_model_txt,fixed = FALSE)
    good_model_txt <- gsub( "I\\(0\\.000001,0\\.999999\\)", " ", good_model_txt,fixed = FALSE)
    writeLines(good_model_txt, curr_model_txt_file)
    
    if(is.null(mcmc_options$jags.dir)){mcmc_options$jags.dir=""}
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
    )
    return(gs)
  }

#' Initialize individual latent status (for `JAGS`)
#' 
#' @details In `JAGS` 3.4.0, if an initial value contradicts the probabilistic specification, e.g.
#' `MSS_1[i,j] ~ dbern(mu_ss_1[i,j])`, where `MSS_1[i,j]=1` but `mu_ss_1[i,j]=0`,
#' then `JAGS` cannot understand it. In PERCH application, this is most likely used when the specificity of the 
#' silver-standard data is `1`. Note: this is not a problem in `WinBUGS`.
#' 
#' @param MSS_list A list of silver-standard measurement data, possibly with more than one 
#' slices; see `data_nplcm` argument in [nplcm()]
#' @param cause_list See `model_options` arguments in [nplcm()]
#' @param patho A vector of measured pathogen name for `MSS`; 
#' default is `colnames(MSS)`
#' 
#' @return a list of numbers, indicating categories of individual latent causes.
#' 
#' @family initialization functions
#'    
init_latent_jags_multipleSS <- function(MSS_list,cause_list,
                                        patho=unlist(lapply(MSS_list,colnames))){
  # <--- revising for multiple silver-standard data.
  #table(apply(MSS,1,function(v) paste(v,collapse="")))
  MSS <- do.call(cbind,MSS_list)
  ind_positive <- which(apply(MSS,1,sum,na.rm=TRUE)>0)
  res <- sample.int(length(cause_list),size = nrow(MSS), replace=TRUE)
  if (length(ind_positive)>0){
    vec <- sapply(ind_positive, function(i) paste(unique(patho[which(MSS[i,]==1)]),collapse="+"))
    res[ind_positive] <- match_cause(cause_list,vec)
  }
  if (sum(is.na(res[ind_positive]))>0){ # <--- corrected to add res[].
    ind_NA <- which(is.na(res))
    stop(paste0("==[baker] Case(s) The ",paste(ind_NA,collapse=", "), "-th subject(s) have positive silver-standard
                measurements on causative agent combinations that are not specified in the 'cause_list' of 
                'model_options$likelihood'! Please consider 1) deleting these cases,
                or 2) adding these combinations into 'cause_list'.==\n"))
  }
  res
}

# MSS_list <- data_nplcm$Mobs$MSS
# cause_list <- model_options$likelihood$cause_list
# patho <- unlist(lapply(MSS_list,colnames))
# 
# init_latent_jags_multipleSS(MSS_list,cause_list)
# 





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
#' @seealso [write_model_NoReg] for automatically generate `.bug` model
#'file; This present function store it in location: `mcmc_options$bugsmodel.dir`.
#' 
#' @family model fitting functions 
#' 
#'    
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
        xx <- as.matrix_or_vec(MBS_list[[i]]);  attr(xx,"dimnames") <- NULL # remove dimnames.
        assign(paste("MBS", i, sep = "_"), xx)
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
    here <- environment()
    
    ##JAGS
    in_data.list <- lapply(as.list(in_data), get, envir= here)
    names(in_data.list) <- in_data
    #lapply(names(in_data.list), dump, append = TRUE, envir = here,
    #       file = file.path(mcmc_options$result.folder,"jagsdata.txt"))
    curr_data_txt_file <- file.path(mcmc_options$result.folder,"jagsdata.txt")
    if(file.exists(curr_data_txt_file)){file.remove(curr_data_txt_file)}
    dump(names(in_data.list), append = FALSE, envir = here,
         file = curr_data_txt_file)
    
    ## fixed some problems of JAGS 4.3.2 not having cut function; and I(a,b) functions weirdly, even though
    ## the two elements are already constants (errors says that are not constant).
    curr_model_txt_file <- file.path(mcmc_options$result.folder,model_bugfile_name)
    bad_model_txt <- readLines(curr_model_txt_file)
    good_model_txt <- gsub( "cut\\(", "(", bad_model_txt,fixed = FALSE)
    good_model_txt <- gsub( "I\\(0\\.000001,0\\.999999\\)", " ", good_model_txt,fixed = FALSE)
    writeLines(good_model_txt, curr_model_txt_file)
    
    if(is.null(mcmc_options$jags.dir)){mcmc_options$jags.dir=""}
    gs <- jags2_baker(data   = curr_data_txt_file,
                      inits  = xxx,
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





#' Fit nested partially-latent class model with regression (low-level)
#'
#' @details This function prepares data, specifies hyperparameters in priors 
#' (true positive rates and CSCFs), initializes the posterior
#' sampling chain, writes the model file (for JAGS or WinBUGS with slight
#' differences in syntax), and fits the model. Features:
#' \itemize{
#' \item regression (not all discrete covariates);
#' \item no nested subclasses, i.e. conditional independence of 
#' multivariate measurements given disease class and covariates;
#' \item multiple BrS + multiple SS.
#' }
#'
#' @inheritParams nplcm
#' @return BUGS fit results from `JAGS`.
#' 
#' @seealso [write_model_NoReg] for constructing `.bug` model file; This function
#' then puts it in the folder `mcmc_options$bugsmodel.dir`.
#' 
#' @family model fitting functions 
#' 
#'    
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
    
    
    if(is.null(prior$Eti_hyper_pflex)){
      prior$Eti_hyper_pflex = c(1,0.5)
      # beta prior for selecting constant vs non-constant additive component
    }
    
    if(is.null(prior$FPR_coef_prior)){
      prior$FPR_coef_prior = c(2,2)
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
        as.numeric(0+(1:ncol(Z_Eti)))
      } else{
        as.numeric((1:ncol(Z_Eti))[-grep("^s_",dimnames(Z_Eti)[[2]])])
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
          as.numeric(0+(1:ncol(Z)))
        } else{
          as.numeric((1:ncol(Z))[-grep("^s_",dimnames(Z)[[2]])])
        }
      })
      
      #Z_FPR_list_orth <- list()
      for(i in seq_along(JBrS_list)){
        attributes(Z_FPR_list[[i]])[names(attributes(Z_FPR_list[[i]]))!="dim"] <- NULL
        
        assign(paste("JBrS", i, sep = "_"), JBrS_list[[i]])    
        xx <- as.matrix_or_vec(MBS_list[[i]]);  attr(xx,"dimnames") <- NULL # remove dimnames.
        assign(paste("MBS", i, sep = "_"), xx)
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
      
      ## fixed some problems of JAGS 4.3.2 not having cut function; and I(a,b) functions weirdly, even though
      ## the two elements are already constants (errors says that are not constant).
      curr_model_txt_file <- file.path(mcmc_options$result.folder,model_bugfile_name)
      bad_model_txt <- readLines(curr_model_txt_file)
      good_model_txt <- gsub( "cut\\(", "(", bad_model_txt,fixed = FALSE)
      good_model_txt <- gsub( "I\\(0\\.000001,0\\.999999\\)", " ", good_model_txt,fixed = FALSE)
      writeLines(good_model_txt, curr_model_txt_file)
      
      if(is.null(mcmc_options$jags.dir)){mcmc_options$jags.dir=""}
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
#' Called by [nplcm()] upon being assigned to this nested regression by 
#' [assign_model()]
#'
#' @details This function prepares data, specifies hyperparameters in priors 
#' (true positive rates and etiology fractions), initializes the posterior
#' sampling chain, writes the model file (for `JAGS`), and fits the model. Features:
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
#' @seealso [write_model_Reg_Nest] for constructing `.bug` model file; This function
#' then put it in the folder `mcmc_options$bugsmodel.dir`.
#' 
#' @family model fitting functions 
#' 
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
  
  
  if(is.null(prior$half_nu_s2)){
    prior$half_nu_s2 = c(1/2,100/2) # hyper-parameters 
  }
  
  if(is.null(prior$Eti_hyper_pflex)){
    prior$Eti_hyper_pflex = c(1,0.5)
    # beta prior for selecting constant vs non-constant additive component
  }
  
  if(is.null(prior$FPR_hyper_pflex)){
    prior$FPR_hyper_pflex = c(0.5,1) # same as above, but for FPR regression in the controls.
  }
  
  if(is.null(prior$case_FPR_hyper_pflex)){
    prior$case_FPR_hyper_pflex = c(0.5,1) # same as above, but for FPR regression in the cases.
  }
  
  if(is.null(prior$FPR_coef_prior)){
    prior$FPR_coef_prior = c(2,2)
  }
  
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
      as.numeric(0+(1:ncol(Z_Eti)))
    } else{
      as.numeric((1:ncol(Z_Eti))[-grep("^s_",dimnames(Z_Eti)[[2]])])
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
        as.numeric( 0+(1:ncol(Z)))
      } else{
        as.numeric( (1:ncol(Z))[-grep("^s_",dimnames(Z)[[2]])])
      }
    })
    
    #Z_FPR_list_orth <- list()
    for(i in seq_along(JBrS_list)){
      attributes(Z_FPR_list[[i]])[names(attributes(Z_FPR_list[[i]]))!="dim"] <- NULL
      
      assign(paste("JBrS", i, sep = "_"), JBrS_list[[i]])    
      xx <- as.matrix_or_vec(MBS_list[[i]]);  attr(xx,"dimnames") <- NULL # remove dimnames.
      assign(paste("MBS", i, sep = "_"), xx)
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
        #cat("==[baker] Running etiology regression with nested subclasses...==\n")
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
  ctrl_flex_alpha    <- prior$FPR_hyper_pflex[1]
  ctrl_flex_beta     <- prior$FPR_hyper_pflex[2]
  case_flex_alpha    <- prior$case_FPR_hyper_pflex[1]
  case_flex_beta     <- prior$case_FPR_hyper_pflex[2]
  
  for (s in seq_along(Mobs$MBS)){
    if(!is.null(prior$half_nu_s2)){
      assign(paste("half_nu_ctrl", i, sep = "_"), prior$half_nu_s2[1])
      assign(paste("half_s2_ctrl", i, sep = "_"), prior$half_nu_s2[2])
    }
  }
  in_data <- unique(c(in_data,
                      "ER_basis_id"[ER_has_basis],
                      "ER_n_basis"[ER_has_basis],
                      "ER_non_basis_id"[ER_has_non_basis],
                      "sd_betaEti_basis"[ER_has_basis],
                      "ER_alpha"[ER_has_basis],
                      "ER_beta"[ER_has_basis],
                      "sd_betaEti_nonbasis"[ER_has_non_basis],
                      "sd_betaFPR_basis"[any(unlist(has_basis_list))],
                      "sd_betaFPR_nonbasis"[any(unlist(has_non_basis_list))],
                      "ctrl_flex_alpha"[any(unlist(has_basis_list))],
                      "ctrl_flex_beta"[any(unlist(has_basis_list))],
                      "case_flex_alpha"[any(unlist(has_basis_list))],
                      "case_flex_beta"[any(unlist(has_basis_list))],
                      paste("half_nu_ctrl",1:length(JBrS_list),sep="_"),
                      paste("half_s2_ctrl",1:length(JBrS_list),sep="_")#,
                      #paste("half_s2_case",1:length(JBrS_list),sep="_")[unlist(has_basis_list)]
                      
                      
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
  here <- environment()
  ##JAGS
  in_data.list <- lapply(as.list(in_data), get, envir= here)
  names(in_data.list) <- in_data
  #lapply(names(in_data.list), dump, append = TRUE, envir = here,
  #       file = file.path(mcmc_options$result.folder,"jagsdata.txt"))
  curr_data_txt_file <- file.path(mcmc_options$result.folder,"jagsdata.txt")
  if(file.exists(curr_data_txt_file)){file.remove(curr_data_txt_file)}
  dump(names(in_data.list), append = FALSE, envir = here,
       file = curr_data_txt_file)

  
  ## fixed some problems of JAGS 4.3.2 not having cut function; and I(a,b) functions weirdly, even though
  ## the two elements are already constants (errors says that are not constant).
  curr_model_txt_file <- file.path(mcmc_options$result.folder,model_bugfile_name)
  bad_model_txt <- readLines(curr_model_txt_file)
  good_model_txt <- gsub( "cut\\(", "(", bad_model_txt,fixed = FALSE)
  good_model_txt <- gsub( "I\\(0\\.000001,0\\.999999\\)", " ", good_model_txt,fixed = FALSE)
  writeLines(good_model_txt, curr_model_txt_file)
  
  if(is.null(mcmc_options$jags.dir)){mcmc_options$jags.dir=""}
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

