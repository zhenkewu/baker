#' Read data and other model information from a folder that stores model results.
#'
#' @param DIR_NPLCM File path to the folder containing posterior samples
#' 
#' @return A list with data, options and posterior samples.
#' \itemize{
#' \item \code{bugs.dat}
#' \item \code{model_options}
#' \item \code{clean_otions}
#' \item \code{Nd}; \code{Nu}; \code{Y}; \code{Mobs}; 
#' \item \code{res_nplcm}.
#' }
#'
#' @export

nplcm_read_folder <- function(DIR_NPLCM){
  use_jags      <- is_jags_folder(DIR_NPLCM)
  #
  # Read data from DIR_NPLCM:
  #
  if (use_jags){
    new_env <- new.env()
    source(file.path(DIR_NPLCM,"jagsdata.txt"),local=new_env)
    bugs.dat <- as.list(new_env)
    rm(new_env)
    res_nplcm <- coda::read.coda(file.path(DIR_NPLCM,"CODAchain1.txt"),
                                file.path(DIR_NPLCM,"CODAindex.txt"),
                                quiet=TRUE)
    for (bugs.variable.name in names(bugs.dat)) {
      assign(bugs.variable.name, bugs.dat[[bugs.variable.name]])
    }
    
  }else{
    bugs.dat <- dget(file.path(DIR_NPLCM,"data.txt"))  
    res_nplcm <- coda::read.coda(file.path(DIR_NPLCM,"coda1.txt"),
                                 file.path(DIR_NPLCM,"codaIndex.txt"),
                                 quiet=TRUE)
    # WinBUGS stores objects with transposition, so need to rev or aperm back:
    for (bugs.variable.name in names(bugs.dat)) {
        if (!is.null(dim(bugs.dat[[bugs.variable.name]]))) {
          dim(bugs.dat[[bugs.variable.name]]) <- rev(dim(bugs.dat[[bugs.variable.name]]))
          bugs.dat[[bugs.variable.name]] <- aperm(bugs.dat[[bugs.variable.name]])
        }
      assign(bugs.variable.name, bugs.dat[[bugs.variable.name]])
    }
  }

  model_options  <- dget(file.path(DIR_NPLCM,"model_options.txt"))
  if (!file.exists(file.path(DIR_NPLCM,"data_clean_options.txt"))){
    stop("=='data_clean_options.txt' does not exist in the result folder. Please 'dput' the clean_options in the result folder. ==")
  } 
  clean_options  <- dget(file.path(DIR_NPLCM,"data_clean_options.txt"))
  #some data preparation:
  Nd <- bugs.dat$Nd
  Y  <- c(rep(1,Nd))
  Nu <- 0
  if (!is.null(bugs.dat$Nu)){
    Nu <- bugs.dat$Nu
    Y  <- c(rep(1,Nd),rep(0,Nu))
  }
  
  
  get_MBS <- function(){
      MBS_nm <- (names(bugs.dat)[grep("MBS_",names(bugs.dat))])
      MBS_nm <- MBS_nm[order(MBS_nm)]
      res <- list()
      for (s in seq_along(MBS_nm)){
        res[[s]] <- as.data.frame(bugs.dat[[MBS_nm[s]]])
        colnames(res[[s]]) <- clean_options$BrS_objects[[s]]$patho
      }
      names(res) <- unlist(lapply(clean_options$BrS_objects,"[[","nm_spec_test"))
      res
  }
  
  get_MSS <- function(){
    MSS_nm <- (names(bugs.dat)[grep("MSS_",names(bugs.dat))])
    MSS_nm <- MSS_nm[order(MSS_nm)]
    res <- list()
    for (s in seq_along(MSS_nm)){
      res[[s]] <- as.data.frame(bugs.dat[[MSS_nm[s]]])
      colnames(res[[s]]) <- clean_options$SS_objects[[s]]$patho
    }
    names(res) <- unlist(lapply(clean_options$SS_objects,"[[","nm_spec_test"))
    res
  }
  
  MBS <- MSS <- MGS <- NULL
  
  if ("BrS" %in% model_options$use_measurements){MBS <- get_MBS()}
  if ("SS" %in% model_options$use_measurements){MSS <- get_MSS()}
  
  Mobs <- make_list(MBS, MSS, MGS)
  

  res <- list(bugs.dat = bugs.dat,
              model_options = model_options,
              clean_options = clean_options,
              Nd = Nd,
              Nu = Nu,
              Y  = Y,
              Mobs = Mobs,
              res_nplcm = res_nplcm)
}


#' get individual prediction (Bayesian posterior)
#' 
#' @param read_res a list read from a folder that stores fitted model results; 
#' Commonly, it is the returned value from function \code{\link{nplcm_read_folder}}
#' 
#' @return a matrix of individual predictions; rows for cases, columns for causes 
#' specified in \code{model_options$likelihood$cause_list}; See \code{\link{nplcm}}
#' 
#' @examples 
#' \dontrun{
#' get_individual_prediction(nplcm_read_folder("C://2015_09_02_02GAM"))
#' }
#' 
#' @export
#' 
get_individual_prediction <- function(read_res){
  res_nplcm     <- read_res$res_nplcm
  model_options <- read_res$model_options
  NSAMP         <- nrow(res_nplcm)
  NCAUSE        <- length(model_options$likelihood$cause_list)
  ind_pred_mat  <- matrix(0,nrow=read_res$Nd,ncol=NCAUSE)
  Nd            <- read_res$Nd
  
  for (i in 1:Nd){
    curr_Icat <- res_nplcm[,grep(paste0("^Icat\\[",i,"\\]"),colnames(res_nplcm))]
    curr_ind_pred_mat <- matrix(0,nrow=NSAMP,ncol=NCAUSE)
    for (iter in 1:NSAMP){
      curr_ind_pred_mat[iter,curr_Icat[iter]] <- 1  
    }
    ind_pred_mat[i,] <- colMeans(curr_ind_pred_mat)
  }
  colnames(ind_pred_mat) <- model_options$likelihood$cause_list
  ind_pred_mat
}

#' get individual data
#' 
#' @param i index of individual as appeared in \code{data_nplcm}
#' @param data_nplcm the data for nplcm; see \code{\link{nplcm}}
#' 
#' @return a list of the same structure as \code{data_nplcm}; just with one row of values
#' 
#' @export
#' 
get_individual_data <- function(i, data_nplcm){
  MBS <- NULL
  if (!is.null(data_nplcm$Mobs$MBS)){
    MBS <- vector("list",length = length(data_nplcm$Mobs$MBS))
    names(MBS) <- names(data_nplcm$Mobs$MBS)
    for (s in seq_along(data_nplcm$Mobs$MBS)){
      MBS[[s]] <- data_nplcm$Mobs$MBS[[s]][i,,drop=FALSE]
    }
  }
  
  MSS <- NULL
  if (!is.null(data_nplcm$Mobs$MSS)){
    MSS <- vector("list",length = length(data_nplcm$Mobs$MSS))
    names(MSS) <- names(data_nplcm$Mobs$MSS)
    for (s in seq_along(data_nplcm$Mobs$MSS)){
      MSS[[s]] <- data_nplcm$Mobs$MSS[[s]][i,,drop=FALSE]
    }
  }
  
  MGS <- NULL
  if (!is.null(data_nplcm$Mobs$MGS)){
    MGS <- vector("list",length = length(data_nplcm$Mobs$MGS))
    names(MGS) <- names(data_nplcm$Mobs$MGS)
    for (s in seq_along(data_nplcm$Mobs$MGS)){
      MGS[[s]] <- data_nplcm$Mobs$MGS[[s]][i,,drop=FALSE]
    }
  }
  
  X <- NULL
  if (!is.null(data_nplcm$X)){
    X <- data_nplcm$X[i,,drop=FALSE]
  }
  
  Y <- NULL
  if (!is.null(data_nplcm$Y)){
    Y <- data_nplcm$Y[i]
  }
  list(Mobs = list(MBS = MBS, MSS=MSS, MGS=MGS), X=X,Y=Y)
}


#' See if a result folder is obtained by JAGS
#' 
#' 
#' @param DIR_NPLCM directory to the folder with results. 
#' "mcmc_options.txt" must be in the folder.
#' 
#' @return TRUE for from JAGS; FALSE otherwise.
#' 
#' @export

is_jags_folder <- function(DIR_NPLCM){
  mcmc_options  <- dget(file.path(DIR_NPLCM,"mcmc_options.txt"))
  !is.null(mcmc_options$use_jags) && mcmc_options$use_jags
}

