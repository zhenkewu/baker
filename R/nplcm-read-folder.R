#' Read data and other model information from a folder that stores model results.
#'
#' @param DIR_NPLCM File path to the folder containing posterior samples
#' 
#' @return A list with data, options and posterior samples.
#'
#' @export

nplcm_read_folder <- function(DIR_NPLCM){
  #
  # Read data from DIR_NPLCM:
  #
  bugs.dat <- dget(paste(DIR_NPLCM,"data.txt",sep="/"))
  for (bugs.variable.name in names(bugs.dat)) {
    if (!is.null(dim(bugs.dat[[bugs.variable.name]]))) {
      dim(bugs.dat[[bugs.variable.name]]) <- rev(dim(bugs.dat[[bugs.variable.name]]))
      bugs.dat[[bugs.variable.name]] <- aperm(bugs.dat[[bugs.variable.name]])
    }
    assign(bugs.variable.name, bugs.dat[[bugs.variable.name]])
  }
  
  model_options  <- dget(file.path(DIR_NPLCM,"model_options.txt"))
  if (!file.exists(file.path(DIR_NPLCM,"data_clean_options.txt"))){
    stop("=='data_clean_options.txt' does not exist in the result folder. Please 'dput' the clean_options in the result folder. ==")
  } 
  clean_options  <- dget(file.path(DIR_NPLCM,"data_clean_options.txt"))
  #some data preparation:
  Nd <- bugs.dat$Nd
  Nu <- bugs.dat$Nu
  Y  <- c(rep(1,Nd),rep(0,Nu))
  
  
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
  
  res_nplcm <- coda::read.coda(file.path(DIR_NPLCM,"coda1.txt"),
                         file.path(DIR_NPLCM,"codaIndex.txt"),
                         quiet=TRUE)
  res <- list(bugs.dat = bugs.dat,
              model_options = model_options,
              clean_options = clean_options,
              Nd = Nd,
              Nu = Nu,
              Y  = Y,
              Mobs = Mobs,
              res_nplcm = res_nplcm)
}