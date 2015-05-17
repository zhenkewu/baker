#' Read data and other model information from a folder that stores model results.
#'
#' @param DIR_NPLCM File path to the folder containing posterior samples
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
  
  Mobs <- list(MBS = bugs.dat$MBS,
               MSS = bugs.dat$MSS,
               MGS = bugs.dat$MGS)
  
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