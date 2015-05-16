#' Obtain MSE, bias-squared and variance of posterior from nplcm fitted folders.
#'
#' The result is equivalent to Euclidean-type calculation after the
#' compositional vector (e.g., etiologic fraction) is centered-logratio tranformed.
#' For simulation only.
#'
#' @param DIR_NPLCM File path where Bayesian results are stored
#' @param truth Default is \code{NULL}. True etiologic fraction vector (must sum to 1) used to generate data
#' @param check.truth Default is \code{FALSE}. If \code{truth} is a non-null vector,
#' it asks the function to check the user specified truth against the one used
#' in the actual simulation.
#' @importFrom robCompositions cenLR
#' @importFrom coda read.coda
#' @return a vector of (MSE, bias-squared, variance,truth)
#' @export
#'


get_mse <- function(DIR_NPLCM,truth=NULL,check.truth=FALSE){

  my.mse <- function(A,B){
    if (nrow(A)!=nrow(B)){
      stop("not equal number of rows!")
    }else{
      aa=cenLR(A)$x.clr
      bb=cenLR(B)$x.clr
      res <- rowSums((aa-bb)^2)
      mean(res)

    }
  }

  my.bias.sq <- function(A,B){
    if (nrow(A)!=nrow(B)){
      stop("not equal number of rows!")
    }else{
      aa=cenLR(A)$x.clr
      bb=cenLR(B)$x.clr
      sum((colMeans(aa-bb))^2)
    }
  }

  #A = as.data.frame(true_pEti)
  my.var <- function(A){

    res <- cenLR(A)$x.clr
    sum(diag(cov(res)))

  }

  if (!file.exists(DIR_NPLCM)){
    stop("==No such folder!==")
  }else{
    #read in data from the folder directory:
    bugs.dat <- dget(paste(DIR_NPLCM,"data.txt",sep="/"))
    for (bugs.variable.name in names(bugs.dat)) {
      if (!is.null(dim(bugs.dat[[bugs.variable.name]]))) {
        dim(bugs.dat[[bugs.variable.name]]) <- rev(dim(bugs.dat[[bugs.variable.name]]))
        bugs.dat[[bugs.variable.name]] <- aperm(bugs.dat[[bugs.variable.name]])
      }
      assign(bugs.variable.name, bugs.dat[[bugs.variable.name]])
    }

    model_options  <- dget(paste(DIR_NPLCM,"model_options.txt",sep="/"))

    if (is.null(truth)){
      print("==No truth specified; Looking in the result folder!==")
      if (!file.exists(paste(DIR_NPLCM,"set_parameter.txt",sep="/"))){
        stop("==No set_parameter.txt available in result folder to check
                  equality of specified and used true etiology!==")
      }else{
        set_parameter <- dget(paste(DIR_NPLCM,"set_parameter.txt",sep="/"))
        use.truth <- set_parameter$etiology
      }
    }else{
        if (check.truth){
          if (!file.exists(paste(DIR_NPLCM,"set_parameter.txt",sep="/"))){
            stop("==No set_parameter.txt available in result folder to check
                  equality of specified and used true etiology!==")
          } else{
            set_parameter <- dget(paste(DIR_NPLCM,"set_parameter.txt",sep="/"))
            if (sum(set_parameter$etiology-truth)>0){
              stop("==Specified true etiology is different than the actual one that
               generates the simulated data!==")
            }
          }
        }
      use.truth <- truth
    }
    #compatibility checking:
    if (length(model_options$M_use)!=length(model_options$TPR_prior)){
      stop("The number of measurement source(s) is different from
             the number of TPR prior option!
             Make them equal, and match with order!")
    }

    #some data preparation:
    Nd <- bugs.dat$Nd
    Nu <- bugs.dat$Nu

    Y = c(rep(1,Nd),rep(0,Nu))

    pathogen_list     <- model_options$pathogen_list

    Jfull_BrS         <- length(pathogen_list)

    #reading nplcm outputs:
    res_nplcm <- read.coda(paste(DIR_NPLCM,"coda1.txt",sep="/"),
                           paste(DIR_NPLCM,"codaIndex.txt",sep="/"),
                           quiet=TRUE)

    Jfull <- length(grep("pEti",colnames(res_nplcm)))

    # extract and process some data and posterior samples:
    SubVarName <- rep(NA,Jfull)
    for (j in 1:Jfull){
      SubVarName[j] = paste("pEti","[",j,"]",sep="")
    }

    #get etiology fraction MCMC samples:
    pEti_mat   <- as.matrix(res_nplcm[,SubVarName])


    true_pEti <- as.data.frame(matrix(use.truth,
                                      nrow(pEti_mat),ncol(pEti_mat),
                                      byrow=TRUE))
    # Aitchison distances:
    #aDist(true_pEti,pEti_mat)
    mse <- my.mse(true_pEti,pEti_mat)
    bias.sq <- my.bias.sq(true_pEti,pEti_mat)
    vv <- my.var(as.data.frame(pEti_mat))

    res <- list(mse=mse,bias.sq=bias.sq,vv=vv,truth=use.truth)
    res
  }

}





