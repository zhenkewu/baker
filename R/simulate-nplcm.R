#' Simulation data from nested partially-latent class model (npLCM) family
#'
#' Use different case and control subclass mixing weights. Eta is of
#' dimension J times K. NB: document the elements in \code{set_parameter}. Also, current
#' function is written in a way to facilitate adding more measurement components.
#' 
#' @param set_parameter True model parameters in the npLCM specification
#' 
#'  
#' @return A list of measurements, true latent statues:
#' \itemize{ 
#'  \item{\code{template}} a matrix: rows for causes, columns for measurements; 
#'  generated as a lookup table to match mixture component parmeters for every type
#'   (a particular cause) of indiviuals. 
#'  \item{\code{data_nplcm}} a list of structured data (see \code{\link{nplcm}} for 
#'  description) for use in visualization 
#'  e.g., \code{\link{plot_logORmat}} or model fitting, e.g., \code{\link{nplcm}}. 
#'  The pathogen taxonomy is set to default "B".
#'  \item{\code{latent_cat}} integer values to indicate the latent category. The integer
#'  code corresponds to the order specified in \code{set_parameter$etiology}. 
#'  Controls are coded as \code{length(set_parameter$etiology)+1}.)
#'  }
#'  
#' @examples 
#' \dontrun{
#' K.true  <- 2   # no. of latent subclasses in actual simulation.
#' J       <- 5   # no. of pathogens.
#' N      <- 10000
#' eta_seq      <- seq(0,0.5,by=0.05) 
#' 
#' eta <- 0
#' 
#' set_parameter <- list(
#'   pathogen_BrS    = LETTERS[1:J],
#'   cause_list      = c(LETTERS[1:J],"A+B","D+E"),
#'   etiology        = c(0.5,0.1,0.15,0.1,0.05,0.05,0.05), #same length as cause_list
#'   Lambda          = c(1-eta,eta), #ctrl mix
#'   Eta             = rbind(c(1-eta,eta),
#'                           c(1-eta,eta),
#'                           c(1-eta,eta),
#'                           c(1-eta,eta),
#'                           c(1-eta,eta),
#'                           c(1-eta,eta),
#'                           c(1-eta,eta)), #case mix, row number equal to Jcause.
#'   PsiBS 		  =    matrix(c(0.1,0.3,
#'                           0.1,0.3,
#'                           0.1,0.3,
#'                           0.1,0.3,
#'                           0.1,0.3),nrow=J,ncol=K.true,byrow=TRUE),
#'  ThetaBS         =  matrix(c(0.9,0.1,
#'                               0.9,0.1,
#'                               0.9,0.1,
#'                               0.9,0.1,
#'                               0.9,0.1),nrow=J,ncol=K.true,byrow=TRUE),
#'   Nu      =     N, # control size.
#'   Nd      =     N  # case size.
#' )
#' 
#'  pathogen_display <- data_nplcm$Mname$Mname_BrS
#'  plot_logORmat(data_nplcm,pathogen_display)
#' }
#'  
#' @export

simulate_nplcm <- function(set_parameter){
  pathogen_BrS <- set_parameter$pathogen_BrS
  
  # simulate latent status  
  latent <- simulate_latent(set_parameter)
  # simulate BrS measurements:
  out_brs   <- simulate_brs(set_parameter,latent)
  
  MBS_list <- list(MBS_1 = out_brs[,-grep("Y",colnames(out_brs)),drop=FALSE])
  # construct a list for visualization and model fitting:
  Mobs <- list(MBS = MBS_list ,MSS=NULL,MGS=NULL)
  Y    <- MBS$Y
  X    <- NULL
  Mname<- list(Mname_BrS    = pathogen_BrS,
               Mname_SSonly = NULL)
  taxonomy <- list(taxo_BrS = rep("B",length(pathogen_BrS)),
                   taxo_SSonly = NULL)
  
  data_nplcm <- make_list(Mobs, Y, X, Mname, taxonomy )
  template  <- latent$template
  latent_cat <- latent$iLcatAllnumeric
  make_list(template, data_nplcm,latent_cat)
}

#' Simulate Latent Status:
#' @param set_parameter parameters for measurements
#' 
#' @return a list of latent status samples for use in sampling measurements. It
#' also includes a template to look up measurement parameters for each type of causes.
#' @export
#' 
simulate_latent <- function(set_parameter){
  # etiology common to all measurements:
  cause_list <- set_parameter$cause_list
  etiology   <- set_parameter$etiology
  Jcause     <- length(cause_list)
  
  # sample size:
  Nd       <- set_parameter$Nd
  Nu       <- set_parameter$Nu
  
  # simulate latent status (common to all measurements):
  iLcat <- rep(NA,Nd)
  iLall <- matrix(NA,nrow=Nd+Nu,ncol=Jcause)
  etiologyMat <- matrix(NA,nrow=Nd,ncol=Jcause)
  
  # sample cause for cases:
  for (i in 1:Nd){
    etiologyMat[i,] <- etiology
    iLcat[i]        <- sample(cause_list,1,prob = etiologyMat[i,])
  }
  
  
  pathogen_BrS    <- set_parameter$pathogen_BrS
  J_BrS  <- length(pathogen_BrS)
  
  # convert categorical to template (cases):
  iL    <- symb2I(iLcat,pathogen_BrS)
  # convert back to categorical (cases):
  iLcat.case.numeric <- Imat2cat(iL,cause_list,pathogen_BrS)
  # create 
  iLall <- rbind(iL,matrix(0,nrow=Nu,ncol=J_BrS))
  iLcatAllnumeric    <- c(iLcat.case.numeric,rep(Jcause+1,Nu))
  
  template <- as.matrix(rbind(symb2I(cause_list,pathogen_BrS),rep(0,J_BrS)))
  colnames(template) <- pathogen_BrS
  rownames(template) <- c(cause_list,"control")
  
  make_list(iLall,iLcatAllnumeric,iLcat.case.numeric,iL,template)
}

#' Simulate Bronze-Standard Data
#' 
#' 
#' simulate BrS measurements:
#' @param set_parameter parameters for BrS measurements
#' @param latent_samples sampled latent status for all the subjects, for use in simulate
#' BrS measurements.
#' 
#' @return a data frame with first column being case-control status (case at top) and
#' columns of bronze-standard measurements
#' @export
simulate_brs <- function(set_parameter,latent_samples){
  
  pathogen_BrS    <- set_parameter$pathogen_BrS
  J_BrS           <- length(pathogen_BrS)
  PsiBS           <- set_parameter$PsiBS
  ThetaBS         <- set_parameter$ThetaBS
  Lambda          <- set_parameter$Lambda
  Eta             <- set_parameter$Eta
  
  iLall              <- latent_samples$iLall
  iLcatAllnumeric    <- latent_samples$iLcatAllnumeric
  iLcat.case.numeric <- latent_samples$iLcat.case.numeric
  iL <- latent_samples$iL
  # sample size:
  Nd       <- set_parameter$Nd
  Nu       <- set_parameter$Nu
  
  Zd <- rep(NA,Nd)
  Md <- matrix(NA,nrow=Nd,ncol=J_BrS)
  MdP <- Md
  for (i in 1:Nd){
    Zd[i] = sample(1:ncol(Eta),1,prob = Eta[iLcat.case.numeric[i],])
    for (j in 1:J_BrS){
      MdP[i,j]  = PsiBS[j,Zd[i]]*(1-iL[i,j])+iL[i,j]*ThetaBS[j,Zd[i]]
    }
  }
  Md <- rvbern(MdP)
  
  Zu  <- rep(NA,Nu)
  Mu  <- matrix(NA,nrow=Nu,ncol=J_BrS)
  MuP <- matrix(NA,nrow=Nu,ncol=J_BrS)
  for (i in 1:Nu){
    Zu[i]     <- sample(1:length(Lambda),1,prob = Lambda)
    for (j in 1:J_BrS){
      MuP[i,j]  <- PsiBS[j,Zu[i]]
    }
  }
  Mu <- rvbern(MuP)
  
  ## organize case/control status, iL, BS, GS data into dataframes
  datacolnames    <- c("Y", paste("MBS",pathogen_BrS,sep="_"))
  #   datres <- data.frame(Y = c(rep(1,Nd),rep(0,Nu)),
  #                        iLcat = iLcatAllnumeric,
  #                        iL = iLall,
  #                        MBS = rbind(Md,Mu))
  #  colnames(datres) <- datacolnames
  
  # dat_meas <- datres[,rev(rev(1:ncol(datres))[1:J_BrS]),drop=FALSE]
  #dat_case <- as.matrix(dat_meas[(1:set_parameter$Nd),])
  #dat_ctrl <- as.matrix(dat_meas[-(1:set_parameter$Nd),])
  
  datres <- data.frame(Y = c(rep(1,Nd),rep(0,Nu)), MBS = rbind(Md,Mu))
  colnames(datres) <- datacolnames
  datres
}


