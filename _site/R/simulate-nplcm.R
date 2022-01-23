#' Simulate data from nested partially-latent class model (npLCM) family
#'
#' 
#' 
#' 
#' @param set_parameter True model parameters in an npLCM specification:
#'  \itemize{
#'     \item{`cause_list`} a vector of disease class names among cases (since
#'     the causes could be multi-agent (e.g., multiple pathogens may cause an individual case's
#'     pneumonia), so its length could be longer than the total number of unique
#'     causative agents)
#'     \item{`etiology`} a vector of proportions that sum to 100 percent
#'     \item{`pathogen_BrS`} a vector of putative causative agents' names measured in bronze-standard (BrS) data.
#'     This function simulates only one slice defined by {specimen}{test}{pathogen}
#'     \item{`pathogen_SS`} a vector of pathogen names measured in silver-standard (SS) data.
#'     \item{`meas_nm`} a list of {specimen}{test} names e.g., `list(MBS = c("NPPCR"),MSS="BCX")`
#'     for nasopharyngeal (NP) specimen tested by polymerase chain reaction (PCR) - `NPPCR` and
#'      blood (B) tested by culture (Cx) - `BCX`
#'     \item{`Lambda`} controls' subclass weights \eqn{\nu_1, \nu_2, \ldots, \nu_K}
#'     a vector of `K` probabilities that sum to 1.
#'     \item{`Eta`}    a matrix of dimension `length(cause_list)` by `K`;
#'     each row represents a disease class (among cases); the values in that row
#'     are subclass weights \eqn{\eta_1, \eta_2, \ldots, \eta_K} for that disease class,
#'     so needs to sum to one. In Wu et al. 2016 (JRSS-C), the subclass weights are the same across disease
#'     classes across rows. But when simulating data, one can specify rows with distinct
#'     subclass weights - it is a matter whether we can recover these parameters (possible when
#'     some cases' true disease classes are observed)
#'     \item{`PsiBS/PsiSS`} False positive rates for Bronze-Standard data and 
#'     for Silver-Standard data. For example, the rows of `PsiBS` correspond to the dimension of the particular
#'     slice of BrS measures, e.g., `10` for 10 causative agents measured by NPPCR; the
#'     columns correspond to `K` subclasses; generically, the dimension is `J` by `K`
#'     `PsiSS` is supposed to be a vector of all zeros (perfect specificity in silver-standard measures).
#'     \item{`ThetaBS/ThetaSS`}  True positive rates \eqn{\Theta} for Bronze-Standard data and 
#'     for Silver-Standard data. Dimension is `J` by `K` (can contain `NA` if the total number of 
#'     causative agents measured by BrS or SS exceeds the measured causative agents in SS. For example,
#'     in PERCH study, nasopharyngeal polymerase chain reaction (NPPCR; bronze-standard) may target 30 distinct pathogens, but blood culture (BCX; silver-standard) may only target a subset of the 30,
#'     so we have to specify `NA` in `ThetaSS`for those pathogens not targeted by BCX).
#'     \item{`Nu`} the number of control subjects
#'     \item{`Nd`} the number of case subjects
#'  }
#' 
#' @return A list of diagnostic test measurements, true latent statues:
#' \itemize{
#'  \item{`data_nplcm`} a list of structured data (see [nplcm()] for
#'  description). 
#'  \item{`template`} a matrix: rows for causes (may comprise a single or multiple causative agents), 
#'  columns for measurements; generated as a lookup table to match disease-class specific 
#'  parameters (true and false positive rates)
#'  \item{`latent_cat`} integer values to indicate the latent category. The integer
#'  code corresponds to the order specified in `set_parameter$etiology`.
#'  Controls are coded as `length(set_parameter$etiology)+1`.)
#'  }
#'
#' @seealso [simulate_latent] for simulating discrete latent status, given
#' which [simulate_brs] simulates bronze-standard data.
#'
#' @examples
#' K.true  <- 2   # no. of latent subclasses in actual simulation. 
#'                # If eta = c(1,0), effectively, it is K.true=1.
#' J       <- 21   # no. of pathogens.
#' N       <- 600 # no. of cases/controls.
#' 
#' eta <- c(1,0) 
#' # if it is c(1,0),then it is conditional independence model, and
#' # only the first column of parameters in PsiBS, ThetaBS matter!
#' 
#' seed_start <- 20150202
#' print(eta)
#' 
#' # set fixed simulation sequence:
#' set.seed(seed_start)
#' 
#' ThetaBS_withNA <- c(.75,rep(c(.75,.75,.75,NA),5))
#' PsiBS_withNA <- c(.15,rep(c(.05,.05,.05,NA),5))
#' 
#' ThetaSS_withNA <- c(NA,rep(c(0.15,NA,0.15,0.15),5))
#' PsiSS_withNA <- c(NA,rep(c(0,NA,0,0),5))
#' 
#' set_parameter <- list(
#'   cause_list      = c(LETTERS[1:J]),
#'   etiology        = c(c(0.36,0.1,0.1,0.1,0.1,0.05,0.05,0.05,
#'                  0.05,0.01,0.01,0.01,0.01),rep(0.00,8)), 
#'                  #same length as cause_list.
#'   pathogen_BrS    = LETTERS[1:J][!is.na(ThetaBS_withNA)],
#'   pathogen_SS     = LETTERS[1:J][!is.na(ThetaSS_withNA)],
#'   meas_nm         = list(MBS = c("MBS1"),MSS="MSS1"),
#'   Lambda          = eta, #ctrl mix
#'   Eta             = t(replicate(J,eta)), #case mix, row number equal to Jcause.
#'   PsiBS           = cbind(PsiBS_withNA[!is.na(PsiBS_withNA)],
#'                           rep(0,sum(!is.na(PsiBS_withNA)))),
#'   ThetaBS         = cbind(ThetaBS_withNA[!is.na(ThetaBS_withNA)],
#'                           rep(0,sum(!is.na(ThetaBS_withNA)))),
#'   PsiSS           = PsiSS_withNA[!is.na(PsiSS_withNA)],
#'   ThetaSS         = ThetaSS_withNA[!is.na(ThetaSS_withNA)],
#'   Nu      =     N, # control size.
#'   Nd      =     N  # case size.
#' )
#'  simu_out <- simulate_nplcm(set_parameter)
#'  data_nplcm <- simu_out$data_nplcm
#'  
#'  pathogen_display <- rev(set_parameter$pathogen_BrS)
#'  plot_logORmat(data_nplcm,pathogen_display)
#'  # more examples are provided in the vignette, including settings with 
#'  # covariates.
#' @family simulation functions
#' @export
simulate_nplcm <- function(set_parameter) {
  # simulate latent status
  latent <- simulate_latent(set_parameter)
  # simulate BrS measurements:
  out_brs   <- simulate_brs(set_parameter,latent)
  # organize  bronze-standard data:
  MBS_list <-
    list(out_brs$datres[,-grep("case",colnames(out_brs$datres)),drop = FALSE])
  names(MBS_list) <- set_parameter$meas_nm$MBS
  Mobs <- list(MBS = MBS_list, MSS=NULL, MGS = NULL)
  
  if (!is.null(set_parameter$meas_nm$MSS)){
      # simulate SS measurements:
      out_ss    <- simulate_ss(set_parameter,latent)
      # silver-standard data:
      MSS_list <-
        list(out_ss$datres[,-grep("case",colnames(out_ss$datres)),drop = FALSE])
      names(MSS_list) <- set_parameter$meas_nm$MSS
      
      Mobs <- list(MBS = MBS_list, MSS = MSS_list, MGS = NULL)
  }
  
  Y    <- out_brs$datres$case
  X    <- NULL
  
  data_nplcm <- make_list(Mobs, Y, X)
  #template   <- out_brs$template
  latent_cat <- latent$iLcat
  make_list(data_nplcm,latent_cat)
  
}

#' Simulate Latent Status:
#' 
#' @inheritParams simulate_nplcm
#'
#' @return a list of latent status samples for use in simulating measurements. It
#' also includes a template to look up measurement parameters for each disease class.
#' 
#' @family internal simulation functions
simulate_latent <- function(set_parameter) {
  # etiology common to all measurements:
  cause_list <- set_parameter$cause_list
  etiology   <- set_parameter$etiology
  Jcause     <- length(cause_list)
  
  # sample size:
  Nd       <- set_parameter$Nd
  Nu       <- set_parameter$Nu
  
  # simulate latent status (common to all measurements):
  iLcat <- rep(NA,Nd)
  #iLall <- matrix(NA,nrow = Nd + Nu,ncol = Jcause)
  etiologyMat <- matrix(NA,nrow = Nd,ncol = Jcause)
  
  # sample cause for cases:
  for (i in 1:Nd) {
    etiologyMat[i,] <- etiology
    iLcat[i]        <- sample(1:length(cause_list),1,prob = etiologyMat[i,])
  }
  
  iLnm <- cause_list[iLcat]
  # pathogen_BrS    <- set_parameter$pathogen_BrS
  # J_BrS  <- length(pathogen_BrS)
  
  # convert categorical to template (cases):
  #iL    <- symb2I(iLcat,cause_list)
  # convert back to categorical (cases):
  #iLcat.case.numeric <- Imat2cat(iL,cause_list,cause_list)
  # create
  #iLall <- rbind(iL,matrix(0,nrow = Nu,ncol = Jcause))
  #iLcatAllnumeric    <- c(iLcat.case.numeric,rep(Jcause + 1,Nu))
  
  #make_list(iLall,iLcatAllnumeric,iLcat.case.numeric,iL)
  make_list(iLcat,iLnm)
}

#' Simulate Bronze-Standard (BrS) Data
#'
#'
#' @inheritParams simulate_nplcm
#' @param latent_samples simulated latent status for all the subjects, for use in simulating
#' BrS measurements.
#'
#' @return a data frame with first column being case-control status (case at top) and
#' columns of bronze-standard measurements
#' @family internal simulation functions
simulate_brs <- function(set_parameter,latent_samples) {
  pathogen_BrS    <- set_parameter$pathogen_BrS
  cause_list      <- set_parameter$cause_list
  template    <- make_template(pathogen_BrS,cause_list)
  
  J_BrS           <- length(pathogen_BrS)
  PsiBS           <- set_parameter$PsiBS
  ThetaBS         <- set_parameter$ThetaBS
  Lambda          <- set_parameter$Lambda
  Eta             <- set_parameter$Eta
  
  iLcat <- latent_samples$iLcat
  # sample size:
  Nd       <- set_parameter$Nd
  Nu       <- set_parameter$Nu
  
  Zd <- rep(NA,Nd)
  Md <- matrix(NA,nrow = Nd,ncol = J_BrS)
  MdP <- Md
  for (i in 1:Nd) {
    Zd[i] = sample(1:ncol(Eta),1,prob = Eta[iLcat[i],])
    for (j in 1:J_BrS) {
      tmp <- template[iLcat[i],j]
      MdP[i,j]  = PsiBS[j,Zd[i]] * (1 - tmp) + tmp * ThetaBS[j,Zd[i]]
    }
  }
  Md <- rvbern(MdP)
  
  Zu  <- rep(NA,Nu)
  Mu  <- matrix(NA,nrow = Nu,ncol = J_BrS)
  MuP <- matrix(NA,nrow = Nu,ncol = J_BrS)
  for (i in 1:Nu) {
    Zu[i]     <- sample(1:length(Lambda),1,prob = Lambda)
    for (j in 1:J_BrS) {
      MuP[i,j]  <- PsiBS[j,Zu[i]]
    }
  }
  Mu <- rvbern(MuP)
  
  ## organize case/control status, iL, BS, GS data into dataframes
  datacolnames    <- c("case", pathogen_BrS)
  #   datres <- data.frame(Y = c(rep(1,Nd),rep(0,Nu)),
  #                        iLcat = iLcatAllnumeric,
  #                        iL = iLall,
  #                        MBS = rbind(Md,Mu))
  #  colnames(datres) <- datacolnames
  
  # dat_meas <- datres[,rev(rev(1:ncol(datres))[1:J_BrS]),drop=FALSE]
  #dat_case <- as.matrix(dat_meas[(1:set_parameter$Nd),])
  #dat_ctrl <- as.matrix(dat_meas[-(1:set_parameter$Nd),])
  
  # template <- make_template(pathogen_BrS, cause_list)
  
  datres <-
    data.frame(case = c(rep(1,Nd),rep(0,Nu)), MBS = rbind(Md,Mu))
  colnames(datres) <- datacolnames
  make_list(datres,template)
}

#' Simulate Silver-Standard (SS) Data
#'
#' @inheritParams simulate_nplcm
#' @param latent_samples simulated latent status for all the subjects, 
#' for use in simulating SS measurements.
#'
#' @return a data frame with first column being case-control status (case at top) and
#' columns of silver-standard measurements
#' @family internal simulation functions
simulate_ss <- function(set_parameter,latent_samples) {
  pathogen_SS    <- set_parameter$pathogen_SS
  cause_list      <- set_parameter$cause_list
  template   <- make_template(pathogen_SS,cause_list)
  J_SS           <- length(pathogen_SS)
  PsiSS           <- set_parameter$PsiSS
  ThetaSS         <- set_parameter$ThetaSS
  
  iLcat <- latent_samples$iLcat
  
  # sample size:
  Nd       <- set_parameter$Nd
  Nu       <- set_parameter$Nu
  
  Md <- matrix(NA,nrow = Nd,ncol = J_SS)
  MdP <- Md
  for (i in 1:Nd) {
    for (j in 1:J_SS) {
      tmp <- template[iLcat[i],j]
      MdP[i,j]  = PsiSS[j] * (1 - tmp) + tmp * ThetaSS[j]
    }
  }
  Md <- rvbern(MdP)
  
  Mu  <- matrix(NA,nrow = Nu,ncol = J_SS)
  MuP <- matrix(NA,nrow = Nu,ncol = J_SS)
  for (i in 1:Nu) {
    for (j in 1:J_SS) {
      MuP[i,j]  <- PsiSS[j]
    }
  }
  Mu <- rvbern(MuP)
  
  ## organize case/control status, iL, BS, GS data into dataframes
  datacolnames    <- c("case", pathogen_SS)
  #   datres <- data.frame(Y = c(rep(1,Nd),rep(0,Nu)),  n 
  #                        iLcat = iLcatAllnumeric,
  #                        iL = iLall,
  #                        MBS = rbind(Md,Mu))
  #  colnames(datres) <- datacolnames
  
  # dat_meas <- datres[,rev(rev(1:ncol(datres))[1:J_BrS]),drop=FALSE]
  #dat_case <- as.matrix(dat_meas[(1:set_parameter$Nd),])
  #dat_ctrl <- as.matrix(dat_meas[-(1:set_parameter$Nd),])
  
  #template <- make_template(pathogen_SS, cause_list)
  
  datres <-
    data.frame(case = c(rep(1,Nd),rep(0,Nu)), MSS = rbind(Md,Mu))
  colnames(datres) <- datacolnames
  make_list(datres,template)
}

