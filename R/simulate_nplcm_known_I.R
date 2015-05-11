#' Simulation case data from nested partially-latent class model (npLCM) family 
#' with every individual from the same known disease class
#'
#'
#' @param set_parameter True model parameters in the npLCM specification
#' @param I_known The name of cause that must be in \code{set_parameter$cause_list}
#'
#' @export

simulate_nplcm_known_I <- function(set_parameter,I_known){
  
  pathogen_BrS    <- set_parameter$pathogen_BrS
  cause_list      <- set_parameter$cause_list
  PsiBS    <- set_parameter$PsiBS
  ThetaBS  <- set_parameter$ThetaBS
  Nd       <- set_parameter$Nd
  Nu       <- set_parameter$Nu
  Lambda   <- set_parameter$Lambda
  Eta      <- set_parameter$Eta
  etiology <- set_parameter$etiology
  
  J_BrS  <- length(pathogen_BrS)
  Jcause <- length(cause_list)
  
  iLcat <- rep(NA,Nd)
  iLall <- matrix(NA,nrow=Nd+Nu,ncol=J_BrS)
  etiologyMat <- matrix(NA,nrow=Nd,ncol=Jcause)
  for (i in 1:Nd){
    etiologyMat[i,] <- etiology
    iLcat[i]        <- I_known #sample(cause_list,1,prob = etiologyMat[i,])
  }
  
  iL    <- symb2I(iLcat,pathogen_BrS)
  iLall <- rbind(iL,matrix(0,nrow=Nu,ncol=J_BrS))
  iLcat.case.numeric <- Imat2cat(iL,cause_list,pathogen_BrS)
  iLcatAllnumeric    <- c(iLcat.case.numeric,rep(Jcause+1,Nu))
  
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
  datacolnames    <- c("Y","iLcat",
                       paste("iL",pathogen_BrS,sep="_"),
                       paste("MBS",pathogen_BrS,sep="_"))
  datres <- data.frame(Y = c(rep(1,Nd),rep(0,Nu)),
                       iLcat = iLcatAllnumeric,
                       iL = iLall,
                       MBS = rbind(Md,Mu))
  colnames(datres) <- datacolnames
  
  template <- as.matrix(rbind(symb2I(cause_list,pathogen_BrS),rep(0,J_BrS)))
  colnames(template) <- pathogen_BrS
  rownames(template) <- c(cause_list,"control")
  return(list(template = template,
              dat     = datres[1:Nd,]))
}
