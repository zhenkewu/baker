#' insert bronze-standard measurements code chunks into .bug file
#' 
#' @param Mobs measurement data in the form of \code{data_nplcm}
#' @param cause_list a list of latent status (crucial for building templates)
#' @param meas "BrS", or "SS"
#' 
#' @return a long character string to be inserted into target .bug model file
#' @export
insert_bugfile_chunk_noreg_meas <- function(Mobs,cause_list,meas = "BrS"){
  
  
  if ("BrS" %in% meas ){
        #
        # 1. BrS data:
        #
        
        # mapping template (by `make_template` function):
        patho_BrS_list <- lapply(Mobs$MBS,colnames)
        template_BrS_list <- lapply(patho_BrS_list,make_template,cause_list) # key.
        
        #
        # create variable names:
        #
        BrS_nm   <- names(Mobs$MBS)
        
        # index measurement slices by numbers:
        JBrS_nm  <- paste("JBrS",seq_along(BrS_nm),sep = "_")
        MBS_nm  <- paste("MBS",seq_along(BrS_nm),sep = "_")
        mu_bs_nm   <- paste("mu_bs",seq_along(BrS_nm),sep = "_")
        thetaBS_nm <- paste("thetaBS",seq_along(BrS_nm),sep = "_")
        psiBS.cut_nm <- paste("psiBS.cut",seq_along(BrS_nm),sep = "_")
        psiBS_nm   <-  paste("psiBS",seq_along(BrS_nm),sep = "_")
        
        alphaB_nm     <- paste("alphaB",seq_along(BrS_nm),sep = "_")
        betaB_nm     <- paste("betaB",seq_along(BrS_nm),sep = "_")
        templateBS_nm  <- paste("templateBS",seq_along(BrS_nm),sep = "_")
        indBS_nm  <- paste("indBS",seq_along(BrS_nm),sep = "_")
        
        # only for this function to determine how many BrS measurements are there:
        nslice_BrS <- length(BrS_nm) 
        
        #
        # create chunk:
        #
        add_meas_BrS_case <- function(nslice){
          res <- rep(NA,nslice)
          for (s in 1:nslice){
            if (length(patho_BrS_list[[s]]) > 1 ){
              res[s]<-
                paste0("
                for (j in 1:",JBrS_nm[s],"){
                  ",indBS_nm[s],"[i,j] <- equals(1,",templateBS_nm[s],"[Icat[i],j])
                  ",MBS_nm[s],"[i,j] ~ dbern(",mu_bs_nm[s],"[i,j])
                  ",mu_bs_nm[s],"[i,j]<-", indBS_nm[s],"[i,j]*",thetaBS_nm[s],"[j]+(1-", indBS_nm[s],"[i,j])*",psiBS.cut_nm[s],"[j]
                }","\n")
            } else{
              res[s]<-
                paste0("\n",indBS_nm[s],"[i] <- equals(1,",templateBS_nm[s],"[Icat[i]])
                  ",MBS_nm[s],"[i] ~ dbern(",mu_bs_nm[s],"[i])
                  ",mu_bs_nm[s],"[i]<-", indBS_nm[s],"[i]*",thetaBS_nm[s],"+(1-", indBS_nm[s],"[i])*",psiBS.cut_nm[s],"\n")      
              }
          }
          
          paste0(res,collapse="")
        }  
        
        add_meas_BrS_ctrl <- function(nslice){
          res <- rep(NA,nslice)
          for (s in 1:nslice){
            if (length(patho_BrS_list[[s]]) > 1 ){
              res[s]<-paste0("
                for (j in 1:",JBrS_nm[s],"){
                  ",MBS_nm[s],"[i,j] ~ dbern(",mu_bs_nm[s],"[i,j])
                  ",mu_bs_nm[s],"[i,j]<- ",psiBS_nm[s],"[j]
                }
               ")
            } else{
                res[s]<-paste0("\n",MBS_nm[s],"[i] ~ dbern(",mu_bs_nm[s],"[i])
                ",mu_bs_nm[s],"[i]<- ",psiBS_nm[s],"\n")
              }
          }
          
          paste0(res,collapse="")
        } 
        
        add_meas_BrS_param <- function(nslice){
            res <- rep(NA,nslice)
            for (s in 1:nslice){
              if (length(patho_BrS_list[[s]]) > 1 ){
              res[s] <- paste0("
                for (j in 1:",JBrS_nm[s],"){
                  ",thetaBS_nm[s],"[j]~dbeta(",alphaB_nm[s],"[j],",betaB_nm[s],"[j])
                  ",psiBS_nm[s],"[j]~dbeta(1,1)
                  ",psiBS.cut_nm[s],"[j]<-cut(",psiBS_nm[s],"[j])
                }"
                )
              } else{
                res[s] <- paste0("\n",thetaBS_nm[s],"~dbeta(",alphaB_nm[s],",",betaB_nm[s],")
                                 ",psiBS_nm[s],"~dbeta(1,1)
                                 ",psiBS.cut_nm[s],"<-cut(",psiBS_nm[s],")\n")
              }
            }
            
            paste0(res,collapse="")
        }
  }
  
  if ("SS" %in% meas ){
      #
      # 2. SS data:
      #
      
      # mapping template (by `make_template` function):
      patho_SS_list <- lapply(Mobs$MSS,colnames)
      template_SS_list <- lapply(patho_SS_list,make_template,cause_list) # key.
      
      #
      # create variable names:
      #
      SS_nm    <- names(Mobs$MSS)
      
      # index measurement slices by numbers:
      JSS_nm  <- paste("JSS",seq_along(SS_nm),sep = "_")
      MSS_nm  <- paste("MSS",seq_along(SS_nm),sep = "_")
      mu_ss_nm   <- paste("mu_ss",seq_along(SS_nm),sep = "_")
      thetaSS_nm <- paste("thetaSS",seq_along(SS_nm),sep = "_")
      psiSS_nm   <-  paste("psiSS",seq_along(SS_nm),sep = "_")
      
      alphaS_nm     <- paste("alphaS",seq_along(SS_nm),sep = "_")
      betaS_nm     <- paste("betaS",seq_along(SS_nm),sep = "_")
      templateSS_nm  <- paste("templateSS",seq_along(SS_nm),sep = "_")
      indSS_nm  <- paste("indSS",seq_along(SS_nm),sep = "_")
      
      # only for this function to determine how many BrS measurements are there:
      nslice_SS <- length(SS_nm) 
      
      
      #
      # create chunk:
      #
      add_meas_SS_case <- function(nslice){
        res <- rep(NA,nslice)
        for (s in 1:nslice){
          if (length(patho_SS_list[[s]]) > 1 ){
            res[s]<-
              paste0("
                     for (j in 1:",JSS_nm[s],"){
                     ",indSS_nm[s],"[i,j] <- equals(1,",templateSS_nm[s],"[Icat[i],j])
                     ",MSS_nm[s],"[i,j] ~ dbern(",mu_ss_nm[s],"[i,j])
                     ",mu_ss_nm[s],"[i,j]<-", indSS_nm[s],"[i,j]*",thetaSS_nm[s],"[j]+(1-", indSS_nm[s],"[i,j])*",psiSS_nm[s],"[j]
                     }","\n")
          } else{
            res[s]<-
              paste0("\n",indSS_nm[s],"[i] <- equals(1,",templateSS_nm[s],"[Icat[i]])
                     ",MSS_nm[s],"[i] ~ dbern(",mu_ss_nm[s],"[i])
                     ",mu_ss_nm[s],"[i]<-", indSS_nm[s],"[i]*",thetaSS_nm[s],"+(1-", indSS_nm[s],"[i])*",psiSS_nm[s],"\n")      
        }
      }
      
      paste0(res,collapse="")
      }  
    
      add_meas_SS_param <- function(nslice){
        res <- rep(NA,nslice)
        for (s in 1:nslice){
          if (length(patho_SS_list[[s]]) > 1 ){
            res[s] <- paste0("
                             for (j in 1:",JSS_nm[s],"){
                             ",thetaSS_nm[s],"[j]~dbeta(",alphaS_nm[s],"[j],",betaS_nm[s],"[j])
                             ",psiSS_nm[s],"[j]<-0
                             ","
                             }"
            )
          } else{
            res[s] <- paste0("\n",thetaSS_nm[s],"~dbeta(",alphaS_nm[s],",",betaS_nm[s],")
                             ",psiSS_nm[s],"<- 0
                             ","\n")
          }
        }
        
        paste0(res,collapse="")
      }
  }
  # generate file:
  if ("BrS" %in% meas & !("SS" %in% meas)){
    chunk <- paste0(
      "# BrS measurements:
      for (i in 1:Nd){
      ",add_meas_BrS_case(nslice_BrS),"
      }
      for (i in (Nd+1):(Nd+Nu)){
      ",add_meas_BrS_ctrl(nslice_BrS),"
      }
      
      # bronze-standard measurement characteristics:
      ",add_meas_BrS_param(nslice_BrS)
    )
  }
  
  if (!("BrS" %in% meas) & ("SS" %in% meas)){
    chunk <- paste0(
      "# BrS measurements:
      for (i in 1:Nd){
      ",add_meas_SS_case(nslice_SS),"
      }
      
      # bronze-standard measurement characteristics:
      ",add_meas_SS_param(nslice_SS)
    )
  }
  
  if ("BrS" %in% meas & ("SS" %in% meas)){
    chunk <- paste0(
      "# BrS measurements:
      for (i in 1:Nd){
      ",add_meas_BrS_case(nslice_BrS),"
      ",add_meas_SS_case(nslice_SS),"
      }
      for (i in (Nd+1):(Nd+Nu)){
      ",add_meas_BrS_ctrl(nslice_BrS),"
      }
      
      # bronze-standard measurement characteristics:
      ",add_meas_BrS_param(nslice_BrS),
      add_meas_SS_param(nslice_SS)
    )
  }
  
  paste0(chunk,"\n")
} 


















#' insert etiology code chunks into .bug file
#' 
#' @return a long character string to be inserted into target .bug model file
#' @export

insert_bugfile_chunk_noreg_etiology <- function(){

  chunk_etiology <- paste0("
  # priors
  for (i in 1:Nd){
    Icat[i] ~ dcat(pEti[1:Jcause])
    
  }
  pEti[1:Jcause]~ddirch(alpha[])")
  
  paste0(chunk_etiology,"\n")
}



  
