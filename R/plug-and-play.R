#
# 1. Bronze-standard data: conditional independence model (no nested structure).
#

#' add a likelihood component for a BrS measurement slice among cases (conditional independence)
#' 
#' @param s the slice
#' @param Mobs See \code{data_nplcm} described in \code{\link{nplcm}}
#' @param cause_list the list of causes in \code{data_nplcm} described in \code{\link{nplcm}}
#' @param ppd Default is NULL; Set to TRUE for enabling posterior predictive checking.
#' @return a list of two elements: the first is \code{plug}, the .bug code; 
#' the second is \code{parameters} that stores model parameters introduced by this 
#' plugged measurement slice
#' @family likelihood specification functions
#' @family plug-and-play functions
#' @export
add_meas_BrS_case_NoNest_Slice <- function(s,Mobs,cause_list,ppd=NULL) {
  # mapping template (by `make_template` function):
  patho_BrS_list <- lapply(Mobs$MBS,colnames)
  template_BrS_list <-
    lapply(patho_BrS_list,make_template,cause_list) # key.
  
  # create variable names:
  BrS_nm   <- names(Mobs$MBS)
  # index measurement slices by numbers:
  JBrS_nm  <- paste("JBrS",seq_along(BrS_nm),sep = "_") #
  MBS_nm   <- paste("MBS",seq_along(BrS_nm),sep = "_")#
  mu_bs_nm   <- paste("mu_bs",seq_along(BrS_nm),sep = "_")#
  thetaBS_nm <- paste("thetaBS",seq_along(BrS_nm),sep = "_")#
  psiBS.cut_nm <- paste("psiBS.cut",seq_along(BrS_nm),sep = "_")#
  templateBS_nm  <- paste("templateBS",seq_along(BrS_nm),sep = "_")#
  indBS_nm  <- paste("indBS",seq_along(BrS_nm),sep = "_")#
  Icat_nm   <- "Icat"
  
  if (length(patho_BrS_list[[s]]) > 1) {
    plug <-
      paste0(
        "
          # case BrS measurement; non-nested:
          for (j in 1:",JBrS_nm[s],"){
          ",indBS_nm[s],"[i,j] <- equals(1,",templateBS_nm[s],"[",Icat_nm,"[i],j])
          ",MBS_nm[s],"[i,j] ~ dbern(",mu_bs_nm[s],"[i,j])
          ",mu_bs_nm[s],"[i,j]<-", indBS_nm[s],"[i,j]*",thetaBS_nm[s],"[j]+(1-", indBS_nm[s],"[i,j])*",psiBS.cut_nm[s],"[j]
          }","\n"
      )
  } else{
    plug <-
      paste0(
        "
           
          # case BrS measurement; non-nested (with only one column):
          ",indBS_nm[s],"[i] <- equals(1,",templateBS_nm[s],"[",Icat_nm,"[i]])
          ",MBS_nm[s],"[i] ~ dbern(",mu_bs_nm[s],"[i])
          ",mu_bs_nm[s],"[i]<-", indBS_nm[s],"[i]*",thetaBS_nm[s],"+(1-", indBS_nm[s],"[i])*",psiBS.cut_nm[s],"\n"
      )
  }
  parameters <- c(Icat_nm, thetaBS_nm[s],psiBS.cut_nm[s])
  # if posterior predictive distribution is requested:
  if (!is.null(ppd) && ppd){
    MBS_nm.new     <- paste("MBS.new",seq_along(BrS_nm),sep = "_")#
    mu_bs_nm.new   <- paste("mu_bs.new",seq_along(BrS_nm),sep = "_")#
    indBS_nm.new   <- paste("indBS.new",seq_along(BrS_nm),sep = "_")#
    Icat_nm.new    <- "Icat.new"  
    
    if (length(patho_BrS_list[[s]]) > 1) {
      plug <-
        paste0(
          "
          # case BrS measurement; non-nested:
          for (j in 1:",JBrS_nm[s],"){
          ",indBS_nm[s],"[i,j] <- equals(1,",templateBS_nm[s],"[",Icat_nm,"[i],j])
          ",MBS_nm[s],"[i,j] ~ dbern(",mu_bs_nm[s],"[i,j])
          ",mu_bs_nm[s],"[i,j]<-", indBS_nm[s],"[i,j]*",thetaBS_nm[s],"[j]+(1-", indBS_nm[s],"[i,j])*",psiBS.cut_nm[s],"[j]
          # posterior predictive distribution:
          ",indBS_nm.new[s],"[i,j] <- equals(1,",templateBS_nm[s],"[",Icat_nm.new,"[i],j])
          ",MBS_nm.new[s],"[i,j] ~ dbern(",mu_bs_nm.new[s],"[i,j])
          ",mu_bs_nm.new[s],"[i,j]<-", indBS_nm.new[s],"[i,j]*",thetaBS_nm[s],"[j]+(1-", indBS_nm.new[s],"[i,j])*",psiBS.cut_nm[s],"[j]
          }","\n"
        )
    } else{
      plug <-
        paste0(
          "
           
          # case BrS measurement; non-nested (with only one column):
          ",indBS_nm[s],"[i] <- equals(1,",templateBS_nm[s],"[",Icat_nm,"[i]])
          ",MBS_nm[s],"[i] ~ dbern(",mu_bs_nm[s],"[i])
          ",mu_bs_nm[s],"[i]<-", indBS_nm[s],"[i]*",thetaBS_nm[s],"+(1-", indBS_nm[s],"[i])*",psiBS.cut_nm[s],"
          # posterior predictive distribution:
          ",indBS_nm.new[s],"[i] <- equals(1,",templateBS_nm[s],"[",Icat_nm.new,"[i]])
          ",MBS_nm.new[s],"[i] ~ dbern(",mu_bs_nm.new[s],"[i])
          ",mu_bs_nm.new[s],"[i]<-", indBS_nm.new[s],"[i]*",thetaBS_nm[s],"+(1-", indBS_nm.new[s],"[i])*",psiBS.cut_nm[s],"
          
          
          \n"
        )
    }
    parameters <- c(Icat_nm,Icat_nm.new, thetaBS_nm[s],psiBS.cut_nm[s])
  }
  
  make_list(plug,parameters)
}


#' add a likelihood component for a BrS measurement slice among controls (conditional independence)
#' 
#' 
#' @inheritParams add_meas_BrS_case_NoNest_Slice
#' 
#' @return a list of two elements: the first is \code{plug}, the .bug code; the second is \code{parameters}
#' that stores model parameters introduced by this plugged measurement slice
#' @family likelihood specification functions
#' @family plug-and-play functions
#' @export
add_meas_BrS_ctrl_NoNest_Slice <- function(s, Mobs,cause_list,ppd=NULL) {
  # mapping template (by `make_template` function):
  patho_BrS_list <- lapply(Mobs$MBS,colnames)
  template_BrS_list <-
    lapply(patho_BrS_list,make_template,cause_list) # key.
  
  # create variable names:
  BrS_nm   <- names(Mobs$MBS)
  # index measurement slices by numbers:
  JBrS_nm  <- paste("JBrS",seq_along(BrS_nm),sep = "_")#
  MBS_nm   <- paste("MBS",seq_along(BrS_nm),sep = "_")#
  mu_bs_nm   <- paste("mu_bs",seq_along(BrS_nm),sep = "_")#
  psiBS_nm   <-  paste("psiBS",seq_along(BrS_nm),sep = "_")#
  
  if (length(patho_BrS_list[[s]]) > 1) {
    plug <- paste0(
      "       
              ## control BrS measurements; no subclass:
              for (j in 1:",JBrS_nm[s],"){
                ",MBS_nm[s],"[i,j] ~ dbern(",mu_bs_nm[s],"[i,j])
                ",mu_bs_nm[s],"[i,j]<- ",psiBS_nm[s],"[j]
              }
              "
    )
  } else{
    plug <- paste0(
      "
            ## control BrS measurements; no subclass (only one column):
            ",MBS_nm[s],"[i] ~ dbern(",mu_bs_nm[s],"[i])
            ",mu_bs_nm[s],"[i]<- ",psiBS_nm[s],"\n"
    )
  }
  
  if (!is.null(ppd) && ppd){
    
    MBS_nm.new   <- paste("MBS.new",seq_along(BrS_nm),sep = "_")#
    mu_bs_nm.new   <- paste("mu_bs.new",seq_along(BrS_nm),sep = "_")#
    
    if (length(patho_BrS_list[[s]]) > 1) {
      plug <- paste0(
        "       
              ## control BrS measurements; no subclass:
              for (j in 1:",JBrS_nm[s],"){
                ",MBS_nm[s],"[i,j] ~ dbern(",mu_bs_nm[s],"[i,j])
                ",mu_bs_nm[s],"[i,j]<- ",psiBS_nm[s],"[j]
              ## posterior predictive distribution
                ",MBS_nm.new[s],"[i,j] ~ dbern(",mu_bs_nm.new[s],"[i,j])
                ",mu_bs_nm.new[s],"[i,j]<- ",psiBS_nm[s],"[j]
              }
              "
      )
    } else{
      plug <- paste0(
        "
            ## control BrS measurements; no subclass (only one column):
            ",MBS_nm[s],"[i] ~ dbern(",mu_bs_nm[s],"[i])
            ",mu_bs_nm[s],"[i]<- ",psiBS_nm[s],"
            ## posterior predictive distribution:
            ",MBS_nm.new[s],"[i] ~ dbern(",mu_bs_nm.new[s],"[i])
            ",mu_bs_nm.new[s],"[i]<- ",psiBS_nm[s],"
             \n"
      )
    }
  }
  
  parameters <- c(psiBS_nm[s])
  make_list(plug,parameters)
}


#' add parameters for a BrS measurement slice among cases and controls (conditional independence)
#' 
#' 
#' @inheritParams add_meas_BrS_case_NoNest_Slice
#' 
#' @return a list of two elements: the first is \code{plug}, the .bug code; the second is \code{parameters}
#' that stores model parameters introduced by this plugged measurement slice
#' @family likelihood specification functions
#' @family plug-and-play functions 
#' @export

add_meas_BrS_param_NoNest_Slice <- function(s,Mobs,cause_list) {
  # mapping template (by `make_template` function):
  patho_BrS_list <- lapply(Mobs$MBS,colnames)
  template_BrS_list <-
    lapply(patho_BrS_list,make_template,cause_list) # key.
  
  # create variable names:
  BrS_nm   <- names(Mobs$MBS)
  # index measurement slices by numbers:
  JBrS_nm  <- paste("JBrS",seq_along(BrS_nm),sep = "_")#
  thetaBS_nm <- paste("thetaBS",seq_along(BrS_nm),sep = "_")#
  psiBS.cut_nm <- paste("psiBS.cut",seq_along(BrS_nm),sep = "_")#
  psiBS_nm   <-  paste("psiBS",seq_along(BrS_nm),sep = "_")#
  alphaB_nm     <- paste("alphaB",seq_along(BrS_nm),sep = "_")#
  betaB_nm     <- paste("betaB",seq_along(BrS_nm),sep = "_")#
  
  if (length(patho_BrS_list[[s]]) > 1) {
    plug <- paste0(
      "
              # BrS measurement characteristics - non-nested:
              for (j in 1:",JBrS_nm[s],"){
                ",thetaBS_nm[s],"[j]~ dbeta(",alphaB_nm[s],"[j],",betaB_nm[s],"[j])
                ",psiBS_nm[s],"[j]  ~ dbeta(1,1)
                ",psiBS.cut_nm[s],"[j]<-cut(",psiBS_nm[s],"[j])
              }"
    )
  } else{
    plug <-
      paste0(
        "
              # BrS measurement characteristics - non-nested (only one column):
              ",thetaBS_nm[s],"~  dbeta(",alphaB_nm[s],",",betaB_nm[s],")
              ",psiBS_nm[s],"  ~  dbeta(1,1)
              ",psiBS.cut_nm[s],"<-cut(",psiBS_nm[s],")\n"
      )
  }
  parameters <- c(thetaBS_nm[s],psiBS_nm[s],alphaB_nm[s],betaB_nm[s])
  make_list(plug,parameters)
}

#
# 2. Bronze-standard data: conditional dependence model (nested structure).
#


#' add likelihood for a BrS measurement slice among cases (conditional dependence)
#' 
#' 
#' @inheritParams add_meas_BrS_case_NoNest_Slice
#' 
#' @return a list of two elements: the first is \code{plug}, the .bug code; the second is \code{parameters}
#' that stores model parameters introduced by this plugged measurement slice
#' @family likelihood specification functions
#' @family plug-and-play functions
#' @export
#' 
add_meas_BrS_case_Nest_Slice <- function(s,Mobs,cause_list,ppd=NULL){
  # mapping template (by `make_template` function):
  patho_BrS_list <- lapply(Mobs$MBS,colnames)
  template_BrS_list <-
    lapply(patho_BrS_list,make_template,cause_list) # key.
  
  # create variable names:
  BrS_nm   <- names(Mobs$MBS)
  # index measurement slices by numbers:
  JBrS_nm  <- paste("JBrS",seq_along(BrS_nm),sep = "_")#
  MBS_nm   <- paste("MBS",seq_along(BrS_nm),sep = "_")#
  mu_bs.bound_nm   <- paste("mu_bs.bound",seq_along(BrS_nm),sep = "_")#
  mu_bs_nm   <- paste("mu_bs",seq_along(BrS_nm),sep = "_")#
  PR_BS_nm   <- paste("PR_BS",seq_along(BrS_nm),sep = "_")#
  PsiBS.cut_nm   <- paste("PsiBS.cut",seq_along(BrS_nm),sep = "_")#
  ThetaBS_nm <- paste("ThetaBS",seq_along(BrS_nm),sep="_")#
  Z_nm <- paste("Z",seq_along(BrS_nm),sep="_")#
  K_nm <- paste("K",seq_along(BrS_nm),sep="_")#
  
  templateBS_nm  <-
    paste("templateBS",seq_along(BrS_nm),sep = "_")#
  indBS_nm  <- paste("indBS",seq_along(BrS_nm),sep = "_")#
  Icat_nm   <- "Icat"
  
  if (length(patho_BrS_list[[s]]) == 1){stop("==cannot do nested modeling for BrS measurement with 1 dimension!==")} 
  
  plug <-
    paste0(
      "
          ## case BrS measurements; with subclasses:
          for (j in 1:",JBrS_nm[s],"){
            ",indBS_nm[s],"[i,j] <- equals(1,",templateBS_nm[s],"[",Icat_nm,"[i],j])
            ",MBS_nm[s],"[i,j]~dbern(",mu_bs.bound_nm[s],"[i,j])
            ",mu_bs.bound_nm[s],"[i,j]<-max(0.000001,min(0.999999,",mu_bs_nm[s],"[i,j]))
            ",mu_bs_nm[s],"[i,j]<-",PR_BS_nm[s],"[i,j,",Z_nm[s],"[i]]
            
            for (s in 1:",K_nm[s],"){
              ",PR_BS_nm[s],"[i,j,s]<-",PsiBS.cut_nm[s],"[j,s]*(1-",indBS_nm[s],"[i,j])+",ThetaBS_nm[s],"[j,s]*",indBS_nm[s],"[i,j]
            }
          }
      "
    )
  parameters <- c(Icat_nm,ThetaBS_nm[s])
  
  if (!is.null(ppd) && ppd){
    MBS_nm.new   <- paste("MBS.new",seq_along(BrS_nm),sep = "_")#
    mu_bs.bound_nm.new   <- paste("mu_bs.bound.new",seq_along(BrS_nm),sep = "_")#
    mu_bs_nm.new   <- paste("mu_bs.new",seq_along(BrS_nm),sep = "_")#
    PR_BS_nm.new   <- paste("PR_BS.new",seq_along(BrS_nm),sep = "_")#
    Z_nm.new <- paste("Z.new",seq_along(BrS_nm),sep="_")#
    Icat_nm.new <- "Icat.new"
    indBS_nm.new  <- paste("indBS.new",seq_along(BrS_nm),sep = "_")#
    
    plug <-
      paste0(
        "
          ## case BrS measurements; with subclasses:
          for (j in 1:",JBrS_nm[s],"){
            ",indBS_nm[s],"[i,j] <- equals(1,",templateBS_nm[s],"[",Icat_nm,"[i],j])
            ",MBS_nm[s],"[i,j]~dbern(",mu_bs.bound_nm[s],"[i,j])
            ",mu_bs.bound_nm[s],"[i,j]<-max(0.000001,min(0.999999,",mu_bs_nm[s],"[i,j]))
            ",mu_bs_nm[s],"[i,j]<-",PR_BS_nm[s],"[i,j,",Z_nm[s],"[i]]
            
            for (s in 1:",K_nm[s],"){
              ",PR_BS_nm[s],"[i,j,s]<-",PsiBS.cut_nm[s],"[j,s]*(1-",indBS_nm[s],"[i,j])+",ThetaBS_nm[s],"[j,s]*",indBS_nm[s],"[i,j]
            }
          ## posterior predictive distribution:
            ",indBS_nm.new[s],"[i,j] <- equals(1,",templateBS_nm[s],"[",Icat_nm.new,"[i],j])
            ",MBS_nm.new[s],"[i,j]~dbern(",mu_bs.bound_nm.new[s],"[i,j])
            ",mu_bs.bound_nm.new[s],"[i,j]<-max(0.000001,min(0.999999,",mu_bs_nm.new[s],"[i,j]))
            ",mu_bs_nm.new[s],"[i,j]<-",PR_BS_nm.new[s],"[i,j,",Z_nm.new[s],"[i]]
            
            for (s in 1:",K_nm[s],"){
              ",PR_BS_nm.new[s],"[i,j,s]<-",PsiBS.cut_nm[s],"[j,s]*(1-",indBS_nm.new[s],"[i,j])+",ThetaBS_nm[s],"[j,s]*",indBS_nm.new[s],"[i,j]
            }
          }
      "
      )  
    parameters <- c(Icat_nm, Icat_nm.new, ThetaBS_nm[s])
  }
  
  
  make_list(plug,parameters)
}


#' add likelihood for a BrS measurement slice among controls (conditional independence)
#' 
#' 
#' @inheritParams add_meas_BrS_case_NoNest_Slice
#' 
#' @return a list of two elements: the first is \code{plug}, the .bug code; the second is \code{parameters}
#' that stores model parameters introduced by this plugged measurement slice
#' @family likelihood specification functions
#' @family plug-and-play functions
#' @export
#' 
add_meas_BrS_ctrl_Nest_Slice <- function(s, Mobs,cause_list,ppd=NULL) {
  # mapping template (by `make_template` function):
  patho_BrS_list <- lapply(Mobs$MBS,colnames)
  template_BrS_list <-
    lapply(patho_BrS_list,make_template,cause_list) # key.
  
  # create variable names:
  BrS_nm   <- names(Mobs$MBS)
  # index measurement slices by numbers:
  JBrS_nm  <- paste("JBrS",seq_along(BrS_nm),sep = "_")#
  MBS_nm   <- paste("MBS",seq_along(BrS_nm),sep = "_")#
  mu_bs.bound_nm   <- paste("mu_bs.bound",seq_along(BrS_nm),sep = "_")#
  mu_bs_nm   <- paste("mu_bs",seq_along(BrS_nm),sep = "_")#
  Z_nm <- paste("Z",seq_along(BrS_nm),sep="_")#
  PsiBS_nm <- paste("PsiBS",seq_along(BrS_nm),sep="_")#
  
  templateBS_nm  <-
    paste("templateBS",seq_along(BrS_nm),sep = "_")
  
  if (length(patho_BrS_list[[s]]) == 1){stop("==cannot do nested modeling for BrS measurement with 1 dimension!==")} 
  plug <- paste0(
    "   
          ## control BrS measurements;  with subclasses:
          for (j in 1:",JBrS_nm[s],"){
              ",MBS_nm[s],"[i,j]~dbern(",mu_bs.bound_nm[s],"[i,j])
              ",mu_bs.bound_nm[s],"[i,j] <-max(0.000001,min(0.999999,",mu_bs_nm[s],"[i,j]))
              ",mu_bs_nm[s],"[i,j]<-",PsiBS_nm[s],"[j,",Z_nm[s],"[i]]
          }
          "
  )
  
  
  
  if (!is.null(ppd) && ppd){
    MBS_nm.new   <- paste("MBS.new",seq_along(BrS_nm),sep = "_")#
    mu_bs.bound_nm.new   <- paste("mu_bs.bound.new",seq_along(BrS_nm),sep = "_")#
    mu_bs_nm.new   <- paste("mu_bs.new",seq_along(BrS_nm),sep = "_")#
    Z_nm.new <- paste("Z.new",seq_along(BrS_nm),sep="_")#
    
    plug <- paste0(
      "   
          ## control BrS measurements;  with subclasses:
          for (j in 1:",JBrS_nm[s],"){
              ",MBS_nm[s],"[i,j]~dbern(",mu_bs.bound_nm[s],"[i,j])
              ",mu_bs.bound_nm[s],"[i,j] <-max(0.000001,min(0.999999,",mu_bs_nm[s],"[i,j]))
              ",mu_bs_nm[s],"[i,j]<-",PsiBS_nm[s],"[j,",Z_nm[s],"[i]]
              ## posterior predictive distribution:
              ",MBS_nm.new[s],"[i,j]~dbern(",mu_bs.bound_nm.new[s],"[i,j])
              ",mu_bs.bound_nm.new[s],"[i,j] <-max(0.000001,min(0.999999,",mu_bs_nm.new[s],"[i,j]))
              ",mu_bs_nm.new[s],"[i,j]<-",PsiBS_nm[s],"[j,",Z_nm.new[s],"[i]]
          }
          "
    )    
  }
  
  parameters <- c(PsiBS_nm[s])
  make_list(plug,parameters)
}

#' add parameters for a BrS measurement slice among cases and controls (conditional dependence)
#' 
#' 
#' @inheritParams add_meas_BrS_case_NoNest_Slice
#' 
#' @return a list of two elements: the first is \code{plug}, the .bug code; the second is \code{parameters}
#' that stores model parameters introduced by this plugged measurement slice
#' @family likelihood specification functions
#' @family plug-and-play functions 
#' @export
#' 


add_meas_BrS_param_Nest_Slice <- function(s,Mobs,cause_list) { #note: has separated case and controls subclass weights.
  # mapping template (by `make_template` function):
  patho_BrS_list <- lapply(Mobs$MBS,colnames)
  template_BrS_list <-
    lapply(patho_BrS_list,make_template,cause_list) # key.
  
  # create variable names:
  BrS_nm   <- names(Mobs$MBS)
  # index measurement slices by numbers:
  JBrS_nm  <- paste("JBrS",seq_along(BrS_nm),sep = "_")#
  PsiBS.cut_nm   <- paste("PsiBS.cut",seq_along(BrS_nm),sep = "_")#
  ThetaBS_nm <- paste("ThetaBS",seq_along(BrS_nm),sep="_")#
  PsiBS_nm <- paste("PsiBS",seq_along(BrS_nm),sep="_")#
  Lambda_nm <- paste("Lambda",seq_along(BrS_nm),sep="_")#
  Lambda0_nm <- paste("Lambda0",seq_along(BrS_nm),sep="_")#
  r0_nm <- paste("r0",seq_along(BrS_nm),sep="_")#
  alphadp0_nm <- paste("alphadp0",seq_along(BrS_nm),sep="_")#
  alphadp0_case_nm <- paste("alphadp0_case",seq_along(BrS_nm),sep="_")# <--- cases' subclass weights.
  
  Eta0_nm <- paste("Eta0",seq_along(BrS_nm),sep="_")#
  r1_nm <- paste("r1",seq_along(BrS_nm),sep="_")#
  Eta_nm <- paste("Eta",seq_along(BrS_nm),sep="_")#
  alphaB_nm     <- paste("alphaB",seq_along(BrS_nm),sep = "_")#
  betaB_nm     <- paste("betaB",seq_along(BrS_nm),sep = "_")#
  K_nm <- paste("K",seq_along(BrS_nm),sep="_")#
  
  templateBS_nm  <-
    paste("templateBS",seq_along(BrS_nm),sep = "_")
  indBS_nm  <- paste("indBS",seq_along(BrS_nm),sep = "_")
  
  if (length(patho_BrS_list[[s]]) == 1){stop("==cannot do nested modeling for BrS measurement with 1 dimension!==")} 
  plug <- paste0(
    "
    ## cut the feedback from case model to FPR:
    for (j in 1:",JBrS_nm[s],"){
      for (s in 1:",K_nm[s],"){
      ",PsiBS.cut_nm[s],"[j,s]<-cut(",PsiBS_nm[s],"[j,s])
      }
    }
    
    ####################################
    ### stick-breaking prior specification
    ####################################
    # control subclass mixing weights:
    ",Lambda0_nm[s],"[1]<-",r0_nm[s],"[1]
    ",r0_nm[s],"[",K_nm[s],"]<-1
    for(j in 2:",K_nm[s],") {",Lambda0_nm[s],"[j]<-",r0_nm[s],"[j]*(1-",r0_nm[s],"[j-1])*",Lambda0_nm[s],"[j-1]/",r0_nm[s],"[j-1]}
    for(k in 1:(",K_nm[s],"-1)){
      ",r0_nm[s],"[k]~dbeta(1,",alphadp0_nm[s],")I(0.000001,0.999999)
    }
    
    for (k in 1:(",K_nm[s],"-1)){",Lambda_nm[s],"[k]<-max(0.000001,min(0.999999,",Lambda0_nm[s],"[k]))}
    ",Lambda_nm[s],"[",K_nm[s],"]<-1-sum(",Lambda_nm[s],"[1:(",K_nm[s],"-1)])
    
    # case subclass mixing weights:
    ",Eta0_nm[s],"[1]<-",r1_nm[s],"[1]
    ",r1_nm[s],"[",K_nm[s],"]<-1
    for(j in 2:",K_nm[s],") {",Eta0_nm[s],"[j]<-",r1_nm[s],"[j]*(1-",r1_nm[s],"[j-1])*",Eta0_nm[s],"[j-1]/",r1_nm[s],"[j-1]}
    for(k in 1:(",K_nm[s],"-1)){
          ",r1_nm[s],"[k]~dbeta(1,",alphadp0_case_nm[s],")I(0.000001,0.999999)
    }
    
    for (k in 1:(",K_nm[s],"-1)){",Eta_nm[s],"[k]<-max(0.000001,min(0.999999,",Eta0_nm[s],"[k]))}
    ",Eta_nm[s],"[",K_nm[s],"]<-1-sum(",Eta_nm[s],"[1:(",K_nm[s],"-1)])
    
    ",alphadp0_nm[s],"~dgamma(.25,.25)I(0.001,20)
    ",alphadp0_case_nm[s],"~dgamma(.25,.25)I(0.001,20)
    
    #########################
    ## priors on TPR and FPR:
    #########################
    
    for (j in 1:",JBrS_nm[s],"){
      for (s in 1:",K_nm[s],"){
        ",PsiBS_nm[s],"[j,s]~dbeta(1,1)
        #ThetaBS[j,s]~dbeta(1,1)
        ",ThetaBS_nm[s],"[j,s]~dbeta(",alphaB_nm[s],"[j],",betaB_nm[s],"[j])
      }
    }
    "
  )
  
  parameters <- c(PsiBS_nm[s], ThetaBS_nm[s], Lambda_nm[s],Eta_nm[s],alphaB_nm[s],betaB_nm[s])
  make_list(plug,parameters)
}


#' add subclass indicators for a BrS measurement slice among cases and controls (conditional independence)
#' 
#' 
#' @inheritParams add_meas_BrS_case_NoNest_Slice
#' 
#' @return a list of two elements: the first is \code{plug}, the .bug code; the second is \code{parameters}
#' that stores model parameters introduced by this plugged measurement slice
#' @family likelihood specification functions
#' @family plug-and-play functions
#' @export

add_meas_BrS_subclass_Nest_Slice <- function(s,Mobs,cause_list,ppd=NULL){
  # mapping template (by `make_template` function):
  patho_BrS_list <- lapply(Mobs$MBS,colnames)
  template_BrS_list <-
    lapply(patho_BrS_list,make_template,cause_list) # key.
  
  # create variable names:
  BrS_nm   <- names(Mobs$MBS)
  # index measurement slices by numbers:
  Z_nm <- paste("Z",seq_along(BrS_nm),sep="_")#
  Eta_nm <- paste("Eta",seq_along(BrS_nm),sep="_")#
  Lambda_nm <- paste("Lambda",seq_along(BrS_nm),sep="_")#
  K_nm <- paste("K",seq_along(BrS_nm),sep="_")#
  
  plug <- paste0(
    " 
      # cases' subclass indicators:
      for (i in 1:Nd){
        ",Z_nm[s],"[i] ~ dcat(",Eta_nm[s],"[1:",K_nm[s],"])
      }
      # controls' subclass indicators:
      for (i in (Nd+1):(Nd+Nu)){
        ",Z_nm[s],"[i] ~ dcat(",Lambda_nm[s],"[1:",K_nm[s],"])
      }
    "
  )
  
  if (!is.null(ppd) && ppd){
    Z_nm.new <- paste("Z.new",seq_along(BrS_nm),sep="_")#
    plug <- paste0(
      " 
      # cases' subclass indicators:
      for (i in 1:Nd){
        ",Z_nm[s],"[i] ~ dcat(",Eta_nm[s],"[1:",K_nm[s],"])
        ",Z_nm.new[s],"[i] ~ dcat(",Eta_nm[s],"[1:",K_nm[s],"])
      }
      # controls' subclass indicators:
      for (i in (Nd+1):(Nd+Nu)){
        ",Z_nm[s],"[i] ~ dcat(",Lambda_nm[s],"[1:",K_nm[s],"])
        ",Z_nm.new[s],"[i] ~ dcat(",Lambda_nm[s],"[1:",K_nm[s],"])
      }
    "
    )  
  }
  
  parameters <- c(Eta_nm[s],Lambda_nm[s])
  make_list(plug,parameters)
}

#
# 3. Silver standard data:
#

#' add likelihood for a SS measurement slice among cases (conditional independence)
#' 
#' 
#' @param nslice the total number of SS measurement slices
#' @param Mobs see \code{data_nplcm} described in \code{\link{nplcm}}
#' @param prior see \code{model_options} described in \code{\link{nplcm}}
#' @param cause_list the list of causes in \code{model_options} described in \code{\link{nplcm}}
#' 
#' @return a list of two elements: the first is \code{plug}, the .bug code; the second is \code{parameters}
#' that stores model parameters introduced by this plugged measurement slice
#' @family likelihood specification functions
#' @family plug-and-play functions
#' @export

add_meas_SS_case <- function(nslice,Mobs,prior,cause_list) {
  # mapping template (by `make_template` function):
  patho_SS_list <- lapply(Mobs$MSS,colnames)
  template_SS_list <-
    lapply(patho_SS_list,make_template,cause_list) # key.
  
  # create variable names:
  SS_nm    <- names(Mobs$MSS)
  # index measurement slices by numbers:
  JSS_nm  <- paste("JSS",seq_along(SS_nm),sep = "_")#
  MSS_nm  <- paste("MSS",seq_along(SS_nm),sep = "_")#
  mu_ss_nm   <- paste("mu_ss",seq_along(SS_nm),sep = "_")#
  thetaSS_nm <- paste("thetaSS",seq_along(SS_nm),sep = "_")#
  psiSS_nm   <-  paste("psiSS",seq_along(SS_nm),sep = "_")#
  
  templateSS_nm  <-
    paste("templateSS",seq_along(SS_nm),sep = "_")#
  indSS_nm  <- paste("indSS",seq_along(SS_nm),sep = "_")#
  
  # for SS TPR across groups (currently only allows uniform grouping across dimensions, i.e.,
  # only allow the same way of splitting cases for every pathogen):
  SS_TPR_strat <- FALSE
  
  SS_TPR_grp_nm <- paste("SS_TPR_grp",seq_along(SS_nm),sep = "_")
  GSS_TPR_nm    <- paste("GSS_TPR",seq_along(SS_nm),sep = "_") # level of groups within each slice.
  prior_SS      <- prior$TPR_prior$SS
  
  if (!is.null(prior_SS$grp) && length(unique(prior_SS$grp)) >1 ){
    SS_TPR_strat <- TRUE
  }
  
  res <- rep(NA,nslice)
  for (s in 1:nslice) {
    if (!SS_TPR_strat){# without stratified TPR in SS:
      if (length(patho_SS_list[[s]]) > 1) { # if more dim>1.
        res[s] <-
          paste0(
          "
          # cases' SS measurements:
          for (j in 1:",JSS_nm[s],"){
              ",indSS_nm[s],"[i,j] <- equals(1,",templateSS_nm[s],"[Icat[i],j])
              ",MSS_nm[s],"[i,j] ~ dbern(",mu_ss_nm[s],"[i,j])
              ",mu_ss_nm[s],"[i,j]<-", indSS_nm[s],"[i,j]*",thetaSS_nm[s],"[j]+(1-", indSS_nm[s],"[i,j])*",psiSS_nm[s],"[j]
          }
              
          ")
      } else{
        res[s] <-
          paste0(
          "
          # cases' SS measurements (only one column):
          ",indSS_nm[s],"[i] <- equals(1,",templateSS_nm[s],"[Icat[i]])
          ",MSS_nm[s],"[i] ~ dbern(",mu_ss_nm[s],"[i])
          ",mu_ss_nm[s],"[i]<-", indSS_nm[s],"[i]*",thetaSS_nm[s],"+(1-",indSS_nm[s],"[i])*",psiSS_nm[s],"\n"
          )
      }
    } else{# WITH stratified TPR in SS:
      if (length(patho_SS_list[[s]]) > 1) {
        res[s] <-
          paste0(
            "
          # cases' SS measurements:
          for (j in 1:",JSS_nm[s],"){
                ",indSS_nm[s],"[i,j] <- equals(1,",templateSS_nm[s],"[Icat[i],j])
                ",MSS_nm[s],"[i,j] ~ dbern(",mu_ss_nm[s],"[i,j])
                ",mu_ss_nm[s],"[i,j]<-", indSS_nm[s],"[i,j]*",thetaSS_nm[s],"[",SS_TPR_grp_nm[s],"[i], j]+(1-", indSS_nm[s],"[i,j])*",psiSS_nm[s],"[j]
          }
          ")
      } else{
        res[s] <-
          paste0(
          "
          # cases' SS measurements (only one column):
          ",indSS_nm[s],"[i] <- equals(1,",templateSS_nm[s],"[Icat[i]])
          ",MSS_nm[s],"[i] ~ dbern(",mu_ss_nm[s],"[i])
          ",mu_ss_nm[s],"[i]<-", indSS_nm[s],"[i]*",thetaSS_nm[s],"[",SS_TPR_grp_nm[s],"[i]]+(1-", indSS_nm[s],"[i])*",psiSS_nm[s],"
          "
          )
      }
    }
  }
  
  plug <- paste0(res,collapse = "")
  
  parameters <- c("Icat",thetaSS_nm[s],psiSS_nm[s])
  make_list(plug,parameters)
}

#' add parameters for a SS measurement slice among cases (conditional independence)
#' 
#' 
#' @inheritParams add_meas_SS_case
#' 
#' @return a list of two elements: the first is \code{plug}, the .bug code; the second is \code{parameters}
#' that stores model parameters introduced by this plugged measurement slice
#' @family likelihood specification functions
#' @family plug-and-play functions
#' @export

add_meas_SS_param <- function(nslice,Mobs,prior,cause_list) {
  # mapping template (by `make_template` function):
  patho_SS_list <- lapply(Mobs$MSS,colnames)
  template_SS_list <-
    lapply(patho_SS_list,make_template,cause_list) # key.
  
  # create variable names:
  SS_nm    <- names(Mobs$MSS)
  # index measurement slices by numbers:
  JSS_nm  <- paste("JSS",seq_along(SS_nm),sep = "_")#
  thetaSS_nm <- paste("thetaSS",seq_along(SS_nm),sep = "_")#
  psiSS_nm   <-  paste("psiSS",seq_along(SS_nm),sep = "_")#
  
  alphaS_nm     <- paste("alphaS",seq_along(SS_nm),sep = "_")#
  betaS_nm     <- paste("betaS",seq_along(SS_nm),sep = "_")#
  
  # for SS TPR across groups (currently only allows uniform grouping across dimensions, i.e.,
  # only allow the same way of splitting cases for every pathogen):
  SS_TPR_strat <- FALSE
  
  SS_TPR_grp_nm <- paste("SS_TPR_grp",seq_along(SS_nm),sep = "_")
  GSS_TPR_nm    <- paste("GSS_TPR",seq_along(SS_nm),sep = "_") # level of groups within each slice.
  prior_SS      <- prior$TPR_prior$SS
  
  if (!is.null(prior_SS$grp) && length(unique(prior_SS$grp)) >1 ){
    SS_TPR_strat <- TRUE
  }
  
  res <- rep(NA,nslice)
  for (s in 1:nslice) {
    if (!SS_TPR_strat){# without stratified TPR in SS:
      if (length(patho_SS_list[[s]]) > 1) {
        res[s] <- paste0(
          "
          for (j in 1:",JSS_nm[s],"){
              ",thetaSS_nm[s],"[j]~dbeta(",alphaS_nm[s],"[j],",betaS_nm[s],"[j])
              ",psiSS_nm[s],"[j]<-0
              ","
          }"
        )
      } else{
        res[s] <-
          paste0(
          "
          ",thetaSS_nm[s],"~dbeta(",alphaS_nm[s],",",betaS_nm[s],")
          ",psiSS_nm[s],"<- 0
          "
          )
      }
    }else{# WITH stratified TPR in SS:
      if (length(patho_SS_list[[s]]) > 1) {
        res[s] <- paste0(
          "
          for (j in 1:",JSS_nm[s],"){
              for (g in 1:",GSS_TPR_nm[s],"){
                  ",thetaSS_nm[s],"[g,j]~dbeta(",alphaS_nm[s],"[g,j],",betaS_nm[s],"[g,j])
                  ","
              }
          ",psiSS_nm[s],"[j]<-0
          }"
        )
      } else{
        res[s] <-
          paste0(
          "
          for (g in 1:",GSS_TPR_nm[s],"){
              ",thetaSS_nm[s],"[g]~dbeta(",alphaS_nm[s],"[g],",betaS_nm[s],"[g])
              ","
          }
          ",psiSS_nm[s],"<- 0
          "
          )
      }
      
    }
  }
  
  plug <- paste0(res,collapse = "")
  parameters <- c(thetaSS_nm[s],psiSS_nm[s], alphaS_nm[s],betaS_nm[s])
  make_list(plug,parameters)
}




##
##
##
##
##
##
##
##
##
##
##           FOR JAGS:
##
##
##
##
##
##
##
##
##
##
##
##
##
##
##
##
##
##
##
##
##
##
##
##
##
##
##
##
##
##
##
##
##
##
##
##
##
##
##
##
##
##
##
##
##
##
##
##
##
##
##
##
##
##
##
##
##

#
# 1. Bronze-standard data: conditional independence model (no nested structure).
#

#' add a likelihood component for a BrS measurement slice among cases (conditional independence)
#' 
#' @param s the slice
#' @param Mobs See \code{data_nplcm} described in \code{\link{nplcm}}
#' @param cause_list the list of causes in \code{data_nplcm} described in \code{\link{nplcm}}
#' @param ppd Default is NULL; Set to TRUE for enabling posterior predictive checking.
#' @return a list of two elements: the first is \code{plug}, the .bug code; 
#' the second is \code{parameters} that stores model parameters introduced by this 
#' plugged measurement slice
#' @family likelihood specification functions
#' @family plug-and-play functions
#' @export
add_meas_BrS_case_NoNest_Slice_jags <- function(s,Mobs,cause_list,ppd=NULL) {
  # mapping template (by `make_template` function):
  patho_BrS_list <- lapply(Mobs$MBS,colnames)
  template_BrS_list <-
    lapply(patho_BrS_list,make_template,cause_list) # key.
  
  # create variable names:
  BrS_nm   <- names(Mobs$MBS)
  # index measurement slices by numbers:
  JBrS_nm  <- paste("JBrS",seq_along(BrS_nm),sep = "_") #
  MBS_nm   <- paste("MBS",seq_along(BrS_nm),sep = "_")#
  mu_bs_nm   <- paste("mu_bs",seq_along(BrS_nm),sep = "_")#
  thetaBS_nm <- paste("thetaBS",seq_along(BrS_nm),sep = "_")#
  psiBS_nm <- paste("psiBS",seq_along(BrS_nm),sep = "_")#
  #psiBS.cut_nm <- paste("psiBS.cut",seq_along(BrS_nm),sep = "_")#
  templateBS_nm  <- paste("templateBS",seq_along(BrS_nm),sep = "_")#
  indBS_nm  <- paste("indBS",seq_along(BrS_nm),sep = "_")#
  Icat_nm   <- "Icat"
  
  if (length(patho_BrS_list[[s]]) > 1) {
    plug <-
      paste0(
        "
        # case BrS measurement; non-nested:
        for (j in 1:",JBrS_nm[s],"){
        ",indBS_nm[s],"[i,j] <- equals(1,",templateBS_nm[s],"[",Icat_nm,"[i],j])
        ",MBS_nm[s],"[i,j] ~ dbern(",mu_bs_nm[s],"[i,j])
        ",mu_bs_nm[s],"[i,j]<-", indBS_nm[s],"[i,j]*",thetaBS_nm[s],"[j]+(1-", indBS_nm[s],"[i,j])*",psiBS_nm[s],"[j]
        }","\n"
      )
  } else{
    plug <-
      paste0(
        "
        
        # case BrS measurement; non-nested (with only one column):
        ",indBS_nm[s],"[i] <- equals(1,",templateBS_nm[s],"[",Icat_nm,"[i]])
        ",MBS_nm[s],"[i] ~ dbern(",mu_bs_nm[s],"[i])
        ",mu_bs_nm[s],"[i]<-", indBS_nm[s],"[i]*",thetaBS_nm[s],"+(1-", indBS_nm[s],"[i])*",psiBS_nm[s],"\n"
      )
  }
  parameters <- c(Icat_nm, thetaBS_nm[s],psiBS_nm[s])
  # if posterior predictive distribution is requested:
  if (!is.null(ppd) && ppd){
    MBS_nm.new     <- paste("MBS.new",seq_along(BrS_nm),sep = "_")#
    mu_bs_nm.new   <- paste("mu_bs.new",seq_along(BrS_nm),sep = "_")#
    indBS_nm.new   <- paste("indBS.new",seq_along(BrS_nm),sep = "_")#
    Icat_nm.new    <- "Icat.new"  
    
    if (length(patho_BrS_list[[s]]) > 1) {
      plug <-
        paste0(
          "
          # case BrS measurement; non-nested:
          for (j in 1:",JBrS_nm[s],"){
          ",indBS_nm[s],"[i,j] <- equals(1,",templateBS_nm[s],"[",Icat_nm,"[i],j])
          ",MBS_nm[s],"[i,j] ~ dbern(",mu_bs_nm[s],"[i,j])
          ",mu_bs_nm[s],"[i,j]<-", indBS_nm[s],"[i,j]*",thetaBS_nm[s],"[j]+(1-", indBS_nm[s],"[i,j])*",psiBS_nm[s],"[j]
          # posterior predictive distribution:
          ",indBS_nm.new[s],"[i,j] <- equals(1,",templateBS_nm[s],"[",Icat_nm.new,"[i],j])
          ",MBS_nm.new[s],"[i,j] ~ dbern(",mu_bs_nm.new[s],"[i,j])
          ",mu_bs_nm.new[s],"[i,j]<-", indBS_nm.new[s],"[i,j]*",thetaBS_nm[s],"[j]+(1-", indBS_nm.new[s],"[i,j])*",psiBS_nm[s],"[j]
          }","\n"
        )
    } else{
      plug <-
        paste0(
          "
          
          # case BrS measurement; non-nested (with only one column):
          ",indBS_nm[s],"[i] <- equals(1,",templateBS_nm[s],"[",Icat_nm,"[i]])
          ",MBS_nm[s],"[i] ~ dbern(",mu_bs_nm[s],"[i])
          ",mu_bs_nm[s],"[i]<-", indBS_nm[s],"[i]*",thetaBS_nm[s],"+(1-", indBS_nm[s],"[i])*",psiBS_nm[s],"
          # posterior predictive distribution:
          ",indBS_nm.new[s],"[i] <- equals(1,",templateBS_nm[s],"[",Icat_nm.new,"[i]])
          ",MBS_nm.new[s],"[i] ~ dbern(",mu_bs_nm.new[s],"[i])
          ",mu_bs_nm.new[s],"[i]<-", indBS_nm.new[s],"[i]*",thetaBS_nm[s],"+(1-", indBS_nm.new[s],"[i])*",psiBS_nm[s],"
          
          
          \n"
        )
    }
    parameters <- c(Icat_nm,Icat_nm.new, thetaBS_nm[s],psiBS_nm[s])
  }
  
  make_list(plug,parameters)
}

#' add parameters for a BrS measurement slice among cases and controls (conditional independence)
#' 
#' 
#' @inheritParams add_meas_BrS_case_NoNest_Slice
#' 
#' @return a list of two elements: the first is \code{plug}, the .bug code; the second is \code{parameters}
#' that stores model parameters introduced by this plugged measurement slice
#' @family likelihood specification functions
#' @family plug-and-play functions
#' @export

add_meas_BrS_param_NoNest_Slice_jags <- function(s,Mobs,cause_list) {
  # mapping template (by `make_template` function):
  patho_BrS_list <- lapply(Mobs$MBS,colnames)
  template_BrS_list <-
    lapply(patho_BrS_list,make_template,cause_list) # key.
  
  # create variable names:
  BrS_nm   <- names(Mobs$MBS)
  # index measurement slices by numbers:
  JBrS_nm  <- paste("JBrS",seq_along(BrS_nm),sep = "_")#
  thetaBS_nm <- paste("thetaBS",seq_along(BrS_nm),sep = "_")#
  #psiBS.cut_nm <- paste("psiBS.cut",seq_along(BrS_nm),sep = "_")#
  psiBS_nm   <-  paste("psiBS",seq_along(BrS_nm),sep = "_")#
  alphaB_nm     <- paste("alphaB",seq_along(BrS_nm),sep = "_")#
  betaB_nm     <- paste("betaB",seq_along(BrS_nm),sep = "_")#
  
  if (length(patho_BrS_list[[s]]) > 1) {
    plug <- paste0(
      "
      # BrS measurement characteristics - non-nested:
      for (j in 1:",JBrS_nm[s],"){
      ",thetaBS_nm[s],"[j]~ dbeta(",alphaB_nm[s],"[j],",betaB_nm[s],"[j])
      ",psiBS_nm[s],"[j]  ~ dbeta(1,1)
      }"
    )
  } else{
    plug <-
      paste0(
        "
        # BrS measurement characteristics - non-nested (only one column):
        ",thetaBS_nm[s],"~  dbeta(",alphaB_nm[s],",",betaB_nm[s],")
        ",psiBS_nm[s],"  ~  dbeta(1,1)
        \n"
      )
  }
  parameters <- c(thetaBS_nm[s],psiBS_nm[s],alphaB_nm[s],betaB_nm[s])
  make_list(plug,parameters)
}

#
# 2. Bronze-standard data: conditional dependence model (nested structure).
#


#' add likelihood for a BrS measurement slice among cases (conditional dependence)
#' 
#' 
#' @inheritParams add_meas_BrS_case_NoNest_Slice
#' 
#' @return a list of two elements: the first is \code{plug}, the .bug code; the second is \code{parameters}
#' that stores model parameters introduced by this plugged measurement slice
#' @family likelihood specification functions
#' @family plug-and-play functions
#' @export
#' 
add_meas_BrS_case_Nest_Slice_jags <- function(s,Mobs,cause_list,ppd=NULL){
  # mapping template (by `make_template` function):
  patho_BrS_list <- lapply(Mobs$MBS,colnames)
  template_BrS_list <-
    lapply(patho_BrS_list,make_template,cause_list) # key.
  
  # create variable names:
  BrS_nm   <- names(Mobs$MBS)
  # index measurement slices by numbers:
  JBrS_nm  <- paste("JBrS",seq_along(BrS_nm),sep = "_")#
  MBS_nm   <- paste("MBS",seq_along(BrS_nm),sep = "_")#
  mu_bs.bound_nm   <- paste("mu_bs.bound",seq_along(BrS_nm),sep = "_")#
  mu_bs_nm   <- paste("mu_bs",seq_along(BrS_nm),sep = "_")#
  PR_BS_nm   <- paste("PR_BS",seq_along(BrS_nm),sep = "_")#
  #PsiBS.cut_nm   <- paste("PsiBS.cut",seq_along(BrS_nm),sep = "_")#
  PsiBS_nm   <- paste("PsiBS",seq_along(BrS_nm),sep = "_")#
  ThetaBS_nm <- paste("ThetaBS",seq_along(BrS_nm),sep="_")#
  Z_nm <- paste("Z",seq_along(BrS_nm),sep="_")#
  K_nm <- paste("K",seq_along(BrS_nm),sep="_")#
  
  templateBS_nm  <-
    paste("templateBS",seq_along(BrS_nm),sep = "_")#
  indBS_nm  <- paste("indBS",seq_along(BrS_nm),sep = "_")#
  Icat_nm   <- "Icat"
  
  if (length(patho_BrS_list[[s]]) == 1){stop("==cannot do nested modeling for BrS measurement with 1 dimension!==")} 
  
  plug <-
    paste0(
      "
      ## case BrS measurements; with subclasses:
      for (j in 1:",JBrS_nm[s],"){
      ",indBS_nm[s],"[i,j] <- equals(1,",templateBS_nm[s],"[",Icat_nm,"[i],j])
      ",MBS_nm[s],"[i,j]~dbern(",mu_bs.bound_nm[s],"[i,j])
      ",mu_bs.bound_nm[s],"[i,j]<-max(0.000001,min(0.999999,",mu_bs_nm[s],"[i,j]))
      ",mu_bs_nm[s],"[i,j]<-",PR_BS_nm[s],"[i,j,",Z_nm[s],"[i]]
      
      for (s in 1:",K_nm[s],"){
      ",PR_BS_nm[s],"[i,j,s]<-",PsiBS_nm[s],"[j,s]*(1-",indBS_nm[s],"[i,j])+",ThetaBS_nm[s],"[j,s]*",indBS_nm[s],"[i,j]
      }
      }
      "
    )
  parameters <- c(Icat_nm,ThetaBS_nm[s])
  
  if (!is.null(ppd) && ppd){
    MBS_nm.new   <- paste("MBS.new",seq_along(BrS_nm),sep = "_")#
    mu_bs.bound_nm.new   <- paste("mu_bs.bound.new",seq_along(BrS_nm),sep = "_")#
    mu_bs_nm.new   <- paste("mu_bs.new",seq_along(BrS_nm),sep = "_")#
    PR_BS_nm.new   <- paste("PR_BS.new",seq_along(BrS_nm),sep = "_")#
    Z_nm.new <- paste("Z.new",seq_along(BrS_nm),sep="_")#
    Icat_nm.new <- "Icat.new"
    indBS_nm.new  <- paste("indBS.new",seq_along(BrS_nm),sep = "_")#
    
    plug <-
      paste0(
        "
        ## case BrS measurements; with subclasses:
        for (j in 1:",JBrS_nm[s],"){
        ",indBS_nm[s],"[i,j] <- equals(1,",templateBS_nm[s],"[",Icat_nm,"[i],j])
        ",MBS_nm[s],"[i,j]~dbern(",mu_bs.bound_nm[s],"[i,j])
        ",mu_bs.bound_nm[s],"[i,j]<-max(0.000001,min(0.999999,",mu_bs_nm[s],"[i,j]))
        ",mu_bs_nm[s],"[i,j]<-",PR_BS_nm[s],"[i,j,",Z_nm[s],"[i]]
        
        for (s in 1:",K_nm[s],"){
        ",PR_BS_nm[s],"[i,j,s]<-",PsiBS_nm[s],"[j,s]*(1-",indBS_nm[s],"[i,j])+",ThetaBS_nm[s],"[j,s]*",indBS_nm[s],"[i,j]
        }
        ## posterior predictive distribution:
        ",indBS_nm.new[s],"[i,j] <- equals(1,",templateBS_nm[s],"[",Icat_nm.new,"[i],j])
        ",MBS_nm.new[s],"[i,j]~dbern(",mu_bs.bound_nm.new[s],"[i,j])
        ",mu_bs.bound_nm.new[s],"[i,j]<-max(0.000001,min(0.999999,",mu_bs_nm.new[s],"[i,j]))
        ",mu_bs_nm.new[s],"[i,j]<-",PR_BS_nm.new[s],"[i,j,",Z_nm.new[s],"[i]]
        
        for (s in 1:",K_nm[s],"){
        ",PR_BS_nm.new[s],"[i,j,s]<-",PsiBS_nm[s],"[j,s]*(1-",indBS_nm.new[s],"[i,j])+",ThetaBS_nm[s],"[j,s]*",indBS_nm.new[s],"[i,j]
        }
        }
        "
      )  
    parameters <- c(Icat_nm, Icat_nm.new, ThetaBS_nm[s])
  }
  
  
  make_list(plug,parameters)
}


#' add parameters for a BrS measurement slice among cases and controls (conditional dependence)
#' 
#' 
#' @inheritParams add_meas_BrS_case_NoNest_Slice
#' 
#' @return a list of two elements: the first is \code{plug}, the .bug code; the second is \code{parameters}
#' that stores model parameters introduced by this plugged measurement slice
#' @family likelihood specification functions
#' @family plug-and-play functions
#' @export
#' 


add_meas_BrS_param_Nest_Slice_jags <- function(s,Mobs,cause_list) {
  # mapping template (by `make_template` function):
  patho_BrS_list <- lapply(Mobs$MBS,colnames)
  template_BrS_list <-
    lapply(patho_BrS_list,make_template,cause_list) # key.
  
  # create variable names:
  BrS_nm   <- names(Mobs$MBS)
  # index measurement slices by numbers:
  JBrS_nm  <- paste("JBrS",seq_along(BrS_nm),sep = "_")#
  #PsiBS.cut_nm   <- paste("PsiBS.cut",seq_along(BrS_nm),sep = "_")#
  ThetaBS_nm <- paste("ThetaBS",seq_along(BrS_nm),sep="_")#
  PsiBS_nm <- paste("PsiBS",seq_along(BrS_nm),sep="_")#
  Lambda_nm <- paste("Lambda",seq_along(BrS_nm),sep="_")#
  Lambda0_nm <- paste("Lambda0",seq_along(BrS_nm),sep="_")#
  r0_nm <- paste("r0",seq_along(BrS_nm),sep="_")#
  alphadp0_nm <- paste("alphadp0",seq_along(BrS_nm),sep="_")#
  alphadp0_case_nm <- paste("alphadp0_case",seq_along(BrS_nm),sep="_")# <--- added case subclass weights.
  
  Eta0_nm <- paste("Eta0",seq_along(BrS_nm),sep="_")#
  r1_nm <- paste("r1",seq_along(BrS_nm),sep="_")#
  Eta_nm <- paste("Eta",seq_along(BrS_nm),sep="_")#
  alphaB_nm     <- paste("alphaB",seq_along(BrS_nm),sep = "_")#
  betaB_nm     <- paste("betaB",seq_along(BrS_nm),sep = "_")#
  K_nm <- paste("K",seq_along(BrS_nm),sep="_")#
  
  templateBS_nm  <-
    paste("templateBS",seq_along(BrS_nm),sep = "_")
  indBS_nm  <- paste("indBS",seq_along(BrS_nm),sep = "_")
  
  if (length(patho_BrS_list[[s]]) == 1){stop("==cannot do nested modeling for BrS measurement with 1 dimension!==")} 
  plug <- paste0(
    "
    ####################################
    ### stick-breaking prior specification
    ####################################
    # control subclass mixing weights:
    ",Lambda0_nm[s],"[1]<-",r0_nm[s],"[1]
    ",r0_nm[s],"[",K_nm[s],"]<-1
    for(j in 2:",K_nm[s],") {",Lambda0_nm[s],"[j]<-",r0_nm[s],"[j]*(1-",r0_nm[s],"[j-1])*",Lambda0_nm[s],"[j-1]/",r0_nm[s],"[j-1]}
    for(k in 1:(",K_nm[s],"-1)){
    ",r0_nm[s],"[k]~dbeta(1,",alphadp0_nm[s],")T(0.000001,0.999999)
    }
    
    for (k in 1:(",K_nm[s],"-1)){",Lambda_nm[s],"[k]<-max(0.000001,min(0.999999,",Lambda0_nm[s],"[k]))}
    ",Lambda_nm[s],"[",K_nm[s],"]<-1-sum(",Lambda_nm[s],"[1:(",K_nm[s],"-1)])
    
    # case subclass mixing weights:
    ",Eta0_nm[s],"[1]<-",r1_nm[s],"[1]
    ",r1_nm[s],"[",K_nm[s],"]<-1
    for(j in 2:",K_nm[s],") {",Eta0_nm[s],"[j]<-",r1_nm[s],"[j]*(1-",r1_nm[s],"[j-1])*",Eta0_nm[s],"[j-1]/",r1_nm[s],"[j-1]}
    for(k in 1:(",K_nm[s],"-1)){
    ",r1_nm[s],"[k]~dbeta(1,",alphadp0_case_nm[s],")T(0.000001,0.999999)
    }
    
    for (k in 1:(",K_nm[s],"-1)){",Eta_nm[s],"[k]<-max(0.000001,min(0.999999,",Eta0_nm[s],"[k]))}
    ",Eta_nm[s],"[",K_nm[s],"]<-1-sum(",Eta_nm[s],"[1:(",K_nm[s],"-1)])
    
    ",alphadp0_nm[s],"~dgamma(.25,.25)T(0.001,20)
    ",alphadp0_case_nm[s],"~dgamma(.25,.25)T(0.001,20)
    
    #########################
    ## priors on TPR and FPR:
    #########################
    
    for (j in 1:",JBrS_nm[s],"){
    for (s in 1:",K_nm[s],"){
    ",PsiBS_nm[s],"[j,s]~dbeta(1,1)
    #ThetaBS[j,s]~dbeta(1,1)
    ",ThetaBS_nm[s],"[j,s]~dbeta(",alphaB_nm[s],"[j],",betaB_nm[s],"[j])
    }
    }
    "
  )
  
  parameters <- c(PsiBS_nm[s], ThetaBS_nm[s], Lambda_nm[s],Eta_nm[s],alphaB_nm[s],betaB_nm[s])
  make_list(plug,parameters)
}

##############################################################################
###############  PLUG-AND-PLAY for REGRESSION MODELS         #################
##############################################################################


#' add likelihood component for a BrS measurement slice among cases 
#' 
#' regression model with no nested subclasses
#' 
#' @param s the slice
#' @param Mobs See \code{data_nplcm} described in \code{\link{nplcm}}
#' @param cause_list the list of causes in \code{data_nplcm} described in \code{\link{nplcm}}
#' @param ppd Default is NULL; Set to TRUE for enabling posterior predictive checking.
#' @return a list of two elements: the first is \code{plug}, the .bug code; 
#' the second is \code{parameters} that stores model parameters introduced by this 
#' plugged measurement slice
#' @family likelihood specification functions
#' @family plug-and-play functions
#' @export
add_meas_BrS_case_NoNest_reg_Slice_jags <- function(s,Mobs,cause_list,ppd=NULL) {
  # mapping template (by `make_template` function):
  patho_BrS_list <- lapply(Mobs$MBS,colnames)
  template_BrS_list <-
    lapply(patho_BrS_list,make_template,cause_list) # key.
  
  # create variable names:
  BrS_nm   <- names(Mobs$MBS)
  # index measurement slices by numbers:
  JBrS_nm  <- paste("JBrS",seq_along(BrS_nm),sep = "_") #
  MBS_nm   <- paste("MBS",seq_along(BrS_nm),sep = "_")#
  mu_bs_nm   <- paste("mu_bs",seq_along(BrS_nm),sep = "_")#
  thetaBS_nm <- paste("thetaBS",seq_along(BrS_nm),sep = "_")#
  #psiBS.cut_nm <- paste("psiBS.cut",seq_along(BrS_nm),sep = "_")#
  templateBS_nm  <- paste("templateBS",seq_along(BrS_nm),sep = "_")#
  indBS_nm  <- paste("indBS",seq_along(BrS_nm),sep = "_")#
  Icat_nm   <- "Icat"
  linpred_psiBS_nm <- paste("linpred_psiBS",seq_along(BrS_nm),sep = "_")
  
  if (length(patho_BrS_list[[s]]) > 1) {
    plug <-
      paste0(
        "
        # case BrS measurements; non-nested:
        for (j in 1:",JBrS_nm[s],"){
        ",indBS_nm[s],"[i,j] <- equals(1,",templateBS_nm[s],"[",Icat_nm,"[i],j])
        ",MBS_nm[s],"[i,j]   ~ dbern(",mu_bs_nm[s],"[i,j])
        ",mu_bs_nm[s],"[i,j] <- ",indBS_nm[s],"[i,j]*",thetaBS_nm[s],"[j]+(1-",indBS_nm[s],"[i,j])*ilogit(",linpred_psiBS_nm[s],"[i,j])
        }","\n"
      )
  } else{
    plug <-
      paste0(
        "
        # case BrS measurement; non-nested (with only one column):
        ",indBS_nm[s],"[i] <- equals(1,",templateBS_nm[s],"[",Icat_nm,"[i]])
        ",MBS_nm[s],"[i]   ~ dbern(",mu_bs_nm[s],"[i])
        ",mu_bs_nm[s],"[i] <- ",indBS_nm[s],"[i]*",thetaBS_nm[s],"+(1-",indBS_nm[s],"[i])*ilogit(",linpred_psiBS_nm[s],"[i,1])
        \n"
      )
  }
  
  parameters <- c(Icat_nm, thetaBS_nm[s],linpred_psiBS_nm[s])
  # if posterior predictive distribution is requested:
  if (!is.null(ppd) && ppd){
    MBS_nm.new     <- paste("MBS.new",seq_along(BrS_nm),sep = "_")#
    mu_bs_nm.new   <- paste("mu_bs.new",seq_along(BrS_nm),sep = "_")#
    indBS_nm.new   <- paste("indBS.new",seq_along(BrS_nm),sep = "_")#
    Icat_nm.new    <- "Icat.new"  
    
    if (length(patho_BrS_list[[s]]) > 1) {
      plug <-
        paste0(
          "
          # case BrS measurement; non-nested:
          for (j in 1:",JBrS_nm[s],"){
          ",indBS_nm[s],"[i,j] <- equals(1,",templateBS_nm[s],"[",Icat_nm,"[i],j])
          ",MBS_nm[s],"[i,j]   ~ dbern(",mu_bs_nm[s],"[i,j])
          ",mu_bs_nm[s],"[i,j] <- ",indBS_nm[s],"[i,j]*",thetaBS_nm[s],"[j]+(1-",indBS_nm[s],"[i,j])*ilogit(",linpred_psiBS_nm[s],"[i,j])
          # posterior predictive distribution:
          ",indBS_nm.new[s],"[i,j] <- equals(1,",templateBS_nm[s],"[",Icat_nm.new[s],"[i],j])
          ",MBS_nm.new[s],"[i,j]   ~ dbern(",mu_bs_nm.new[s],"[i,j])
          ",mu_bs_nm.new[s],"[i,j] <- ",indBS_nm.new[s],"[i,j]*",thetaBS_nm[s],"[j]+(1-",indBS_nm.new[s],"[i,j])*ilogit(",linpred_psiBS_nm[s],"[i,j])
          }","\n"
        )
    } else{
      plug <-
        paste0(
          "
          # case BrS measurement; non-nested (with only one column):
          ",indBS_nm[s],"[i] <- equals(1,",templateBS_nm[s],"[",Icat_nm,"[i]])
          ",MBS_nm[s],"[i]   ~ dbern(",mu_bs_nm[s],"[i])
          ",mu_bs_nm[s],"[i] <- ",indBS_nm[s],"[i]*",thetaBS_nm[s],"+(1-",indBS_nm[s],"[i])*ilogit(",linpred_psiBS_nm[s],"[i,1])
          # posterior predictive distribution:
          ",indBS_nm.new[s],"[i] <- equals(1,",templateBS_nm[s],"[",Icat_nm.new[s],"[i]])
          ",MBS_nm.new[s],"[i]   ~ dbern(",mu_bs_nm.new[s],"[i])
          ",mu_bs_nm.new[s],"[i] <- ",indBS_nm.new[s],"[i]*",thetaBS_nm[s],"+(1-",indBS_nm.new[s],"[i])*ilogit(",linpred_psiBS_nm[s],"[i,1])
          
          \n"
        )
    }
    parameters <- c(Icat_nm,Icat_nm.new, thetaBS_nm[s],linpred_psiBS_nm[s])
  }
  
  make_list(plug,parameters)
}

#' add parameters for a BrS measurement slice among cases and controls
#' 
#' regression model with no nested subclasses
#' 
#' @inheritParams add_meas_BrS_case_NoNest_reg_Slice_jags
#' @param FPR_formula False positive regression formula for slice s of BrS data.
#' Check \code{model_options$likelihood$FPR_formula[[s]]}.
#' 
#' @return a list of two elements: the first is \code{plug}, the .bug code; the second is \code{parameters}
#' that stores model parameters introduced by this plugged measurement slice
#' @family likelihood specification functions
#' @family plug-and-play functions 
#' @export
add_meas_BrS_param_NoNest_reg_Slice_jags <- function(s,Mobs,cause_list,FPR_formula) {
  
  constant_FPR <- is_intercept_only(FPR_formula[[s]])
  # mapping template (by `make_template` function):
  patho_BrS_list <- lapply(Mobs$MBS,colnames)
  template_BrS_list <-
    lapply(patho_BrS_list,make_template,cause_list) # key.
  
  has_basis <- ifelse(length(grep("s_date",FPR_formula[[s]]))==0,FALSE,TRUE)
  
  # create variable names:
  BrS_nm   <- names(Mobs$MBS)
  # index measurement slices by numbers:
  JBrS_nm  <- paste("JBrS",seq_along(BrS_nm),sep = "_")#
  thetaBS_nm <- paste("thetaBS",seq_along(BrS_nm),sep = "_")#
  #psiBS.cut_nm <- paste("psiBS.cut",seq_along(BrS_nm),sep = "_")#
  alphaB_nm     <- paste("alphaB",seq_along(BrS_nm),sep = "_")#
  betaB_nm     <- paste("betaB",seq_along(BrS_nm),sep = "_")#
  betaFPR_nm     <- paste("betaFPR",seq_along(BrS_nm),sep = "_")#
  basis_id_nm     <- paste("basis_id",seq_along(BrS_nm),sep = "_")#
  non_basis_id_nm     <- paste("non_basis_id",seq_along(BrS_nm),sep = "_")#
  n_basis_nm     <- paste("n_basis",seq_along(BrS_nm),sep = "_")#
  prec_first_nm     <- paste("prec_first",seq_along(BrS_nm),sep = "_")#
  taubeta_nm     <- paste("taubeta",seq_along(BrS_nm),sep = "_")#
  taubeta0_nm     <- paste("taubeta0",seq_along(BrS_nm),sep = "_")#
  taubeta_inv_nm     <- paste("taubeta_inv",seq_along(BrS_nm),sep = "_")#
  flexible_select_nm     <- paste("flexible_select",seq_along(BrS_nm),sep = "_")#
  ind_flex_select_nm     <- paste("ind_flex_select",seq_along(BrS_nm),sep = "_")#
  p_flexible_nm     <- paste("p_flexible",seq_along(BrS_nm),sep = "_")#
  linpred_psiBS_nm     <- paste("linpred_psiBS",seq_along(BrS_nm),sep = "_")#
  Z_FPR_nm     <- paste("Z_FPR",seq_along(BrS_nm),sep = "_")#
  
  if (length(patho_BrS_list[[s]]) > 1) {
    if (!constant_FPR){
    plug <- paste0(
      "
      # BrS measurement characteristics - non-nested:
      ",linpred_psiBS_nm[s]," <- ",Z_FPR_nm[s],"%*%",betaFPR_nm[s]," # <--- Z_FPR_1: rows for cases and controls, columns for covariates; betaFPR_1: rows for covariates, columns for 1:JBrS_1, i.e., pathogens. 
      
      for (j in 1:",JBrS_nm[s],"){
          ",thetaBS_nm[s],"[j]~ dbeta(",alphaB_nm[s],"[j],",betaB_nm[s],"[j])\n")
    }else{
      plug <- paste0(
        "
        for (j in 1:",JBrS_nm[s],"){
            ",thetaBS_nm[s],"[j]~ dbeta(",alphaB_nm[s],"[j],",betaB_nm[s],"[j])
            # BrS measurement characteristics - non-nested:
            ",linpred_psiBS_nm[s],"[1:(Nd+Nu),j] <- ",Z_FPR_nm[s],"*",betaFPR_nm[s],"[1,j] # <--- Z_FPR_1: rows for cases and controls, columns for covariates; betaFPR_1: rows for covariates, columns for 1:JBrS_1, i.e., pathogens. 
        \n")      
      
      }
    if(has_basis){
      plug <- paste0(plug,
                     "# B-spline basis coefficients:
          ",betaFPR_nm[s],"[",basis_id_nm[s],"[1],j] ~ dnorm(0,",prec_first_nm[s],")
          for (l in 2:",n_basis_nm[s],"){# iterate over the vector of B-spline basis.
            ",betaFPR_nm[s],"[",basis_id_nm[s],"[l],j] ~ dnorm(",betaFPR_nm[s],"[",basis_id_nm[s],"[l-1],j],",taubeta_nm[s],"[j])
          }
          # select flexible semiparametric regression:
          ",taubeta0_nm[s],"[j,1]      ~ dgamma(3,2)               # <-------- flexible fit.
          ",taubeta_inv_nm[s],"[j]     ~ dpar(1.5,0.0025)          # <--------constant fit.
          ",taubeta0_nm[s],"[j,2]      <- pow(",taubeta_inv_nm[s],"[j],-1)
          ",flexible_select_nm[s],"[j] ~ dbern(",p_flexible_nm[s],")
          ",ind_flex_select_nm[s],"[j] <- 2-",flexible_select_nm[s],"[j]
          ",taubeta_nm[s],"[j]         <- ",taubeta0_nm[s],"[j,",ind_flex_select_nm[s],"[j]]\n")
    }
    plug <- paste0(plug,
          "# non-basis coefficients:
          for (l in ",non_basis_id_nm[s],"){
              ",betaFPR_nm[s],"[l,j] ~ dnorm(0,",prec_first_nm[s],")
          }
        }
        ",prec_first_nm[s]," <- 1/4
        # hyperprior of smoothness:
        ",p_flexible_nm[s]," ~ dbeta(1,1) # <--- useless if no basis functions are used.
        "
    )
  } else{
    if (!constant_FPR){
    plug <-
      paste0(
        "
        ",linpred_psiBS_nm[s],"[1:(Nd+Nu),1] <- ",Z_FPR_nm[s],"%*%",betaFPR_nm[s],"[,1] # <--- Z_FPR_1: rows for cases and controls, columns for covariates; betaFPR_1: rows for covariates, columns for 1:JBrS_1, i.e., pathogens. 
        # BrS measurement characteristics - non-nested (only one column):
        ",thetaBS_nm[s],"~ dbeta(",alphaB_nm[s],",",betaB_nm[s],")
        ")
    }else{
      plug <-
        paste0(
          "
          # BrS measurement characteristics - non-nested (only one column):
          ",thetaBS_nm[s],"~ dbeta(",alphaB_nm[s],",",betaB_nm[s],")
          ",linpred_psiBS_nm[s],"[1:(Nd+Nu),1] <- ",Z_FPR_nm[s],"*",betaFPR_nm[s],"[,1] # <--- Z_FPR_1: rows for cases and controls, columns for covariates; betaFPR_1: rows for covariates, columns for 1:JBrS_1, i.e., pathogens. 
          ")      
      }
    if(has_basis){
      plug <-
        paste0(plug,
        "# B-spline basis coefficients:
        ",betaFPR_nm[s],"[",basis_id_nm[s],"[1],1] ~ dnorm(0,",prec_first_nm[s],")
        for (l in 2:",n_basis_nm[s],"){# iterate over the vector of B-spline basis.
            ",betaFPR_nm[s],"[",basis_id_nm[s],"[l],1] ~ dnorm(",betaFPR_nm[s],"[",basis_id_nm[s],"[l-1],1],",taubeta_nm[s],")
        }
        # select flexible semiparametric regression:
        ",taubeta0_nm[s],"[1]     ~ dgamma(3,2)              # <-------- flexible fit.
        ",taubeta_inv_nm[s],"     ~ dpar(1.5,0.0025)          # <--------constant fit.
        ",taubeta0_nm[s],"[2]      <- pow(",taubeta_inv_nm[s],",-1)
        ",flexible_select_nm[s]," ~ dbern(",p_flexible_nm[s],")
        ",ind_flex_select_nm[s],"   <- 2-",flexible_select_nm[s],"
        ",taubeta_nm[s],"        <- ",taubeta0_nm[s],"[",ind_flex_select_nm[s],"]\n")
    }
    plug <-
      paste0(plug,
        "# non-basis coefficients:
        for (l in ",non_basis_id_nm[s],"){
            ",betaFPR_nm[s],"[l,1] ~ dnorm(0,",prec_first_nm[s],")
        }
        ",prec_first_nm[s]," <- 1/4
        # hyperprior of smoothness:
        ",p_flexible_nm[s]," ~ dbeta(1,1)  # <--- useless if no basis functions are used.
        "
      )
  }
  parameters <- c(thetaBS_nm[s],betaFPR_nm[s],alphaB_nm[s],betaB_nm[s],taubeta_nm[s],p_flexible_nm[s])
  make_list(plug,parameters)
}


#' add a likelihood component for a BrS measurement slice among controls
#' 
#' regression model without nested subclasses
#' 
#' @inheritParams add_meas_BrS_case_NoNest_reg_Slice_jags
#' 
#' @return a list of two elements: the first is \code{plug}, the .bug code; the second is \code{parameters}
#' that stores model parameters introduced by this plugged measurement slice
#' @family likelihood specification functions
#' @family plug-and-play functions
#' @export
add_meas_BrS_ctrl_NoNest_reg_Slice_jags <- function(s, Mobs,cause_list,ppd=NULL) {
  # mapping template (by `make_template` function):
  patho_BrS_list <- lapply(Mobs$MBS,colnames)
  template_BrS_list <-
  lapply(patho_BrS_list,make_template,cause_list) # key.
  
  # create variable names:
  BrS_nm   <- names(Mobs$MBS)
  # index measurement slices by numbers:
  JBrS_nm  <- paste("JBrS",seq_along(BrS_nm),sep = "_")#
  MBS_nm   <- paste("MBS",seq_along(BrS_nm),sep = "_")#
  mu_bs_nm   <- paste("mu_bs",seq_along(BrS_nm),sep = "_")#
  psiBS_nm   <-  paste("psiBS",seq_along(BrS_nm),sep = "_")#
  linpred_psiBS_nm   <-  paste("linpred_psiBS",seq_along(BrS_nm),sep = "_")#
  
  if (length(patho_BrS_list[[s]]) > 1) {
    plug <- paste0(
      "       
      ## control BrS measurements; no subclass:
      for (j in 1:",JBrS_nm[s],"){
      ",MBS_nm[s],"[i,j] ~ dbern(",mu_bs_nm[s],"[i,j])
      logit(",mu_bs_nm[s],"[i,j]) <- ",linpred_psiBS_nm[s],"[i,j] 
      }
      "
    )
  } else{
    plug <- paste0(
      "
      ## control BrS measurements; no subclass (only one column):
      ",MBS_nm[s],"[i] ~ dbern(",mu_bs_nm[s],"[i])
      logit(",mu_bs_nm[s],"[i]) <- ",linpred_psiBS_nm[s],"[i,1] 
      \n"
    )
  }
  
  if (!is.null(ppd) && ppd){
    
    MBS_nm.new   <- paste("MBS.new",seq_along(BrS_nm),sep = "_")#
    mu_bs_nm.new   <- paste("mu_bs.new",seq_along(BrS_nm),sep = "_")#
    
    if (length(patho_BrS_list[[s]]) > 1) {
      plug <- paste0(
        "       
        ## control BrS measurements; no subclass:
        for (j in 1:",JBrS_nm[s],"){
        ",MBS_nm[s],"[i,j] ~ dbern(",mu_bs_nm[s],"[i,j])
        logit(",mu_bs_nm[s],"[i,j]) <- ",linpred_psiBS_nm[s],"[i,j]
        ## posterior predictive distribution
        ",MBS_nm.new[s],"[i,j] ~ dbern(",mu_bs_nm[s],"[i,j])
        }
        "
      )
    } else{
      plug <- paste0(
        "
        ## control BrS measurements; no subclass (only one column):
        ",MBS_nm[s],"[i] ~ dbern(",mu_bs_nm[s],"[i])
        logit(",mu_bs_nm[s],"[i]) <- ",linpred_psiBS_nm[s],"[i,1]
        ## posterior predictive distribution:
        ",MBS_nm.new[s],"[i] ~ dbern(",mu_bs_nm[s],"[i])
        \n"
      )
    }
  }
  
  parameters <- c(psiBS_nm[s])
  make_list(plug,parameters)
}
