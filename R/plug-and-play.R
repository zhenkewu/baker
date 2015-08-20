#
# 1. Bronze-standard data: conditional independence model (no nested structure).
#

#' add a likelihood component for a BrS measurement slice among cases (conditional independence)
#' 
#' @param s the slice
#' @param Mobs See \code{data_nplcm} described in \code{\link{nplcm}}
#' @param cause_list the list of causes in \code{data_nplcm} described in \code{\link{nplcm}}
#' 
#' @return a list of two elements: the first is \code{plug}, the .bug code; 
#' the second is \code{parameters} that stores model parameters introduced by this 
#' plugged measurement slice
#' 
#' @export
add_meas_BrS_case_NoNest_Slice <- function(s,Mobs,cause_list) {
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
  
  if (length(patho_BrS_list[[s]]) > 1) {
    plug <-
      paste0(
        "
          # case BrS measurement; non-nested:
          for (j in 1:",JBrS_nm[s],"){
          ",indBS_nm[s],"[i,j] <- equals(1,",templateBS_nm[s],"[Icat[i],j])
          ",MBS_nm[s],"[i,j] ~ dbern(",mu_bs_nm[s],"[i,j])
          ",mu_bs_nm[s],"[i,j]<-", indBS_nm[s],"[i,j]*",thetaBS_nm[s],"[j]+(1-", indBS_nm[s],"[i,j])*",psiBS.cut_nm[s],"[j]
          }","\n"
      )
  } else{
    plug <-
      paste0(
        "
           
          # case BrS measurement; non-nested (with only one column):
          ",indBS_nm[s],"[i] <- equals(1,",templateBS_nm[s],"[Icat[i]])
          ",MBS_nm[s],"[i] ~ dbern(",mu_bs_nm[s],"[i])
          ",mu_bs_nm[s],"[i]<-", indBS_nm[s],"[i]*",thetaBS_nm[s],"+(1-", indBS_nm[s],"[i])*",psiBS.cut_nm[s],"\n"
      )
  }
  
  parameters <- c("Icat",thetaBS_nm[s],psiBS.cut_nm[s])
  make_list(plug,parameters)
}




#' add a likelihood component for a BrS measurement slice among controls (conditional independence)
#' 
#' 
#' @inheritParams add_meas_BrS_case_NoNest_Slice
#' 
#' @return a list of two elements: the first is \code{plug}, the .bug code; the second is \code{parameters}
#' that stores model parameters introduced by this plugged measurement slice
#' 
#' @export
add_meas_BrS_ctrl_NoNest_Slice <- function(s, Mobs,cause_list) {
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
#' 
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
#' 
#' @export
#' 
add_meas_BrS_case_Nest_Slice <- function(s,Mobs,cause_list){
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
  
  if (length(patho_BrS_list[[s]]) == 1){stop("==cannot do nested modeling for BrS measurement with 1 dimension!==")} 
  
  plug <-
    paste0(
      "
          ## case BrS measurements; with subclasses:
          for (j in 1:",JBrS_nm[s],"){
            ",indBS_nm[s],"[i,j] <- equals(1,",templateBS_nm[s],"[Icat[i],j])
            ",MBS_nm[s],"[i,j]~dbern(",mu_bs.bound_nm[s],"[i,j])
            ",mu_bs.bound_nm[s],"[i,j]<-max(0.000001,min(0.999999,",mu_bs_nm[s],"[i,j]))
            ",mu_bs_nm[s],"[i,j]<-",PR_BS_nm[s],"[i,j,",Z_nm[s],"[i]]
            
            for (s in 1:",K_nm[s],"){
              ",PR_BS_nm[s],"[i,j,s]<-",PsiBS.cut_nm[s],"[j,s]*(1-",indBS_nm[s],"[i,j])+",ThetaBS_nm[s],"[j,s]*",indBS_nm[s],"[i,j]
            }
          }
      "
    )
  
  parameters <- c("Icat",ThetaBS_nm[s])
  make_list(plug,parameters)
}

#' add likelihood for a BrS measurement slice among controls (conditional independence)
#' 
#' 
#' @inheritParams add_meas_BrS_case_NoNest_Slice
#' 
#' @return a list of two elements: the first is \code{plug}, the .bug code; the second is \code{parameters}
#' that stores model parameters introduced by this plugged measurement slice
#' 
#' @export
#' 
add_meas_BrS_ctrl_Nest_Slice <- function(s, Mobs,cause_list) {
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
#' 
#' @export
#' 


add_meas_BrS_param_Nest_Slice <- function(s,Mobs,cause_list) {
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
    for(k in 1:",K_nm[s],"-1){
      ",r0_nm[s],"[k]~dbeta(1,",alphadp0_nm[s],")I(0.000001,0.999999)
    }
    
    for (k in 1:",K_nm[s],"-1){",Lambda_nm[s],"[k]<-max(0.000001,min(0.999999,",Lambda0_nm[s],"[k]))}
    ",Lambda_nm[s],"[",K_nm[s],"]<-1-sum(",Lambda_nm[s],"[1:(",K_nm[s],"-1)])
    
    # case subclass mixing weights:
    ",Eta0_nm[s],"[1]<-",r1_nm[s],"[1]
    ",r1_nm[s],"[",K_nm[s],"]<-1
    for(j in 2:",K_nm[s],") {",Eta0_nm[s],"[j]<-",r1_nm[s],"[j]*(1-",r1_nm[s],"[j-1])*",Eta0_nm[s],"[j-1]/",r1_nm[s],"[j-1]}
    for(k in 1:",K_nm[s],"-1){
          ",r1_nm[s],"[k]~dbeta(1,",alphadp0_nm[s],")I(0.000001,0.999999)
    }
    
    for (k in 1:",K_nm[s],"-1){",Eta_nm[s],"[k]<-max(0.000001,min(0.999999,",Eta0_nm[s],"[k]))}
    ",Eta_nm[s],"[",K_nm[s],"]<-1-sum(",Eta_nm[s],"[1:(",K_nm[s],"-1)])
    
    ",alphadp0_nm[s],"~dgamma(.25,.25)I(0.001,20)
    
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
#' 
#' @export

add_meas_BrS_subclass_Nest_Slice <- function(s,Mobs,cause_list){
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
  
  parameters <- c("Icat",Eta_nm[s],Lambda_nm[s])
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
#' 
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
#' 
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
                    ","\n"
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
                    \n"
          )
      }
      
    }
  }
  
  plug <- paste0(res,collapse = "")
  parameters <- c(thetaSS_nm[s],psiSS_nm[s], alphaS_nm[s],betaS_nm[s])
  make_list(plug,parameters)
}
