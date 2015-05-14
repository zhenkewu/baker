#' Set true positive rate (TPR) prior ranges
#'
#' Current prior assignment let bacteria NPPCR to have uniform range, viral NPPCR
#' to have .5-1.0 range. The PCP (a fungus) NPPCR TPR is also set to be .5-1.0; PCP
#' has no blood culture measurements. Also, not all the bacteria have blood culture
#' measurments. One question is whether to use informative NPPCR range or non-informative
#' NPPCR range (0-100%).
#'
#' DN: 1.make the assignment of prior dependent on the BCX data availability
#' or species (bacteria)?
#'
#' @param model_options See \code{\link{nplcm}} function.
#' @param data_nplcm See \code{\link{assign_model}} function.
#' @return Parameters for the TPR priors, separately for BrS and SS
#'
#' @export

TPR_prior_set <- function(model_options,data_nplcm){
      
      Mobs <- data_nplcm$Mobs
      Y    <- data_nplcm$Y
      X    <- data_nplcm$X
  
      pathogen_BrS_list <- model_options$pathogen_BrS_list
      JBrS         <- length(pathogen_BrS_list)#note that here JBrS=JBrS_BrS in nplcm_fit.

      parsing <- assign_model(data_nplcm,model_options)

      # only BrS data is used:
      if (parsing$measurement$quality=="BrS"){
        temp_cat_ind <- sapply(model_options$pathogen_BrS_list,
                               function(path) {which(model_options$pathogen_BrS_cat$X==path)})
        temp_cat     <- model_options$pathogen_BrS_cat[temp_cat_ind,]

        if (model_options$TPR_prior=="noninformative"){
           alphaB = rep(1,JBrS)
           betaB  = rep(1,JBrS)
        }else if (model_options$TPR_prior=="informative"){
              #virus have BrS measurement sensitivity in the range of 50-100%.
              # currently use beta (6,2). bacterial use 0-100%. currently we use beta(1,1).
              alphaB <- rep(NA,JBrS)
              betaB  <- rep(NA,JBrS)

              for (t in 1:nrow(temp_cat)){
                if (temp_cat$pathogen_type[t]=="B"){
                  alphaB[t] <- 1
                  betaB[t]  <- 1
                }
                if (temp_cat$pathogen_type[t]=="V" | temp_cat$pathogen_type[t]=="F"){
                  alphaB[t] <- 6
                  betaB[t]  <- 2
                }
              }
        }
        res <- list (alphaB = alphaB, betaB = betaB,
                     used_cat = temp_cat)
      }

      # both BrS and SS data are used:
      if (parsing$measurement$quality=="BrS+SS"){
              ## get JSS:----------
              MSS.case <- Mobs$MSS[Y==1,1:JBrS]
              MSS.case <- as.matrix(MSS.case)

              SS_index <- which(colMeans(is.na(MSS.case))<0.9)#.9 is arbitrary; any number <1 will work.
              JSS      <- length(SS_index)
              ##------------------------------

              temp_cat_ind <- sapply(model_options$pathogen_BrS_list,
                                     function(path) {which(model_options$pathogen_BrS_cat$X==path)})
              temp_cat     <- model_options$pathogen_BrS_cat[temp_cat_ind,]

              if (model_options$TPR_prior[1]=="noninformative"){
                   alphaB <- rep(1,JBrS)
                   betaB  <- rep(1,JBrS)
                }else if (model_options$TPR_prior[1]=="informative"){
                  #virus have BrS measurement sensitivity in the range of 50-100%.
                  # currently use beta (6,2). bacterial use 0-100%. currently we use beta(1,1).
                  alphaB <- rep(NA,JBrS)
                  betaB  <- rep(NA,JBrS)

                  for (t in 1:nrow(temp_cat)){
                    if (temp_cat$pathogen_type[t]=="B"){
                      alphaB[t] <- 1
                      betaB[t]  <- 1
                    }
                    if (temp_cat$pathogen_type[t]=="V" | temp_cat$pathogen_type[t]=="F"){
                      alphaB[t] <- 6
                      betaB[t]  <- 2
                    }
                  }
              }

              #need to define JSS here, which is the number of BrS+SS available
              # bacteria. It is not the same as the JSS used in nplcm_three_panel_plot.
              if (model_options$TPR_prior[2]=="noninformative"){
                  alphaS <- rep(1,JSS)
                  betaS  <- rep(1,JSS)
                  if (!is.null(model_options$SSonly) && model_options$SSonly==TRUE){
                          JSS_only    <- length(model_options$pathogen_SSonly_list)
                          alphaS.only <- rep(1,JSS_only)
                          betaS.only  <- rep(1,JSS_only)
                  }
                } else if (model_options$TPR_prior[2]=="informative"){
                  #For bacteria, informative sensitivities lie in the range of
                  # (.05,.15):
                  alphaS <- rep(NA,JSS)
                  betaS  <- rep(NA,JSS)
                  temp_cat_ind <- sapply(model_options$pathogen_BrS_list,
                                         function(path) {which(model_options$pathogen_BrS_cat$X==path)})
                  temp_cat     <- model_options$pathogen_BrS_cat[temp_cat_ind,]
                  for (t in 1:JSS){
                      temp_param <- beta_parms_from_quantiles(c(.05,.15),p=c(0.025,.975),plot=FALSE)
                      alphaS[t] <- temp_param$a
                      betaS[t] <- temp_param$b
                  }
                  if (!is.null(model_options$pathogen_SSonly_list)){
                              JSS_only    <- length(model_options$pathogen_SSonly_list)
                              alphaS.only <- rep(NA,JSS_only)
                              betaS.only  <- rep(NA,JSS_only)
                              for (t in 1:JSS_only){
                                    temp_param_SSonly <- beta_parms_from_quantiles(c(.05,.15),p=c(0.025,.975),plot=FALSE)
                                    alphaS.only[t] <- temp_param_SSonly$a
                                    betaS.only[t] <- temp_param_SSonly$b
                              }
                  }
               }
              if (!is.null(model_options$pathogen_SSonly_list)){
                      res <- list(alphaB = alphaB, betaB = betaB,
                              alphaS = alphaS, betaS = betaS,
                              alphaS.only = alphaS.only,
                              betaS.only  = betaS.only,
                              used_cat = temp_cat)
              }else{
                       res <- list(alphaB = alphaB, betaB = betaB,
                              alphaS = alphaS, betaS = betaS,
                              used_cat = temp_cat)
              }
       }

       res
}

