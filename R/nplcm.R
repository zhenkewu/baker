#' Fit nested partially latent class models (high-level)
#'
#' Uses WinBUGS in Windows system for Bayesian inference
#' (see README file for an instruction to install WinBUGS on windows 7 or 8.).
#' For stratification/regression functionalities, true positive rates are constant
#' across stratum or covariate values.\cr
#' 
#' @details Developer Note (DN):\cr
#' \itemize{
#' \item need to add checking empty cell functionality when performing
#'          stratification
#' \item need to use formula \code{cbind(Y1, Y2)~X1+X2} kind of symbolic
#'  model specification
#' to enable model estimation.
#' \item need to add ppd regression functionality
#'}
#'
#' @param data_nplcm
#' \itemize{
#' \item  \code{Mobs} A list of measurements. The elements of the list
#' should include \code{MBS}, \code{MSS}, and \code{MGS}. If any of the component
#' is not available, please specify it as, e.g., \code{MGS=NA} or \code{MGS=NULL}
#' (effectively deleting \code{MGS} from \code{Mobs}).
#' \itemize{
#' \item \code{MBS} a data frame of bronze-standard (BrS) measurements.
#' Rows are subjects, columns are pathogens. 
#' They have imperfect sensitivity/specificity (e.g. nasalpharyngeal PCR). 
#' \item \code{MSS} a data frame of silver-standard (SS) measurements. 
#' Rows are subjects, columns are pathogens measured in specimen (e.g. blood culture).
#' These measurements have perfect specificity but imperfect sensitivity.
#' \item \code{MGS} a data frame of gold-standard (GS) measurements. 
#' Rows are subject, columns are pathogen measurements. 
#' These measurements have perfect sensitivity and specificity.
#' }
#'
#' \item \code{Y} Vector of disease status: 1 for case, 0 for control.
#' \item \code{X} Covariate matrix for regression modeling. It contains raw covariate
#' data, not design matrix for regression models.
#' }
#' @param model_options A list of model options.
#'
#' \itemize{
#' \item \code{M_use} List of measurements to be used in the model;
#' \item \code{k_subclass}The number of nested subclasses. 1 for conditional independence,
#' >1 for conditional dependence;
#' \item \code{TPR_prior} Description of priors for the measurements 
#' (e.g., informative vs non-informative). 
#' Its length should be the same with \code{M_use};
#' \item \code{Eti_prior} Description of etiology prior 
#' (e.g., \code{overall_uniform} - all hyperparameters are 1; or \code{0_1} - all hyperparameters
#' are 0.1);
#' \item \code{pathogen_BrS_list} The vector of pathogen names with BrS measure;
#' \item \code{cause_list} The vector of causes that are either singleton or 
#' combinations of items in \code{pathogen_BrS_list}; 'NoA' can also be included
#' at the end, meaning 'None-of-Above';
#' \item \code{X_reg_FPR} formula for false positive rates (FPR) regressions; see
#' \code{\link{formula}}. You can use \code{\link{dm_Rdate_FPR}} to specify part
#' of the design matrix for R format enrollment date; it will produce thin-plate
#' regression spline basis for every date (if \code{effect="random"} and \code{num_knots_FPR} is
#' specified to a positive integer, e.g., 10.). If \code{effect="fixed"}, \code{\link{dm_Rdate_FPR}}
#' will just specify a design matrix with appropirately standardized dates. Specify \code{~0} if no
#' regression is intended.
#' \item \code{X_reg_Eti} formula for etiology regressions. You can use 
#' \code{\link{dm_Rdate_Eti}} to specify the design matrix for R format enrollment date;
#' it will produce natural cubic splines for every date. Specify \code{~0} if no
#' regression is intended.
#' \item \code{pathogen_BrS_cat} The two-column data frame that has category of pathogens: virus (V), bacteria (B)
#' and fungi (F);
#' \item \code{pathogen_SSonly_list} The vector of pathogens with only
#' SS measure;
#' \item \code{pathogen_SSonly_cat} The category of pathogens with only SS measure.
#'}
#'
#' @param mcmc_options A list of Markov chain Monte Carlo (MCMC) options.
#'
#' \itemize{
#' \item \code{debugstatus} Logical - whether to pause WinBUGS after it finishes
#' model fitting;
#' \item \code{n.chains} Number of MCMC chains;
#' \item \code{n.burnin} Number of burn-in samples;
#' \item \code{n.thin} To keep every other \code{n.thin} samples after burn-in period;
#' \item \code{individual.pred} whether to perform individual prediction;
#' \item \code{ppd} whether to perform posterior predictive (ppd) checking;
#' \item \code{result.folder} Path to folder storing the results;
#' \item \code{bugsmodel.dir} Path to WinBUGS model files;
#' \item \code{winbugs.dir} Path to where WinBUGS 1.4 is installed.
#' }
#' @return A WinBUGS result, fitted by function \code{bugs()} from
#' the \code{\link{R2WinBUGS}} package. Current implemented models follow the hierarchy below:
#' \itemize{
#' \item no regression: 
#'    \itemize{
#'    \item both BrS and SS measurements:
#' 
#'    \itemize{
#'    \item no pathogens with only SS measure:
#'        \itemize{
#'          \item independence model
#'          \item dependence model
#'        }
#'    \item there are some pathogens with SS only measure:
#'      \itemize{
#'        \item independence model
#'        \item dependence model
#'      }
#'    
#'    }
#'    
#'    \item only BrS measurements:
#'       \itemize{
#'        \item independence model
#'        \item dependence model
#'       }
#'    }
#' }
#'
#' @export

nplcm <- function(data_nplcm,model_options,mcmc_options){
  Mobs <- data_nplcm$Mobs
  Y    <- data_nplcm$Y
  X    <- data_nplcm$X
  parsing <- assign_model(data_nplcm,model_options)
  
  if (!any(unlist(parsing$reg))){
    # if no stratification or regression:
    if (parsing$measurement$quality=="BrS+SS"){
          if (!parsing$measurement$SSonly){
              if (!parsing$measurement$nest){
                # model 1, DONE
                    res <- nplcm_fit_NoReg_BrSandSS_NoNest(data_nplcm,model_options,mcmc_options)
              }else{
                # model 2, DONE
                    res <- nplcm_fit_NoReg_BrSandSS_Nest(data_nplcm,model_options,mcmc_options)
              }
          } else{
              if (!parsing$measurement$nest){
                # model 3, DONE
                res <- nplcm_fit_NoReg_BrSandSS_NoNest_SSonly(data_nplcm,model_options,mcmc_options)
              }else{
                # model 4, DONE
                res <- nplcm_fit_NoReg_BrSandSS_Nest_SSonly(data_nplcm,model_options,mcmc_options)
              }
          }
    }else if (parsing$measurement$quality=="BrS"){
          if (!parsing$measurement$SSonly){
              if (!parsing$measurement$nest){
                # model 5, DONE
                   res <- nplcm_fit_NoReg_BrS_NoNest(data_nplcm,model_options,mcmc_options)
              }else{
                # model 6, DONE
                   res <- nplcm_fit_NoReg_BrS_Nest(data_nplcm,model_options,mcmc_options)
              }
          }
    }
  } else{
      # if stratification or regression:
      if (parsing$measurement$quality=="BrS+SS"){
        if (!parsing$measurement$SSonly){
          if (!parsing$measurement$nest){
            stop("== Done but need to clean code. Please contact maintainer. Thanks.")
          }else{
            # model 2, DONE
            stop("== Done but need to clean code. Please contact maintainer. Thanks.")
            
          }
        } else{
          if (!parsing$measurement$nest){
            stop("== Done but need to clean code. Please contact maintainer. Thanks.")
          }else{
            # model 4, DONE
            stop("== Done but need to clean code. Please contact maintainer. Thanks.")
          }
        }
      }else if (parsing$measurement$quality=="BrS"){
        if (!parsing$measurement$SSonly){
          if (!parsing$measurement$nest){
            #if (is_discrete(X,c(model_options$X_reg_FPR,model_options$X_reg_Eti))){
            #    res <- nplcm_fit_Strat_BrS_NoNest(data_nplcm,model_options,mcmc_options)
            #} else{
            #  stop("==Not done.==")
            #}
            
            res <- nplcm_fit_Reg_BrS_NoNest(data_nplcm,model_options,mcmc_options)
          }else{
            stop("== Done but need to clean code. Please contact maintainer. Thanks.")
          }
        }
      }
  }
  res
}
