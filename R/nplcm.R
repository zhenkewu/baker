#' Fit nested partially latent class models (high-level)
#'
#' Uses WinBUGS/JAGS in Windows system or JAGS in OSX system for Bayesian inference
#' (see README file for an instruction to install WinBUGS or JAGS).
#' For stratification/regression functionalities, true positive rates are constant
#' across strata or covariate values.\cr
#' 
#' @param data_nplcm Must put cases at the top in `Mobs$MBS` and control data at the bottom
#' (this is related to the coding scheme in JAGS).
#' \itemize{
#' \item  \code{Mobs} A list of measurements. The elements of the list
#' should include \code{MBS}, \code{MSS}, and \code{MGS}. If any of the component
#' is not available, please specify it as, e.g., \code{MGS=NULL}
#' (effectively deleting \code{MGS} from \code{Mobs}).
#' \itemize{
#' \item \code{MBS} a list of data frame of bronze-standard (BrS) measurements.
#' Rows are subjects, columns are pathogens. We use list here to accommodate
#' the possibility of multiple sets of BrS data.
#' They have imperfect sensitivity/specificity (e.g. nasopharyngeal PCR). 
#' \item \code{MSS} a list of data frame of silver-standard (SS) measurements. 
#' Rows are subjects, columns are pathogens measured in specimen (e.g. blood culture).
#' These measurements have perfect specificity but imperfect sensitivity.
#' \item \code{MGS} a list of data frame of gold-standard (GS) measurements. 
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
#' 
#' \item \code{likelihood}{
#'     \itemize{
#'          \item{cause_list} The vector of latent status;
#'          \item{k_subclass} The number of nested subclasses. 1 for conditional independence,
#' >1 for conditional dependence; It is a vector of length equal to the number of slices
#' of BrS measurements;
#'          \item{Eti_formula} formula for etiology regressions. You can use 
#' \code{\link{s_date_Eti}} to specify the design matrix for R format enrollment date;
#' it will produce natural cubic spline basis. Specify \code{~ 1} if no
#' regression is intended.
#'          \item{FPR_formula}formula for false positive rates (FPR) regressions; see
#' \code{\link{formula}}. You can use \code{\link{s_date_FPR}} to specify part
#' of the design matrix for R format enrollment date; it will produce penalized-spline
#' basis (based on B-splines). Specify \code{~ 1} if no
#' regression is intended. NB: If \code{effect="fixed"}, \code{\link{dm_Rdate_FPR}}
#' will just specify a design matrix with appropriately standardized dates. 
#'     }
#' }
#' \item \code{use_measurements}{
#'      a vector of characters strings; can be any singleton or combinations of "BrS", "SS", "GS".
#' }
#' 
#' \item \code{prior}{
#'     \itemize{
#'          \item{Eti_prior}Description of etiology prior 
#' (e.g., \code{overall_uniform} - all hyperparameters are 1; or \code{0_1} - all hyperparameters
#' are 0.1);
#'          \item{TPR_prior}Description of priors for the measurements 
#' (e.g., informative vs non-informative). 
#' Its length should be the same with \code{M_use};
#'     }
#' }
#' }
#' @param mcmc_options A list of Markov chain Monte Carlo (MCMC) options.
#'
#' \itemize{
#' \item \code{debugstatus} Logical - whether to pause WinBUGS after it finishes
#' model fitting;
#' \item \code{n.chains} Number of MCMC chains;
#' \item \code{n.burnin} Number of burn-in samples;
#' \item \code{n.thin} To keep every other \code{n.thin} samples after burn-in period;
#' \item \code{individual.pred} \code{TRUE} to perform individual prediction (\code{Icat}
#' variables in the \code{.bug} file);
#' \code{FALSE} otherwise;
#' \item \code{ppd} \code{TRUE} to simulate new data (\code{XXX.new}
#' variables in the \code{.bug} file) from the posterior predictive 
#' distribution (ppd); \code{FALSE} otherwise;
#' \item \code{get.pEti} \code{TRUE} for getting posterior samples of individual etiologic fractions; 
#' \code{FALSE} otherwise. For non-regression, or regression models with all discrete predictors, 
#' this is defaulted to be \code{TRUE}; no need to specify this entry. It is only relevant for regression models
#' with non-discrete covariates. Because individuals have distinct etiology pies at their specific covariate values, 
#' it's easier to just store the posterior samples of the regression coefficients and reconstruct the pies afterwards,
#' rather than storing them through JAGS. 
#' \item \code{result.folder} Path to folder storing the results;
#' \item \code{bugsmodel.dir} Path to \code{.bug} model files;
#' \item \code{winbugs.dir} Path to where WinBUGS 1.4 is installed.
#' }
#' @return A JAGS result, fitted by function \code{\link[R2jags]{jags2}} from \code{R2jags}. 
#' Current implemented models follow the hierarchy below:
#' \itemize{
#' \item no regression:  \link{write_model_NoReg}
#'        
#' \item regression: 
#'  Given disease class (control or a subtype of cases infected with the same sets of pathogens)
#'       \itemize{
#'          \item local independence model for BrS measures: \link{write_model_Reg_NoNest}
#'          \item local dependence model for BrS measures: \link{write_model_Reg_Nest}
#'        }
#' }
#' @export
nplcm <- function(data_nplcm,model_options,mcmc_options){
  Mobs <- data_nplcm$Mobs
  Y    <- data_nplcm$Y
  if (rle(Y)>2 | Y[1]!=1) {stop("==[baker] 'data_nplcm$Y' must have cases on top of controls. Use 'baker::subset_data_nplcm_by_index()' to shuffle the rows. Then retry.==\n")}
  X    <- data_nplcm$X
  
  parsed_model <- assign_model(model_options,data_nplcm)
  
  do_reg_vec <- unlist(parsed_model$regression[grep("^do_reg_",names(parsed_model$regression))])
  do_discrete_vec <- unlist(parsed_model$regression[grep("^is_discrete_predictor",names(parsed_model$regression))])
  do_discrete_and_reg <- do_reg_vec&do_discrete_vec
  do_reg <- any(do_reg_vec)
  do_discrete <- sum(do_discrete_and_reg)==sum(do_reg_vec) 
  # only call nplcm_fit_Reg_discrete_predictor_NoNest
  # when all regressions are using discrete predictors!

  do_nested <- parsed_model$nested
  
  if(!do_reg){
      res <- nplcm_fit_NoReg(data_nplcm,model_options,mcmc_options)
  } 
  if (do_reg & !any(do_nested)){
    if (do_discrete){
      # when every regression is a regression upon discrete variables;
      # when it is not a regression, the fitting function treats it as a single stratum
      # when specifying the model in the .bug file (in the assign_model function, ~1 
      # is considered not a regression, i.e., covariate indendependent - this is the
      # custom in this package.)
      res <- nplcm_fit_Reg_discrete_predictor_NoNest(data_nplcm,model_options,mcmc_options)
      } else{
      # a mix of discrete and continuous regressions:
      res <- nplcm_fit_Reg_NoNest(data_nplcm,model_options,mcmc_options)
      }
  }
  if (do_reg & any(do_nested)){
    res <- nplcm_fit_Reg_Nest(data_nplcm,model_options,mcmc_options)
    cat("==[baker] Regression model with nested subclasses finally done.==\n")
  }
  res
}
