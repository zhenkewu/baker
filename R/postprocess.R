#' get etiology samples by names (no regression)
#' 
#' @inheritParams order_post_eti
#' 
#' @return A list:
#' \itemize{
#'    `pEti_mat`: a matrix of posterior samples (iteration by cause); overall etiology
#'    `latent_nm`: a vector of character strings representing the names of the causes
#' }
#' 
#' @export

get_pEti_samp <- function(res_nplcm,model_options){
  cause_list <- model_options$likelihood$cause_list
  # total no. of causes:
  Jcause     <- length(cause_list)
  # extract and process some data and posterior samples:
  SubVarName <- rep(NA,Jcause)
  for (j in 1:Jcause){
    SubVarName[j] = paste("pEti","[",j,"]",sep="")
  }
  # get etiology fraction MCMC samples:
  pEti_mat   <- as.data.frame(res_nplcm[,SubVarName,drop=FALSE])
  latent_nm  <- model_options$likelihood$cause_list
  make_list(pEti_mat,latent_nm)
}


#' get etiology samples by names under regression models
#' 
#' @param DIR_NPLCM File path to the folder containing posterior samples
#' @param model_options See [nplcm()]
#' @param stratum_bool a vector of TRUE/FALSE with TRUE indicating the rows of subjects to include
#' @param pEti_subject TRUE for getting individual specific samples; Default to `FALSE`
#' @param reg_param TRUE for getting regression parameters; Default to `FALSE`
#' @param truth a list of truths computed from true parameters in simulations; elements: 
#'  Eti, FPR, PR_case,TPR; All default to `NULL` in real data analyses.
#'  Currently only works for one slice of bronze-standard measurements (in a non-nested model).
#'  \itemize{
#'      \item Eti matrix of # of rows = # of subjects, # columns: `length(cause_list)` for Eti
#'      \item FPR matrix of # of rows = # of subjects, # columns: `ncol(data_nplcm$Mobs$MBS$MBS1)`
#'      \item PR_case matrix of # of rows = # of subjects, # columns: `ncol(data_nplcm$Mobs$MBS$MBS1)`
#'      \item TPR a vector of length identical to `PR_case`
#'  }
#' @param return_metric TRUE for showing overall mean etiology, quantiles, s.d., and if `truth$Eti` is supplied, 
#'  coverage, bias, truth and integrated mean squared errors (IMSE).
#' @param RES_NPLCM pre-read res_nplcm; default to NULL to save time.
#' 
#' @return A list:
#' \itemize{
#'    `pEti_mat`: a matrix of posterior samples (iteration by cause); overall etiology
#'    `latent_nm`: a vector of character strings representing the names of the causes
#'    `pEti_subject`: an array of individual specific etiology; cause by subject by iteration
#'    `betaEti`: an array of samples of the etiology regression parameters; coefficient by cause by iteration
#' }
#' 
#' @export
get_pEti_samp_reg <- function(DIR_NPLCM,model_options,
                              stratum_bool=NULL,
                              pEti_subject = FALSE,
                              reg_param = FALSE,
                              truth = NULL,
                              return_metric=FALSE,
                              RES_NPLCM = NULL
){
  cause_list <- model_options$likelihood$cause_list
  # total no. of causes:
  Jcause     <- length(cause_list)
  # extract and process some data and posterior samples:
  SubVarName <- rep(NA,Jcause)
  for (j in 1:Jcause){
    SubVarName[j] = paste("pEti","[",j,"]",sep="")
  }
  
  #
  # Read data from DIR_NPLCM:
  #
  data_nplcm <- dget(file.path(DIR_NPLCM,"data_nplcm.txt"))  
  model_options <- dget(file.path(DIR_NPLCM,"model_options.txt"))
  mcmc_options <- dget(file.path(DIR_NPLCM,"mcmc_options.txt"))
  parsed_model <- assign_model(model_options,data_nplcm)
  
  
  Z_Eti       <- stats::model.matrix(model_options$likelihood$Eti_formula,
                                     data.frame(data_nplcm$X,Y=data_nplcm$Y)[data_nplcm$Y==1,,drop=FALSE])
  
  if (is.null(stratum_bool)){stratum_bool <- 1:length(data_nplcm$Y)}
  
  is_nested    <- parsed_model$nested
  new_env <- new.env()
  source(file.path(DIR_NPLCM,"jagsdata.txt"),local=new_env)
  bugs.dat <- as.list(new_env)
  rm(new_env)
  if (!is.null(RES_NPLCM)){res_nplcm <- RES_NPLCM
  } else {res_nplcm <- coda::read.coda(file.path(DIR_NPLCM,"CODAchain1.txt"),
                                       file.path(DIR_NPLCM,"CODAindex.txt"),
                                       quiet=TRUE)}
  print_res <- function(x) plot(res_nplcm[,grep(x,colnames(res_nplcm))])
  get_res   <- function(x) res_nplcm[,grep(x,colnames(res_nplcm))]
  
  # structure the posterior samples:
  n_samp_kept   <- nrow(res_nplcm)
  ncol_dm_Eti   <- ncol(bugs.dat$Z_Eti)
  Jcause        <- bugs.dat$Jcause
  Nd            <- bugs.dat$Nd
  Nu            <- bugs.dat$Nu
  
  # #####################################################################
  # # add x-axis for dates:
  # X <- data_nplcm$X
  # Y <- data_nplcm$Y
  # # some date transformations:
  # X$date_plot  <- as.Date(X$ENRLDATE)
  # X$date_month_centered <- as.Date(cut(X$date_plot,breaks="2 months"))+30
  # X$date_month <- as.Date(cut(X$date_plot,breaks="2 months"))
  # 
  # dd <-  as.Date(X$ENRLDATE)
  # min_d <- min(dd)
  # min_d_std <- unique(X$std_date[which(X$ENRLDATE==min_d)])
  # min_plot_d <- min_d+days_in_month(month(min_d))-day(min_d)+1
  # 
  # max_d <- max(dd)
  # max_d_std <- unique(X$std_date[which(X$ENRLDATE==max_d)])
  # max_plot_d <- max_d-day(max_d)+1
  # plot_d <- seq.Date(min_plot_d,max_plot_d,by = "quarter")
  # 
  # unit_x <- (max_d_std-min_d_std)/as.numeric(max_d-min_d)
  # plot_d_std <- as.numeric(plot_d - min_d)*unit_x+min_d_std
  # 
  # pred_d <- seq.Date(min_plot_d,max_plot_d,by = "day")
  # pred_d_std <- as.numeric(pred_d - min_d)*unit_x+min_d_std
  # #####################################################################
  
  # pred_dataframe <- data.frame(ENRLDATE=as.POSIXct.Date(pred_d,tz="UTC"),
  #                              t(replicate(length(pred_d),unlist(unique(X[discrete_names])[1,]))))
  # if (nrow(unique(X[discrete_names]))>1){
  #   for (l in 2:nrow(unique(X[discrete_names]))){
  #     pred_dataframe <- rbind(pred_dataframe,
  #                             data.frame(ENRLDATE=as.POSIXct.Date(pred_d,tz="UTC"),
  #                                        t(replicate(length(pred_d),unlist(unique(X[discrete_names])[l,])))))
  #   }
  # }
  # 
  # 
  # pred_dataframe$std_date <- dm_Rdate_FPR(c(pred_dataframe$ENRLDATE,data_nplcm$X$ENRLDATE),
  #                                         c(rep(1,nrow(pred_dataframe)),data_nplcm$Y),
  #                                         effect = "fixed")[-(1:nrow(data_nplcm$X))]
  # pred_dataframe_ok <- cbind(pred_dataframe,Y=rep(1,nrow(pred_dataframe)))
  # 
  # Z_Eti_pred       <- stats::model.matrix(model_options$likelihood$Eti_formula,
  #                                         pred_dataframe_ok)
  # 
  
  if (!is_nested){
    betaEti_samp <- array(t(get_res("^betaEti")),c(ncol_dm_Eti,Jcause,n_samp_kept))
    linpred      <- function(beta,design_matrix){design_matrix%*%beta}
    out_Eti_linpred     <- array(apply(betaEti_samp,3,linpred,design_matrix=bugs.dat$Z_Eti),
                                 c(Nd,Jcause,n_samp_kept))
  } else{
    betaEti_samp <- array(t(get_res("^betaEti")),c(ncol_dm_Eti,Jcause,n_samp_kept),dimnames = list(colnames(Z_Eti),
                                                                                                   model_options$likelihood$cause_list,
                                                                                                   1:n_samp_kept)) #useful in effect estimation.
    linpred      <- function(beta,design_matrix){design_matrix%*%beta}
    
    # out_caseFPR_linpred     <- array(apply(case_betaFPR_samp,3,linpred,design_matrix=bugs.dat[[paste0("Z_FPR_",slice)]]),
    #                              c(Nd+Nu,K_curr,n_samp_kept))
    out_Eti_linpred     <- array(apply(betaEti_samp,3,linpred,design_matrix=bugs.dat$Z_Eti),
                                 c(Nd,Jcause,n_samp_kept)) # can potentially just add pEti to the monitoring.
    #pEti_samp           <- apply(out_Eti_linpred,c(1,3),softmax) # Jcause by Nd by niter.
    pEti_samp           <- abind::abind(aperm(apply(out_Eti_linpred,c(1,3),softmax),c(2,1,3)),
                                        array(0,c(Nu,Jcause,n_samp_kept)),along=1)
  }
  
  # etiology:
  subset_Eti <- data_nplcm$Y==1 & stratum_bool # <--- specifies who to look at.
  plotid_Eti <- which(subset_Eti)[order(data_nplcm$X$std_date[subset_Eti])]
  curr_date_Eti  <- data_nplcm$X$std_date[plotid_Eti]
  
  if (!is_nested){
    Eti_prob_scale <- apply(out_Eti_linpred[plotid_Eti,,],c(1,3),softmax)
  }else{Eti_prob_scale <- aperm(pEti_samp,c(2,1,3))[,plotid_Eti,]}
  Eti_mean <- apply(Eti_prob_scale,c(1,2),mean)
  Eti_q    <- apply(Eti_prob_scale,c(1,2),quantile,c(0.025,0.975))
  Eti_overall <- apply(Eti_prob_scale,c(1,3),mean)
  Eti_overall_mean <- rowMeans(Eti_overall)
  Eti_overall_sd   <- apply(Eti_overall,1,sd)
  Eti_overall_q    <- apply(Eti_overall,1,quantile,c(0.025,0.975))
  
  # get etiology fraction MCMC samples:
  pEti_mat   <- t(Eti_overall)
  colnames(pEti_mat) <- SubVarName
  latent_nm  <- model_options$likelihood$cause_list
  res <- make_list(pEti_mat,latent_nm)
  
  if (pEti_subject){
    res <- append(res,list(pEti_subject=Eti_prob_scale)) # pathogen by subjects by iteration.
  }
  if (reg_param){
    res <- append(res,list(betaEti = betaEti_samp)) # coefficient by cause by iteration.
  }
  
  if (return_metric){
    if (!is.null(truth$Eti)){
      Eti_overall_truth  <- colMeans(truth$Eti[plotid_Eti,])
      # Eti_IMSE <- rep(0,length(cause_list))
      # for (t in 1:n_samp_kept){
      #   Eti_IMSE <- Eti_IMSE*(t-1) + apply(Eti_prob_scale[,,t]-t(truth$Eti[plotid_Eti,]),1,function(v) sum(v^2)/length(v))
      #   Eti_IMSE <- Eti_IMSE/t
      # }
      # compute integrated squared error:
      Eti_ISE <- apply(Eti_mean-t(truth$Eti[plotid_Eti,]),1,function(v) sum(v^2)/length(v))
      Eti_overall_cover  <- sapply(seq_along(Eti_overall_mean),
                                   function(s) (Eti_overall_truth[s]<= Eti_overall_q[2,s]) && 
                                     (Eti_overall_truth[s]>= Eti_overall_q[1,s]))
      Eti_overall_bias  <- Eti_overall_mean -  Eti_overall_truth
      res <- append(res,make_list(Eti_overall_mean,Eti_overall_q,Eti_overall_sd,
                                  Eti_overall_cover, Eti_overall_bias,Eti_overall_truth,Eti_ISE))
    } else{
      res <- append(res,make_list(Eti_overall_mean,Eti_overall_q,Eti_overall_sd))
    }
  }
  res
}