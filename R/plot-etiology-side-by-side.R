if(getRversion() >= "2.15.1") utils::globalVariables(c("cause","probability","DIR"))

##' compare the posterior distribution of population etiologies side-by-side
#'
#' @details \code{plot_etiology_side_by_side} produces what we call "music-sheet plot" 
#' to compare estimates of population distribution of latent health status 
#' (a.k.a population etiology, or population etiology pie). It reads in two or more
#' folders where Bayesian inference results are combined. We implemented a check to make sure
#' that pathogens from one analysis are subset of another. NB: current implementation
#' does not check whether the results are from the same stratum, e.g., same site. 
#' We recommend putting information in the the folder names under comparisons.
#'
#'@param DIR_list The list of directory paths, each storing a model output.
#'
#'@param DIR_pathogen_displayorder_lookup The directory path to the .csv file
#'that stores the display order of pathogens in the combined music sheet plot.
#'
#'@param dodge_val default is 0.5; for width and position of boxplots.
#'@param right_panel default is \code{TRUE}, for bacterial, viral, combo, other groupings. 
#'Set to \code{FALSE} if not wanted. 
#'@param reg_ind A vector of \code{TRUE} or \code{FALSE}, indicating each directory in the DIR_list
#'contains a regression result or not
#' 
#'@import ggplot2
#'@import reshape2
#'
#'@return A figure that compares posterior etiology distribution stored in
#'two or more folders
#'
#'@examples
#'\dontrun{
#' PATH <- "C:\\package_test"
#'
#' dir_list <- list(file.path(PATH,"SouthAfrica"),
#'                  file.path(PATH,"Kenya"))
#' library(baker)
#' plot_etiology_side_by_side(dir_list,
#'   file.path(PATH,"pathogen_displayorder_lookup.csv"))
#'}
#' @family visualization functions
#' @family comparison functions
#'@export
#'
plot_etiology_side_by_side <- function(DIR_list,
                                       DIR_pathogen_displayorder_lookup,
                                       dodge_val = 0.5,
                                       right_panel = TRUE,
                                       reg_ind = NULL){
  old_par <- graphics::par(graphics::par("mfrow", "mar"))
  on.exit(graphics::par(old_par))
  ## read in pathogen display order lookup table:
  pathogen_displayorder_lookup <- utils::read.csv(DIR_pathogen_displayorder_lookup)
  f                            <- pathogen_displayorder_lookup$Cause
  display_order                <- as.character(levels(f))[f] # <--- this is where the display order is specified.
  
  if (any(duplicated(display_order))){stop("== There are duplicated names in the names to specify the displaying order! Please retain only one of the duplicates. ==")}
  # read from folders:
  out_list <- vector("list",length(DIR_list))
  base_nm  <- lapply(DIR_list,basename)
  names(out_list) <- base_nm
  
  # get names of causes (a union of all the supplied result folders):
  union_names <- list()
  pEti_samp_list  <- list()
  for (i in seq_along(DIR_list)){
    out_list[[i]]        <- nplcm_read_folder(DIR_list[[i]])
    union_names          <- c(union_names,out_list[[i]]$model_options$likelihood$cause_list)
    if (is.null(reg_ind)){
      pEti_samp_list[[i]]  <- get_pEti_samp(out_list[[i]]$res_nplcm,out_list[[i]]$model_options)
    }else{
      if (!reg_ind[i]){
        pEti_samp_list[[i]]  <- get_pEti_samp(out_list[[i]]$res_nplcm,out_list[[i]]$model_options)
      } else{
        pEti_samp_list[[i]]  <- get_pEti_samp_reg(DIR_list[[i]],out_list[[i]]$model_options)
      }
    }
  }
  
  union_names <- unlist(union_names)
  
  NDIR        <- length(DIR_list)
  # get unique names of causes:
  uniq_names    <- unique_cause(union_names)
  # order the uniq_names according to the order that we want to display the pathogens
  # in the comparison:
  if (!all(uniq_names%in%display_order)){
    missed <- setdiff(uniq_names,display_order)  
    stop(paste0("==Please add ",paste(missed,collpase=", "), " in the 'DIR_pathogen_displayorder_lookup'!=="))
  }
  order_to_disp <- my_reorder(display_order,uniq_names)
  # the names organized in an order which we want to read from each folder's results:
  read_names    <- uniq_names[order_to_disp]
  # the vector of combo sizes for the unioned causes that are ordered and ready for 
  # reading and displaying results:
  fitted_num    <- unlist(lapply(strsplit(read_names,"\\+"),length))
  # get combined pEti samples listed by the combo sizes:
  res_all_combo <- vector("list",length(unique(fitted_num)))
  
  res   <- vector("list",NDIR)  
  for (d in seq_along(DIR_list)){
    NSAMP       <- nrow(out_list[[d]]$res_nplcm) # no. of retained MCMC iterations.
    tmp           <- matrix(0,nrow = NSAMP,ncol = length(read_names))
    colnames(tmp) <- read_names
    res[[d]]      <- as.data.frame(tmp)
    curr_DIR_cause <- pEti_samp_list[[d]]$latent_nm  
    read_order     <- match_cause(curr_DIR_cause,read_names)
    for (j in seq_along(read_order)){
      if (!is.na(read_order[j])) {res[[d]][,j] <- pEti_samp_list[[d]]$pEti_mat[,read_order[j]] }
    }
    res[[d]]$DIR  <- d
    res[[d]]$ITER <- 1:NSAMP
  }
  
  res_cbind <- do.call(rbind,res)
  
  if (right_panel){
    ind_virus <- get_cause_by_taxo_group(read_names,"virus",pathogen_displayorder_lookup)
    ind_bact  <- get_cause_by_taxo_group(read_names,"bacterium",pathogen_displayorder_lookup)
    ind_combo <- get_cause_by_taxo_group(read_names, "combo",pathogen_displayorder_lookup)
    ind_other <- which(read_names=="other")
    
    res_cbind$Virus <- rowSums(res_cbind[,ind_virus,drop=FALSE])
    res_cbind$Bacteria  <- rowSums(res_cbind[,ind_bact,drop=FALSE])
    if (length(ind_combo)>0){
      res_cbind$Combo <- rowSums(res_cbind[,ind_combo,drop=FALSE])
    }
    if (length(ind_other)>0)
      res_cbind$Other <- rowSums(res_cbind[,ind_other,drop=FALSE])
  }
  
  # first build some functions to summarize posterior distribution 
  # (following ggplot2 syntax):
  f <- function(x) {
    r <- stats::quantile(x, probs = c(0.05, 0.25, 0.5, 0.75, 0.95))
    names(r) <- c("ymin", "lower", "middle", "upper", "ymax")
    r
  }
  mean_with_nm <- function(x){
    r <- rep(mean(x),2)
    names(r)<-c("y","ymax")
    r
  }
  mean_with_nm_txt <- function(x){
    r <- c(ifelse(max(x)-stats::quantile(x,.97)>0.02,stats::quantile(x,.97)+0.02,max(x)),
           round(mean(x),3),round(mean(x),3)*100)
    names(r)<-c("y","ymax","label")
    r
  }
  
  curr_res         <- res_cbind
  data_for_boxplot <- reshape2::melt(curr_res,
                                     id.vars = c("DIR","ITER"),
                                     variable.name="cause",
                                     value.name = "probability")
  ## ggplot2:
  ymax <- max(data_for_boxplot$probability)
  aes_now <- function(...) {
    structure(list(...),  class = "uneval")
  }
  
  gg<-ggplot(data = data_for_boxplot, 
             aes(x = factor(cause), y = probability, fill = factor(DIR))) +
    labs(list(x = "cause", y = "probability"))+theme_bw()+
    stat_summary(fun.data = f, geom="boxplot",aes_now(width=dodge_val),
                 position = position_dodge(dodge_val))+
    stat_summary(fun.data = mean_with_nm,geom="point",aes(size=1.5),
                 position = position_dodge(dodge_val))+scale_size(guide = 'none')+
    stat_summary(fun.data = mean_with_nm_txt,geom="text",
                 aes(angle=90),position = position_dodge(width = dodge_val))+
    scale_fill_discrete("Model\n",labels = c(base_nm))+
    guides(fill=guide_legend(nrow=NDIR,byrow=TRUE))+
    theme(legend.text = element_text(colour="blue",size = 16, face = "bold"),
          legend.title = element_text(size=16,face="bold"),legend.position = "top",
          axis.title   = element_text(size=16,face="bold"),
          axis.text.x = element_text(angle=40, vjust=.8, hjust=1.01,size=16,face="bold"))+
    scale_y_continuous(limits = c(0,ymax))
  if (right_panel){
    gg <- gg+
      geom_vline(xintercept=length(read_names)+.5,linetype = "longdash")
  }
  gg
  
}

#' a function to match each model's result to the desired order:
#' 
#' @param disp_order the desired order
#' @param union_nm the union of all cause names
#' @param nm_list the list of names for each model (i.e. folder)
#' 
#' @return a list of two elements:
#' \itemize{
#' \item \code{model_res} a list of length same as \code{DIR_list} in 
#' \code{\link{plot_etiology_side_by_side}}. Each element is the order
#' for that particular name in the list; the length is the same as union_nm.
#' \item \code{union_res} the names of the causes after union and re-ordering
#' }
#' 
#' @examples  
#'   disp_order <- c("B","E","D","C","F","A")
#'   union_nm   <- c("A","B","C","D","E")
#'   nm_list <-list(c("C","E"),
#'                  c("D","B"),
#'                  c("C","A","E"))
#'   lookup_side_by_side(disp_order,union_nm,nm_list)
#' @export
lookup_side_by_side <- function(disp_order, union_nm, nm_list){
  ord_union <- rep(NA,length(union_nm))
  incre <- 0
  for (j in seq_along(disp_order)){
    if (disp_order[j]%in% union_nm){
      incre <- incre + 1
      ord_union[incre] <- which(union_nm==disp_order[j])
    }
  }
  pick_order <- union_nm[ord_union]
  res <- list()
  for (i in seq_along(nm_list)){
    res[[i]] <- rep(NA,length(pick_order))
    for (j in seq_along(pick_order)){
      if (pick_order[j]%in%nm_list[[i]]){
        res[[i]][j] <- which(nm_list[[i]]==pick_order[j])
      }
    }
  }
  list(model_res=res, union_res = pick_order)
}


#' get etiology samples by names (no regression)
#' 
#' @inheritParams order_post_eti
#' 
#' @return A list:
#' \itemize{
#'    \code{pEti_mat}: a matrix of posterior samples (iteration by cause); overall etiology
#'    \code{latent_nm}: a vector of character strings representing the names of the causes
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


#' get etiology samples by names (no regression)
#' 
#' @param DIR_NPLCM File path to the folder containing posterior samples
#' @param model_options See \code{\link{nplcm}}
#' @param stratum_bool a vector of TRUE/FALSE with TRUE indicating the rows of subjects to include
#' @param pEti_subject TRUE for getting individual specific samples; Default to \code{FALSE}
#' @param reg_param TRUE for getting regression parameters; Default to \code{FALSE}
#' @param truth a list of truths computed from true parameters in simulations; elements: 
#'  Eti, FPR, PR_case,TPR; All default to \code{NULL} in real data analyses.
#'  Currently only works for one slice of bronze-standard measurements (in a non-nested model).
#'  \itemize{
#'      \item Eti matrix of # of rows = # of subjects, # columns: \code{length(cause_list)} for Eti
#'      \item FPR matrix of # of rows = # of subjects, # columns: \code{ncol(data_nplcm$Mobs$MBS$MBS1)}
#'      \item PR_case matrix of # of rows = # of subjects, # columns: \code{ncol(data_nplcm$Mobs$MBS$MBS1)}
#'      \item TPR a vector of length identical to \code{PR_case}
#'  }
#' @param return_metric TRUE for showing overall mean etiology, quantiles, s.d., and if \code{truth$Eti} is supplied, 
#'  coverage, bias, truth and integrated mean squared errors (IMSE).
#' @param RES_NPLCM pre-read res_nplcm; default to NULL to save time.
#' 
#' @return A list:
#' \itemize{
#'    \code{pEti_mat}: a matrix of posterior samples (iteration by cause); overall etiology
#'    \code{latent_nm}: a vector of character strings representing the names of the causes
#'    \code{pEti_subject}: an array of individual specific etiology; cause by subject by iteration
#'    \code{betaEti}: an array of samples of the etiology regression parameters; coefficient by cause by iteration
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

#' Match latent causes that might have the same combo but
#' different specifications
#' 
#'  @details In our cause_list, "A+B" represents the same cause
#'   as "B+A". It is used for plotting side-by-side posterior sample comparisons
#' 
#' @param pattern a vector of latent cause names, e.g., from a particular fit
#' @param vec a vector of latent cause names, e.g., usually a union of cause names
#' from several model fits. Usually, it is also the display order that one wants to 
#' show.
#' 
#' @return A vector of length \code{length(vec)}; \code{NA} means no pattern matches
#' vec; 1 at position 10 means the first element of \code{pattern} matches the 
#' 10th element of \code{vec}.
#' 
#' 
#' @examples 
#' 
#' pattern <- c("X+Y","A+Z","C")
#' vec     <- c(LETTERS[1:26],"Y+Z","Y+X","Z+A")
#' match_cause(pattern,vec)
#' 
#' @export
match_cause <- function(pattern, vec){
  has_plus_sign <- grepl("\\+",pattern)
  vec_split <- strsplit(vec,"\\+")
  res <- rep(NA,length(vec))
  for (p in seq_along(pattern)){
    tmp <- pattern[p]
    if (has_plus_sign[p]) {
      tmp <- strsplit(pattern[p],"\\+")[[1]]
    }
    find_match <- lapply(vec_split,setequal,y=tmp)
    res[which(find_match==TRUE)] <- p
  }
  res
}


#' get the names of causes by the number of pathogens associated with that cause.
#' For now it is specifically for PERCH.
#' 
#' @param cause_list see \code{model_options}
#' @param num the number of pathogens for a cause
#' 
#' @examples 
#' 
#' x <- c("A","B","C","A+B","C+D","AA+BD","A+B+C","HELLO+BYE")
#' num <- 2
#' get_cause_by_num(x,num)
#' 
#' @return a subvector of cause_list
#' 
#' @export
get_cause_by_num <- function(cause_list,num){
  cause_split <- strsplit(cause_list,"\\+")
  ind <- which(lapply(cause_split,length)==num)
  cause_list[ind]
}


#' get unique causes, regardless of the actual order in combo
#' 
#' @param cause_vec a vector of characters with potential combo repetitions
#' written in scrambled orders separated by "+"
#' 
#' @examples 
#' x <- c("A","B","A","CC+DD","DD+CC","E+F+G","B")
#' unique_cause(x)
#' 
#' @return a vector of characters with unique meanings for latent causes
#' 
#' @export

unique_cause <- function(cause_vec){
  cause_split <- strsplit(cause_vec,"\\+")
  num <- lapply(cause_split,length)
  res <- list()
  res[[1]] <- cause_split[[1]] 
  
  if (length(cause_vec)==1){return(cause_vec)}
  
  j <- 1
  for (i in 2:length(cause_vec)){
    got_new <- TRUE
    for (iter in 1:j){
      if (setequal(cause_split[[i]],res[[iter]])){got_new <- FALSE}
    }
    if (got_new){j <- j+1;res[[j]]<-cause_split[[i]]}
  }
  unlist(lapply(res,paste,collapse="+"))
}

#' pick categories from a vector of causes
#' 
#' @param cause_list see \code{model_options}
#' @param taxo_group a character string, can be "virus", "bacterium", "combo".
#' @param lookup_table a data frame with two columns ("Category" and "Cause").
#' 
#' @return a vector of indices indicating if each element belongs to 
#' \code{taxo_group}.
#' 
#' @export

get_cause_by_taxo_group <- function(cause_list,taxo_group,lookup_table){
  if (is.null(lookup_table$Category)){stop("==No Category information!==")}
  if (is.null(lookup_table$Cause)){stop("==No Cause information!==")}
  ind_subset <- list()
  for (e in seq_along(cause_list)){
    ind <- which(lookup_table$Cause==cause_list[e])
    ind_subset[[e]] <- ind
  }
  ind_subset <- unlist(ind_subset)
  if (taxo_group == "combo"){
    res <- grep("\\+",lookup_table$Category[ind_subset]) 
    return(res)
  }
  res <- grep(taxo_group,lookup_table$Category[ind_subset])
}

# cause_list <- colnames(res_cbind)[seq_along(read_names)]
# taxo_group <- "virus"
# lookup_table <- pathogen_displayorder_lookup
# get_cause_by_taxo_group(colnames(res_cbind)[seq_along(read_names)],
#                         "combo",pathogen_displayorder_lookup)






