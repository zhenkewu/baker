if(getRversion() >= "2.15.1") utils::globalVariables(c("pattern","frequency","DIR","..y.."))

#' Posterior predictive checking for the nested partially class models - 
#' frequent patterns in the BrS data. (for multiple folders)
#' 
#' At each MCMC iteration, we generate a new data set based on the model and 
#' parameter values at that iteration. The sample size of the new data set equals
#' that of the actual data set, i.e. the same number of cases and controls.
#' 
#' 
#' 
#' @param DIR_list The list of directory paths, each storing a model output.
#' @param slice_vec Default are 1s, for the first slice of BrS data.
#' @param n_pat Number of the most common BrS measurement pattern among cases and controls.
#' Default is 10.
#' @param dodge_val Default is 0.8; For width of boxplots.
#' 
#' @return A figure of posterior predicted frequencies compared with the observed 
#' frequencies of the most common patterns for the BrS data. 
#' 
#' @import ggplot2
#' @importFrom gridExtra grid.arrange
#' 
#' @examples 
#' \dontrun{
#' DIR_list <- list("C:\\2015_09_04_05SAF_k=1","C:\\2015_09_04_05SAF_k=2")
#' plot_check_common_pattern(DIR_list)
#' plot_check_common_pattern(DIR_list[[1]])
#' }
#' 
#' @export

plot_check_common_pattern <- function(DIR_list,
                                      slice_vec = rep(1,length(DIR_list)),
                                      n_pat     = 10,
                                      dodge_val = 0.8){
  # read in data:
  # names of the measurements for the selected slice:
  name_vec <- vector("list",length(DIR_list))
  # the list of results read from the specified list of directories:
  out      <- vector("list",length(DIR_list))
  is_jags  <- rep(NA,length(DIR_list))
  for (d in seq_along(DIR_list)){
    curr_dir   <- DIR_list[[d]]
    curr_slice <- slice_vec[d]
    # read NPLCM outputs:
    out[[d]]           <- nplcm_read_folder(curr_dir)
    name_vec[[d]]      <- out[[d]]$clean_options$BrS_objects[[curr_slice]]$patho # <-- it means we need to store clean_options in the result folder.
    is_jags[d]         <- is_jags_folder(curr_dir)
  }
  
  if (!(length(unique(name_vec))==1)){
    stop("==The results under comparison have different BrS measurement names! 
         Please use `slice_vec` to match the names.==")
    }
  
  get_top_pattern <- function(curr_out,case_status,slice,n_pat,curr_is_jags){
    # getting data:
    curr_bugs.dat <- curr_out$bugs.dat
    curr_Nd       <- curr_out$Nd
    curr_Nu       <- curr_out$Nu
    curr_slice    <- slice
    # get observed data:
    curr_observed <- curr_out$Mobs$MBS[[curr_slice]]
    # length of a pattern (e.g., 10001 means length is 5)
    len_pat       <- ncol(curr_out$Mobs$MBS[[curr_slice]])
    curr_Y        <- curr_out$Y # case control status.
    curr_res_nplcm <- curr_out$res_nplcm # get posterior samples.
    curr_predicted <- curr_res_nplcm[,grep(paste0("^MBS.new_",curr_slice,"\\["),colnames(curr_res_nplcm)),drop=FALSE]
    NSAMP <- nrow(curr_predicted)
    # organize into an array for easy subsetting (case and control's measurements):
    if (!curr_is_jags){
      curr_predicted_array <- array(curr_predicted,dim=c(NSAMP,len_pat,curr_Nd+curr_Nu)) # <-- first dimension for iterations; second dimension for pathogen measurements; third dimension for individual (cases first and controls).
    } else{
      curr_predicted_array <- aperm(array(curr_predicted,dim=c(NSAMP,curr_Nd+curr_Nu,len_pat)),c(1,3,2))
    }
    # subsetting into case or control based on input `case_status`:
    observed  <- curr_observed[curr_Y==case_status,,drop=FALSE]
    predicted <- curr_predicted_array[,,curr_Y==case_status,drop=FALSE]
    
    # convert numeric vector into a character string: e.g., c(1,0,0,1,1) into "10011"; for faster pattern matching:
    collapse_byrow <- function(mat){
      NA2dot(apply(mat,1,paste,collapse = "" ))
    }
    
    observed_pat  <- collapse_byrow(observed)
    predicted_pat <- as.matrix(apply(predicted,3,collapse_byrow))
    
    # counting patterns:
    pat               <- sort(table(observed_pat),decreasing=TRUE) # observed patten.
    n_pat_used        <- min(n_pat,length(pat)) # actually used pattern number.
    ind_missing       <- grep("\\.",names(pat)) # pick out patterns with missing measurements.
    n_missing         <- sum(pat[ind_missing])  # the total number of individuals with missing measurements.
    pat_high_frac     <- pat[1:n_pat_used]/(length(observed_pat)-n_missing) # divide by no. of individuals with complete measurements.
    pat_high_frac_no_missing     <- pat_high_frac
    if (length(ind_missing)){
      pat_high_frac_no_missing     <- pat_high_frac[-ind_missing] # delete patterns with missingness.
    }
    pat_high_name_no_missing     <- names(pat_high_frac_no_missing) # get names to display in the plot.
    
    n_pat_used_no_missing <- length(pat_high_name_no_missing) # the length of patterns without missingness.
    
    ppd_pat_ct  <- matrix(NA,nrow=NSAMP,ncol=n_pat_used_no_missing)
    exist_other <- (length(pat)-length(ind_missing))>0
    if (exist_other){
      ppd_pat_ct <- matrix(NA,nrow=NSAMP,ncol=n_pat_used_no_missing+1) # the extra column for "other" patterns.
    }
    ppd_pat_ct <- as.data.frame(ppd_pat_ct)
    for (iter in 1:NSAMP){
      ppd_pat_table <- table(predicted_pat[iter,])
      curr_ct <- (sapply(pat_high_name_no_missing,function(x) {
        indtmp <- names(ppd_pat_table)==x
        if (sum(indtmp)==0){
          res <- 0
        } else{
          res <- ppd_pat_table[indtmp]
        }     
        res
      }
      ))
      ppd_pat_ct[iter,1:n_pat_used_no_missing] <- curr_ct/ncol(predicted_pat)
      if (exist_other){
        ppd_pat_ct[iter,n_pat_used_no_missing+1] <- 1-sum(curr_ct)/ncol(predicted_pat)
      }
    }
    colnames(ppd_pat_ct) <- c(1:length(pat_high_name_no_missing))
    pattern_names        <- c(pat_high_name_no_missing)
    obs_pat              <- pat_high_frac_no_missing
    if (exist_other){
       colnames(ppd_pat_ct) <- 1:(length(pat_high_name_no_missing)+1)
       pattern_names  <- c(pat_high_name_no_missing,"other")
       obs_pat        <- c(pat_high_frac_no_missing,1-sum(pat_high_frac_no_missing))
    }
    names(obs_pat) <- pattern_names
    
    make_list(ppd_pat_ct,obs_pat,pattern_names,exist_other)
  }
  
  plot_ppd <- function(DIR_list,case_or_control="case"){
    case_res_list <- vector("list",length=length(DIR_list)) # posterior predictive pattern frequencies.
    case_pat_list <- vector("list",length=length(DIR_list)) # pattern names.
    res_list      <- vector("list",length=length(DIR_list)) # long format results with extra information.
    obs_case_pat_list <- vector("list",length=length(DIR_list)) # observed pattern frequencies.
    
    if (case_or_control=="case"){select <- 1}
    if (case_or_control=="control"){select <- 0}
    
    for (d in seq_along(DIR_list)){
      out_case_pat       <- get_top_pattern(out[[d]],select,slice_vec[d],n_pat,is_jags[d])
      case_res_list[[d]] <- out_case_pat$ppd_pat_ct
      case_res_list[[d]]$DIR <- d 
      case_res_list[[d]]$ITER <- 1:nrow(case_res_list[[d]])
      case_res_list[[d]]$CASE <- rep(select,nrow(case_res_list[[d]]))
      res_list[[d]] <- reshape2::melt(case_res_list[[d]],
                                      id.vars = c("CASE","DIR","ITER"),
                                      variable.name="pattern",
                                      value.name = "frequency")
      case_pat_list[[d]]     <- out_case_pat$pattern_names
      obs_case_pat_list[[d]] <- out_case_pat$obs_pat
    }
    
    if(!(length(unique(case_pat_list))==1)){stop("==Different data sets are used under comparison! Please use the results fitted from the same data set.==")}
    if(!(length(unique(obs_case_pat_list))==1)){stop("==Different data sets are used under comparison! Please use the results fitted from the same data set.==")}
    
    res <- do.call(rbind,res_list) # combine results across directories.
    
    base_nm <- lapply(DIR_list,basename)
    NDIR    <- length(DIR_list)
    # first build some functions to summarize posterior distribution 
    # (following ggplot2 syntax):
    f <- function(x) {
      r <- quantile(x, probs = c(0.05, 0.25, 0.5, 0.75, 0.95))
      names(r) <- c("ymin", "lower", "middle", "upper", "ymax")
      r
    }
    mean_with_nm <- function(x){
      r <- rep(mean(x),2)
      names(r)<-c("y","ymax")
      r
    }
    mean_with_nm_txt <- function(x){
      r <- c(ifelse(max(x)-quantile(x,.97)>0.02,quantile(x,.97)+0.02,max(x)),
             round(mean(x),3),round(mean(x),3)*100)
      names(r)<-c("y","ymax","label")
      r
    }
    
    ## ggplot2:
    ymax <- max(res$frequency)
    aes_now <- function(...) {
      structure(list(...),  class = "uneval")
    }
    
#     case_status_labeller <- function(variable,value){
#       c("Case","Control")[2-value]
#     }
    
    # plot for cases:
    hline.data <- as.data.frame(list(frequency = obs_case_pat_list[[1]],
                                       pattern   = 1:(length(obs_case_pat_list[[1]])),
                                       DIR       = rep(1,length(obs_case_pat_list[[1]]))))
    gg1<-ggplot(data = res, 
                aes(x = factor(pattern), y = frequency, fill = factor(DIR))) +
      #facet_wrap(~ CASE, ncol = 2)+ 
      #facet_grid(~CASE,labeller=case_status_labeller)+
      labs(list(x = "pattern", y = "frequency"))+theme_bw()+
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
            axis.text.x = element_text(angle=40, vjust=.8, hjust=1.01,size=16,face="bold"),
            strip.text.x = element_text(size = 16, colour = "red",face="bold"))+
      scale_x_discrete(labels=case_pat_list[[1]])+
      scale_y_continuous(limits = c(0,ymax))+
      annotate("text", label = case_or_control, 
               x = length(case_pat_list[[1]])/2, 
               y = ymax*0.67, size = 8, colour = "red")+
      geom_errorbar(stat = "hline", 
                    width=0.8,aes(yintercept = frequency,ymax=..y..,ymin=..y..),
                    color="blue",size=0.9,data=hline.data)
    gg1
  }
  
  cat("==Plotting for model checking: frequent BrS measurements patterns. ==")
  gg1 <- plot_ppd(DIR_list,"case")
  gg0 <- plot_ppd(DIR_list,"control")
  grid.arrange(gg1,gg0,ncol=2)
  
#   case_res_list <- vector("list",length=length(DIR_list))
#   ctrl_res_list <- vector("list",length=length(DIR_list))
#   res_list <- vector("list",length=length(DIR_list))
#   case_pat_list <- vector("list",length=length(DIR_list))
#   ctrl_pat_list <- vector("list",length=length(DIR_list))
#   obs_case_pat_list <- vector("list",length=length(DIR_list))
#   obs_ctrl_pat_list <- vector("list",length=length(DIR_list))
#   
#   for (d in seq_along(DIR_list)){
#       # cases:
#       out_case_pat <- get_top_pattern(out[[d]],1,slice_vec[d],n_pat)
#       case_res_list[[d]] <- out_case_pat$ppd_pat_ct
#       case_res_list[[d]]$DIR <- d 
#       case_res_list[[d]]$ITER <- 1:nrow(case_res_list[[d]])
#       case_res_list[[d]]$CASE <- rep(1,nrow(case_res_list[[d]]))
#       case_res_melt <- reshape2::melt(case_res_list[[d]],
#                                       id.vars = c("CASE","DIR","ITER"),
#                                       variable.name="pattern",
#                                       value.name = "frequency")
#       
#       # controls:
#       out_ctrl_pat <- get_top_pattern(out[[d]],0,slice_vec[d],n_pat)
#       ctrl_res_list[[d]] <- out_ctrl_pat$ppd_pat_ct
#       ctrl_res_list[[d]]$DIR <- d 
#       ctrl_res_list[[d]]$ITER <- 1:nrow(ctrl_res_list[[d]])
#       ctrl_res_list[[d]]$CASE <- rep(0,nrow(ctrl_res_list[[d]]))
#       ctrl_res_melt <- reshape2::melt(ctrl_res_list[[d]],
#                                       id.vars = c("CASE","DIR","ITER"),
#                                       variable.name="pattern",
#                                       value.name = "frequency")
#       
#       res_list[[d]] <- rbind(case_res_melt,ctrl_res_melt)
#       case_pat_list[[d]] <- c(out_case_pat$pattern_names)
#       ctrl_pat_list[[d]] <- c(out_ctrl_pat$pattern_names)
#       
#       obs_case_pat_list[[d]] <- out_case_pat$obs_pat
#       obs_ctrl_pat_list[[d]] <- out_ctrl_pat$obs_pat
#   }
#   
#   if(!(length(unique(case_pat_list))==1)){stop("==Different data sets are used under comparison! Please use the results fitted from the same data set.==")}
#   if(!(length(unique(ctrl_pat_list))==1)){stop("==Different data sets are used under comparison! Please use the results fitted from the same data set.==")}
#   if(!(length(unique(obs_case_pat_list))==1)){stop("==Different data sets are used under comparison! Please use the results fitted from the same data set.==")}
#   if(!(length(unique(obs_ctrl_pat_list))==1)){stop("==Different data sets are used under comparison! Please use the results fitted from the same data set.==")}
#   
#   res <- do.call(rbind,res_list)
#   
#   base_nm <- lapply(DIR_list,basename)
#   NDIR <- length(DIR_list)
#   # first build some functions to summarize posterior distribution 
#   # (following ggplot2 syntax):
#   f <- function(x) {
#     r <- quantile(x, probs = c(0.05, 0.25, 0.5, 0.75, 0.95))
#     names(r) <- c("ymin", "lower", "middle", "upper", "ymax")
#     r
#   }
#   mean_with_nm <- function(x){
#     r <- rep(mean(x),2)
#     names(r)<-c("y","ymax")
#     r
#   }
#   mean_with_nm_txt <- function(x){
#     r <- c(ifelse(max(x)-quantile(x,.97)>0.02,quantile(x,.97)+0.02,max(x)),
#            round(mean(x),3),round(mean(x),3)*100)
#     names(r)<-c("y","ymax","label")
#     r
#   }
#   
#   ## ggplot2:
#   ymax <- max(res$frequency)
#   aes_now <- function(...) {
#     structure(list(...),  class = "uneval")
#   }
#   
#   case_status_labeller <- function(variable,value){
#     c("Case","Control")[2-value]
#   }
#   
#   # plot for cases:
#   hline.data <- as.data.frame(list(frequency=obs_case_pat_list[[1]],
#                               pattern= c(1:(length(obs_case_pat_list[[1]])-1),"other"),
#                               DIR    = rep(1,length(obs_case_pat_list[[1]]))))
#   gg1<-ggplot(data = res[res$CASE==1, ], 
#              aes(x = factor(pattern), y = frequency, fill = factor(DIR))) +
#     #facet_wrap(~ CASE, ncol = 2)+ 
#     #facet_grid(~CASE,labeller=case_status_labeller)+
#     labs(list(x = "pattern", y = "frequency"))+theme_bw()+
#     stat_summary(fun.data = f, geom="boxplot",aes_now(width=dodge_val),
#                  position = position_dodge(dodge_val))+
#     stat_summary(fun.data = mean_with_nm,geom="point",aes(size=1.5),
#                  position = position_dodge(dodge_val))+scale_size(guide = 'none')+
#     stat_summary(fun.data = mean_with_nm_txt,geom="text",
#                  aes(angle=90),position = position_dodge(width = dodge_val))+
#     scale_fill_discrete("Model\n",labels = c(base_nm))+
#     guides(fill=guide_legend(nrow=NDIR,byrow=TRUE))+
#     theme(legend.text = element_text(colour="blue",size = 16, face = "bold"),
#           legend.title = element_text(size=16,face="bold"),legend.position = "top",
#           axis.title   = element_text(size=16,face="bold"),
#           axis.text.x = element_text(angle=40, vjust=.8, hjust=1.01,size=16,face="bold"),
#           strip.text.x = element_text(size = 16, colour = "red",face="bold"))+
#     scale_x_discrete(labels=c(case_pat_list[[1]]))+
#     scale_y_continuous(limits = c(0,ymax))+
#     annotate("text", label = "case", x = length(case_pat_list[[1]])/2, 
#              y = ymax*0.67, size = 8, colour = "red")+
#     geom_errorbar(stat = "hline", 
#                   width=0.8,aes(yintercept = frequency,ymax=..y..,ymin=..y..),
#                   color="blue",size=0.9,data=hline.data)
#   #gg1
#   
#   # plot for controls:
#   hline.data <- as.data.frame(list(frequency=obs_ctrl_pat_list[[1]],
#                                    pattern= c(1:(length(obs_ctrl_pat_list[[1]])-1),"other"),
#                                    DIR    = rep(1,length(obs_ctrl_pat_list[[1]]))))
#   gg0<-ggplot(data = res[res$CASE==0, ], #<-- modified to 0
#              aes(x = factor(pattern), y = frequency, fill = factor(DIR))) +
#     #facet_wrap(~ CASE, ncol = 2)+ 
#     #facet_grid(~CASE,labeller=case_status_labeller)+
#     labs(list(x = "pattern", y = "frequency"))+theme_bw()+
#     stat_summary(fun.data = f, geom="boxplot",aes_now(width=dodge_val),
#                  position = position_dodge(dodge_val))+
#     stat_summary(fun.data = mean_with_nm,geom="point",aes(size=1.5),
#                  position = position_dodge(dodge_val))+scale_size(guide = 'none')+
#     stat_summary(fun.data = mean_with_nm_txt,geom="text",
#                  aes(angle=90),position = position_dodge(width = dodge_val))+
#     scale_fill_discrete("Model\n",labels = c(base_nm))+
#     guides(fill=guide_legend(nrow=NDIR,byrow=TRUE))+
#     theme(legend.text = element_text(colour="blue",size = 16, face = "bold"),
#           legend.title = element_text(size=16,face="bold"),legend.position = "top",
#           axis.title   = element_text(size=16,face="bold"),
#           axis.text.x = element_text(angle=40, vjust=.8, hjust=1.01,size=16,face="bold"),
#           strip.text.x = element_text(size = 16, colour = "red",face="bold"))+
#     scale_x_discrete(labels=c(ctrl_pat_list[[1]]))+# <-- modified to ctrl_pat_list
#     scale_y_continuous(limits = c(0,ymax))+
#     annotate("text", label = "control", x = length(ctrl_pat_list[[1]])/2, 
#              y = ymax*0.67, size = 8, colour = "red")+
#     geom_errorbar(stat = "hline", 
#                   width=0.8,aes(yintercept = frequency,ymax=..y..,ymin=..y..),
#                   color="blue",size=0.9,data=hline.data)
#   #gg0
#   
#   

}
