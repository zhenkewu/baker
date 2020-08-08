#' Plot three-panel figures for nested partially-latent model results
#'
#' `plot_panels()` visualizes the model outputs for communicating how the data inform final
#' latent disease status (etiology). It works for singleton or combo etiologies.
#'  
#' @details Missing data for BrS or SS are dropped when calculating observed measurement
#' positive rates
#'
#' @param DIR_NPLCM File path to the folder containing posterior samples
#' @param slices DEFAULT is "all" - to plot all measurements; Otherwise, one can
#' specify a list: `list(MBS=c(1,3),MSS=1)` means to plot the 1st and
#' 3rd slice of BrS measurements and 1st of SS measurement.
#' @param bg_color A list with names "BrS", "SS", "pie" to specify background colors.
#' The current default is `list(BrS = "lavenderblush", SS = "mistyrose", 
#' pie="antiquewhite")`. If no background is intended, specify as NULL or for a particular
#' measurement, e.g., `BrS = NULL`.
#' @param select_latent a vector of character strings representing latent status. It is used for
#' just plotting a subset of latent status. For example, you can specify `select_latent = "HINF"`
#' to plot all latent status information relevant to `"HINF"`.
#' @param exact Default is `TRUE` to use `select_latent` as exact names of causes. If you want to
#' specify a name and plot all single or combo causes with that name, specify it to be `FALSE`.
#' @param SS_upperlimit The upper limit of horizontal bar for the silver-standard
#' subpanel (the middle panel). The default value is .25.
#'
#' @param eti_upperlimit The upper limit of horizontal bar for the etiology
#' posterior subpanel (the rightmost panel). The default value is .4
#' @param silent Default is TRUE to not print any warning messages; FALSE otherwise.
#' @param ref_eti0 reference quantiles and means; a list: pEti_ref_q, pEti_ref_mean_ord
#' @return A figure with two or three columns
#'
#' @family visualization functions
#' @export

plot_panels <- function(DIR_NPLCM,
                        slices = "all",
                        bg_color = list(BrS = "lavenderblush", 
                                        SS  = "mistyrose",
                                        pie = "antiquewhite"),
                        select_latent = NULL,
                        exact = TRUE,
                        SS_upperlimit=1,
                        eti_upperlimit=1,
                        silent=TRUE,
                        ref_eti0 = NULL){#BEGIN function
  old_par <- graphics::par(no.readonly=TRUE)
  on.exit(graphics::par(old_par))
  
  # read NPLCM outputs:
  out           <- nplcm_read_folder(DIR_NPLCM)
  # organize ouputs:
  Mobs          <- out$Mobs
  Y             <- out$Y
  model_options <- out$model_options
  clean_options <- out$clean_options
  res_nplcm     <- out$res_nplcm
  bugs.dat      <- out$bugs.dat
  rm(out)
  
  #data_nplcm <- list(Mobs  = Mobs, Y = Y)
  data_nplcm <- dget(file.path(DIR_NPLCM,"data_nplcm.txt")) 
  # Determine which three-panel plot to draw:
  parsed_model <- assign_model(model_options, data_nplcm)
  # X not needed in the three-panel plot, but because 'assign_model' was designed
  # to distinguish models even with X, so we have to stick to the useage of 
  # assign_model.
  
  template_BrS <- template_SS <- NULL
  check_combo_BrS <- check_combo_SS <- NULL
  
  if ("BrS" %in% model_options$use_measurements){
    template_BrS <- lapply(clean_options$BrS_objects,"[[","template")
    names(template_BrS) <- lapply(clean_options$BrS_objects,"[[","nm_spec_test")
    check_combo_BrS <- any(unlist(lapply(template_BrS,rowSums))>1)
  }
  if ("SS" %in% model_options$use_measurements){
    template_SS  <- lapply(clean_options$SS_objects,"[[","template")
    names(template_SS) <- lapply(clean_options$SS_objects,"[[","nm_spec_test")
    check_combo_SS <- any(unlist(lapply(template_SS,rowSums))>1)
  }
  if (is.null(template_BrS) && is.null(template_SS)){stop("== No BrS or SS used in the fit. ==")}
  
  #
  # Plot - setup layout for panels:
  #
  
  if (length(slices)==1 && slices=="all") {slices <- lapply(Mobs,seq_along);
                n_total_meas <- sum(parsed_model$num_slice)} # <-- converts slices to a list if it is specified as "all".
  
  n_total_meas <- length(unlist(slices))
  
  it <- graphics::layout(matrix(1:(n_total_meas+2),1,n_total_meas+1+1,byrow = TRUE),
               widths=c(1.5,rep(2.5,n_total_meas),3),heights=c(8))
  ## graphics::layout.show(it)

  # the labels on the left margin:
  height_leftmost <- length(model_options$likelihood$cause_list)
  if (!is.null(select_latent)){
    height_leftmost <- length(select_latent)
    if (!exact){
      height_leftmost <- sum(rowSums(make_template(select_latent,model_options$likelihood$cause_list))>0)
    }    
  }
  
  cat("\n == Plotting Panels of Measurements and Marginal Posterior of Etiology Fractions == \n")
  plot_leftmost(model_options,height_leftmost)
  
  if (!is.null(slices$MBS)){
    # bronze-standard
    for (s in slices$MBS){
      plot_BrS_panel(s,data_nplcm,model_options,
                     clean_options,bugs.dat,res_nplcm,bg_color = bg_color, 
                     select_latent, exact, silent=silent)
    }
  }
  
  if (!is.null(slices$MSS)){
    # silver-standard
    for (s in slices$MSS){
      plot_SS_panel(s,data_nplcm,model_options,
                    clean_options,bugs.dat,res_nplcm,bg_color = bg_color,
                    select_latent,exact)
    }
  }
  plot_pie_panel(model_options,res_nplcm,bugs.dat,bg_color = bg_color,select_latent,exact,ref_eti=ref_eti0)
  
  cat("\n == Done. == \n")
  
#   
#   if (!any(unlist(parsing$reg))){
#     # if no stratification or regression:
#     if (parsing$measurement$quality=="BrS+SS"){
#       if (!parsing$measurement$SSonly){
#         if (!parsing$measurement$nest){
#           print("== BrS+SS; no SSonly; subclass number: K = 1. ==")
#         }else{
#           print("== BrS+SS; no SSonly; subclass number: K > 1. ==")
#         }
#         #
#         # Plot - put layout in place for three panels:
#         #
#         graphics::layout(matrix(c(1,2,3),1,3,byrow = TRUE),
#                widths=c(3,2,3),heights=c(8))
# 
#         # BrS panel:
#         plot_BrS_panel(Mobs$MBS,model_options,clean_options,res_nplcm,
#                              bugs.dat,
#                              top_BrS = 1.3,prior_shape = "interval")
#         # SS panel:
#         plot_SS_panel(Mobs$MSS,model_options,clean_options,res_nplcm,
#                             bugs.dat,top_SS = SS_upperlimit)
#         # Etiology panel:
#         plot_pie_panel(model_options,res_nplcm,bugs.dat,top_pie = eti_upperlimit)
#           
#       } else{# with SS-only measured pathogens:
#         if (!parsing$measurement$nest){
#           stop("== BrS+SS; SSonly; subclass number: K = 1: not done.  ==")
# 
#         }else{
#           stop("== BrS+SS; SSonly; subclass number: K > 1: not done.  ==")
#         }
#       }
#     } else if (parsing$measurement$quality=="BrS"){
#       if (!parsing$measurement$SSonly){
#         if (!parsing$measurement$nest){
#           stop("== BrS; no SSonly; subclass number: K = 1: not done.  ==")
#         }else{
#           stop("== BrS; no SSonly; subclass number: K > 1: not done.  ==")
#         }
#       }
#     }
#   } else{
#     stop("== Three panel plot not implemented for stratification or regression settings. Please check back later for updates. Thanks. ==")
#   }
}# END function

#' order latent status by posterior mean
#' 
#' @param res_nplcm result from model fits
#' @param model_options model specification
#' 
#' @return a list with order (`ord`) and ordered posterior samples (by column)
#' 
#' @export
order_post_eti <- function(res_nplcm,model_options){
  cause_list <- model_options$likelihood$cause_list
  # total no. of causes:
  Jcause     <- length(cause_list)
  # extract and process some data and posterior samples:
  SubVarName <- rep(NA,Jcause)
  for (j in 1:Jcause){
    SubVarName[j] = paste("pEti","[",j,"]",sep="")
  }
  
  # get etiology fraction MCMC samples:
  pEti_mat   <- res_nplcm[,SubVarName,drop=FALSE]
  pEti_mean  <- colMeans(pEti_mat)
  # order the causes by posterior mean:
  ord <- order(pEti_mean)
  pEti_mat_ord <- pEti_mat[,ord]
  make_list(ord,pEti_mat_ord)
}

#' check if a list has elements all of length one
#' 
#' @param x a list
#' 
#' @return TRUE or FALSE
#' @export
#' @examples
#' l = list(a = 5, b = 1:2)
#' is_length_all_one(l) # FALSE
#' l = list(a = 5, b = 1)
#' is_length_all_one(l) # TRUE
#' 
is_length_all_one <- function(x){
   len_vec <-  unlist(lapply(x,length))
   all(len_vec==1)
}

#' get a list of measurement index where to look for data
#' 
#' @param template See [nplcm()]
#' 
#' @return a list of index vectors
#' 
#' @export


get_plot_pos <- function(template){
  pos <- vector("list",nrow(template))
  for (i in 1:nrow(template)){
    if (sum(template[i,])>0){
      pos[[i]] <- which(template[i,] == 1)
    }
  }
  pos
}

#' get the plotting positions (numeric) for the fitted means; 3 positions for each cell
#' 
#' @param e Integer index from 1 to length(cause_list)
#' @param height the total number of causes
#' 
#' @return a triple with numerical plotting positions
#' 
#' @export

get_plot_num <- function(e, height){
  x <- seq(0.5,height+0.5,by=1/4)[-(c(1,(1:height)*4+1))]
  tmp <- length(x)/height*e
  x[c(tmp-2,tmp-1,tmp)]
}

#' plotting the labels on the left margin for panels plot
#' 
#' @param model_options See [nplcm()]
#' @param height no. of rows in the panels plot; commonly set as `length(select_latent)`
#' @return a plot
#' @seealso [plot_panels]
#' @export

plot_leftmost <- function(model_options,height){
  
  op <- graphics::par(mar=c(5.1,4,4.1,0))
  graphics::plot(rep(0,3*height),
       c(sapply(1:height,get_plot_num,height)),
       xlim=c(0,0.1),
       ylim=c(0.5, height+0.5),
       xaxt="n",pch="",xlab="",bty="l",axes=FALSE,
       ylab="",yaxt="n")
  #add axis labels on the left:
  #axis(2,at = c(sapply(1:Jcause,get_plot_num,height=Jcause)),
  #     labels=rep(c("","case","control"),Jcause),las=2)
  # axis(2,at=(1:Jcause)-.45,labels=rep("",Jcause),las=2,cex.axis=.5)
  # axis(2,at=(1:Jcause)-.35,labels=rep("",Jcause),las=2,cex.axis=.5)
  
  graphics::text(0.1,c(sapply(1:height,get_plot_num,height)),
       labels=rep(c(expression(paste(symbol("\052"),"(FPR)--",Delta,"(fitted)--+(TPR)")),
                    "case","control"),height),adj=1,cex=c(1,2,2),
       col=c("purple",1,1))
  graphics::text(0.1,c(sapply(1:height,get_plot_num,height))+0.1,
       labels=rep(c(expression(italic("posterior mean:")),"",
                    expression(italic("data:"))),height),col=c("purple",1,1),adj=1)
  graphics::text(0.1,(1:height)-.35,labels=rep("posterior CI: '|'-95%;'[]'-50%",height),col="purple",adj=1)
  graphics::text(0.1,(1:height)-.45,labels=rep("prior: '|'-95%;'[]'-50%",height),adj=1)
}
