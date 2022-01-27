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
#' @param is_plot default to `TRUE` for plotting only; set to `FALSE` if to get summary.
#' @return A figure with two or three columns (if `is_plot=TRUE`); otherwise, it
#' provide posterior summaries of Etiology information to 
#' used by [print.summary.nplcm.no_reg()]
#'
#' @family visualization functions
#'         
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
                        ref_eti0 = NULL,is_plot=TRUE){#BEGIN function
  if(is_plot){
    old_par <- graphics::par(no.readonly=TRUE)
    on.exit(graphics::par(old_par))
  }
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
  
  if (is_plot){
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
    
    cat("== Plotting Panels of Measurements and Marginal Posterior of CSCFs == \n")
    plot_leftmost(model_options,height_leftmost)
    
    if (!is.null(slices$MBS)){
      res_BrS <- vector("list",length(slices$MBS))
      # bronze-standard
      l <- 0
      for (s in slices$MBS){
        l <- l+1
        res_BrS[[l]] <- plot_BrS_panel(s,data_nplcm,model_options,
                                       clean_options,bugs.dat,res_nplcm,bg_color = bg_color, 
                                       select_latent, exact, silent=silent)
      }
    }
    
    if (!is.null(slices$MSS)){
      res_SS <- vector("list",length(slices$MSS))
      # silver-standard
      l <- 0
      for (s in slices$MSS){
        l <- l+1
        res_SS[[l]] <- plot_SS_panel(s,data_nplcm,model_options,
                                     clean_options,bugs.dat,res_nplcm,bg_color = bg_color,
                                     select_latent,exact)
      }
    }
  }
  res_Eti <- plot_pie_panel(model_options,res_nplcm,bugs.dat,bg_color = bg_color,
                            select_latent,exact,ref_eti=ref_eti0,is_plot=is_plot)
  
  if(is_plot){cat("== Done. == \n")}
  if (!is_plot){return(make_list(res_Eti,parsed_model))}
  
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
#' @examples
#' l = list(a = 5, b = 1:2)
#' baker:::is_length_all_one(l) # FALSE
#' l = list(a = 5, b = 1)
#' baker:::is_length_all_one(l) # TRUE
#' @export       
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









#' Plot bronze-standard (BrS) panel
#' 
#' 
#' @param slice the index of measurement slice for BrS.
#' @param data_nplcm See [nplcm()]
#' @param model_options See [nplcm()]
#' @param clean_options See [clean_perch_data()]
#' @param res_nplcm See [nplcm_read_folder()]
#' @param bugs.dat Data input for the model fitting.
#' @param bg_color A list with names "BrS", "SS", "pie" to specify background colors
#' @param select_latent a vector of character strings representing latent status. It is used for
#' just plotting a subset of latent status. For example, you can specify `select_latent = "HINF"`
#' to plot all latent status information relevant to `"HINF"`.
#' @param exact Default is `TRUE` to use `select_latent` as exact names of causes. If you want to
#' specify a name and plot all single or combo causes with that name, specify it to be `FALSE`.
#' @param top_BrS Numerical value to specify the rightmost limit 
#' on the horizontal axis for the BrS panel.
#' @param cexval Default is 1 - size of text of the BrS percentages.
#' @param srtval Default is 0 - the direction of the text for the BrS percentages.
#' @param prior_shape `interval` or `boxplot` - for how to represent
#' prior/posteriors of the TPR/FPRs of measurements.
#' @param silent Default is TRUE to not print any warning messages; FALSE otherwise.
#' 
#' @return plotting function.
#' 
#' @family visualization functions
plot_BrS_panel <- function(slice,data_nplcm,model_options,
                           clean_options,bugs.dat,res_nplcm,
                           bg_color,
                           select_latent = NULL,
                           exact   = TRUE,
                           top_BrS = 1.3, 
                           cexval = 1,
                           srtval = 0,
                           prior_shape="interval",
                           silent=TRUE){
  template_BrS    <- NULL
  check_combo_BrS <- NULL
  if ("BrS" %in% model_options$use_measurements){
    template_BrS <- lapply(clean_options$BrS_objects,"[[","template")
    names(template_BrS) <- lapply(clean_options$BrS_objects,"[[","nm_spec_test")
    check_combo_BrS <- any(unlist(lapply(template_BrS,rowSums))>1)
  }
  # order cause_list by posterior means:
  ord <- order_post_eti(res_nplcm,model_options)$ord
  pEti_mat_ord <- order_post_eti(res_nplcm,model_options)$pEti_mat_ord
  
  template_ord <- template_BrS[[slice]][ord,,drop=FALSE]
  cause_list <- model_options$likelihood$cause_list
  cause_list_ord <- cause_list[ord]
  
  # focus on selected latent status:
  latent_seq <- get_latent_seq(cause_list,ord,select_latent,exact)$latent_seq
  template_ord <- template_ord[latent_seq,,drop = FALSE]
  
  
  # which model was fitted:
  parsed_model <- assign_model(model_options, data_nplcm)
  this_slice_nest <- parsed_model$nested[slice]
  
  if (!any(unlist(parsed_model$regression))){
    if (!this_slice_nest){ marginal_rates <- get_marginal_rates_no_nested(slice,res_nplcm,model_options,data_nplcm)}
    if (this_slice_nest){marginal_rates <- get_marginal_rates_nested(slice,res_nplcm,model_options,data_nplcm)}
    
    # TPR and FPR rates:
    theta_mat <- as.matrix(marginal_rates$res_tpr)
    theta_mean <- colMeans(theta_mat)
    
    psi_mat   <- as.matrix(marginal_rates$res_fpr)
    psi_mean   <- colMeans(psi_mat)
    
    ## model fitted postive rate for each pathogen
    if (!this_slice_nest){fitted_mean <- get_fitted_mean_no_nested(slice,res_nplcm,model_options,data_nplcm,clean_options)}
    if (this_slice_nest){fitted_mean <- get_fitted_mean_nested(slice,res_nplcm,model_options,data_nplcm,clean_options)}
    fittedmean_case    <- fitted_mean$res_case
    fittedmean_ctrl    <- fitted_mean$res_ctrl
    #plot_pos <- get_plot_pos(template_ord)  # 1 at the 5th means for the fifth smallest etiology, we should plot 1st dimension in this slice.
  }
  
  # get observed rates' summaries:
  Nd <- bugs.dat$Nd
  Nu <- bugs.dat$Nu
  
  MBS_curr <- data_nplcm$Mobs$MBS[[slice]]
  
  # positive rates and confidence intervals:
  #cases:
  
  MBS_case_curr <- MBS_curr[1:Nd,,drop=FALSE]
  count    <- as.integer(do.call(cbind,lapply(MBS_case_curr,sum,na.rm=TRUE)))
  NA_count <- apply(MBS_case_curr,2,function(v) sum(is.na(v)))
  tmp.case <- binom::binom.confint(count,Nd-NA_count,conf.level = 0.95, methods = "ac")
  
  #controls:
  MBS_ctrl_curr <- MBS_curr[-(1:Nd),,drop=FALSE]
  count    <- as.integer(do.call(cbind,lapply(MBS_ctrl_curr,sum,na.rm=TRUE)))
  NA_count <- apply(MBS_ctrl_curr,2,function(v) sum(is.na(v)))
  tmp.ctrl <- binom::binom.confint(count, Nu-NA_count, conf.level = 0.95, methods = "ac")
  
  # case and control positive rate, lower and upper limit
  MBS_mean  <- rbind(round(tmp.case$mean,5),round(tmp.ctrl$mean,5))
  MBS_q1 <- rbind(tmp.case[,c("lower")],tmp.ctrl[,c("lower")])
  MBS_q2 <- rbind(tmp.case[,c("upper")],tmp.ctrl[,c("upper")])
  
  # prior parameters:
  alphaB_nm   <- paste0("alphaB_",slice)
  betaB_nm   <- paste0("betaB_",slice)
  alphaB         <- bugs.dat[[alphaB_nm]]
  betaB          <- bugs.dat[[betaB_nm]]
  
  pos_vec <- get_plot_pos(template_ord)
  
  get_COR <- function(pos){
    if (is.null(pos) || is.na(pos)){return(NULL)}
    y <- c(rep(1,Nd), rep(0,Nu))
    brs.data <- as.data.frame(rbind(MBS_case_curr,MBS_ctrl_curr))
    dat.reg  <- as.data.frame(cbind(y,brs.data))
    fit      <- stats::glm(y~.,data=dat.reg,family=stats::binomial,na.action="na.omit")
    
    if (sum(is.na(stats::coef(fit)))==0 & sum(diag(stats::vcov(fit))>100)==0){
      res0 = cbind(exp(suppressMessages(stats::confint(fit))),exp(fit$coef))[-1,]
      res  = list(ORinterval = matrix(res0,nrow = ncol(brs.data),ncol=3)[pos,],
                  label = "conditional OR")
    }else{
      if (!silent){
        print("==[baker] Conditional OR not calculatable. Switch to mariginal OR.\n")
      }
      res0 = matrix(NA,nrow=1,ncol=3)
      res0 = data.frame(res0)
      tb <- table(dat.reg$y,brs.data[,pos])
      fit_tmp  <- stats::glm(y~brs.data[,pos],family=stats::binomial,na.action = "na.omit")
      if (length(stats::vcov(fit_tmp))>1 && stats::vcov(fit_tmp)[2,2]<100 && ncol(tb)==2){
        #print(l)
        res0 <- cbind(exp(suppressMessages(stats::confint(fit_tmp))),exp(fit_tmp$coef))[-1,]
      }
      res = list(ORinterval=res0,label="marginal OR")
    }
    res
  }
  
  plot_BrS_cell  <- function(lat_pos, pos, height,gap = 0){
    plotat <- get_plot_num(lat_pos,height) + gap
    
    graphics::plot(c(fittedmean_case[pos],MBS_mean[,pos]),
                   plotat,
                   xlim=c(0,top_BrS),
                   ylim=c(0.5, height+0.5),
                   xaxt="n",xlab="positive rate",
                   ylab="",yaxt="n",
                   pch = c(2,20,20),
                   col = c("purple","dodgerblue2", "dodgerblue2"),
                   cex = c(1,2,2))
    
    
    graphics::points(c(theta_mean[pos],MBS_q2[,pos]),
                     plotat,
                     pch = c("+","|","|"),
                     col = c("purple",1,1),
                     cex = c(2,1,1))
    graphics::points(c(fittedmean_ctrl[pos],MBS_q1[,pos]),
                     plotat,
                     pch = c("*","|","|"),
                     col = c("purple",1,1),
                     cex = c(2,1,1))
    # connect case and control rates:
    graphics::segments(
      x0 = MBS_mean[1,pos],x1 = MBS_mean[2,pos],
      y0 = plotat[2],y1 = plotat[3],
      lty = 1,
      col = "dodgerblue2",
      lwd = 2
    )
    # case: rates
    graphics::segments(
      x0 = MBS_q1[1,pos],x1 = MBS_q2[1,pos],
      y0 = plotat[2], y1 = plotat[2],
      lty = 1
    )
    tmp.hpos <- ifelse(MBS_q2[1,pos]+0.15>0.95,MBS_q1[1,pos]-0.2,MBS_q2[1,pos]+0.15 )
    graphics::text(tmp.hpos, plotat[2], paste0(round(100*MBS_mean[1,pos],1),"%"),
                   srt=srtval,cex=cexval)
    # control:rates
    graphics::segments(
      x0 = MBS_q1[2,pos],x1 = MBS_q2[2,pos],
      y0 = plotat[3], y1 = plotat[3],
      lty = 1
    )
    tmp.hpos <- ifelse(MBS_q2[2,pos]+0.15>0.95,MBS_q1[2,pos]-0.2,MBS_q2[2,pos]+0.15 )
    graphics::text(tmp.hpos, plotat[3], paste0(round(100*MBS_mean[2,pos],1),"%"),
                   srt=srtval,cex=cexval)
    # poster means of TPR, FPR and fitted marginal rate:
    graphics::segments(
      x0 = theta_mean[pos],x1 = psi_mean[pos],
      y0 = plotat[1], y1 = plotat[1],
      lty = 4,col="gray"
    )
    
    if (!is.null(pos) && !is.na(pos)){#some pos can be NA: because certain cause has no measurements.
      # prior and posterior of TPR:
      if (prior_shape == "interval") {
        # prior of TPR:
        prior_plot_at <- lat_pos - .45 + gap
        tmp = stats::qbeta(c(0.025,0.975,0.25,0.75),alphaB[pos],betaB[pos])
        graphics::points(tmp,rep(prior_plot_at,4),pch = c("|","|","[","]"),col="gray")
        graphics::segments(tmp[1],prior_plot_at,
                           tmp[2],prior_plot_at,lty = 1,col="gray")
        graphics::segments(tmp[3],prior_plot_at,
                           tmp[4],prior_plot_at,lty = 1,col="gray",lwd=2)
        
        # posterior of TPR:
        post_plot_at <- lat_pos - .35 + gap
        tmp.post = as.matrix(theta_mat)[,pos]
        tmp  = stats::quantile(tmp.post, c(0.025,0.975,0.25,0.75))
        graphics::points(tmp,rep(post_plot_at,4),pch = c("|","|","[","]"),col = "purple")
        graphics::segments(tmp[1],post_plot_at,
                           tmp[2],post_plot_at,lty = 1)
        graphics::segments(tmp[3],post_plot_at,
                           tmp[4],post_plot_at,lty = 1,lwd=2)
      } else if (prior_shape == "boxplot") {
        tmp = stats::rbeta(10000,alphaB[pos],betaB[pos])
        graphics::boxplot(
          tmp,at = prior_plot_at, boxwex = 1 / 10 , col = "gray",
          add = TRUE,horizontal = TRUE,outline = FALSE,xaxt =
            "n"
        )
        tmp.post = as.matrix(theta_mat)[,pos]
        graphics::boxplot(
          tmp.post,at = post_plot_at,boxwex = 1 / 10,add = TRUE,
          horizontal = TRUE,outline = FALSE,xaxt = "n"
        )
      }
      
      #plot conditional odds ratio on the right:
      tmp0 <- get_COR(pos)
      tmp  <- as.numeric(tmp0$ORinterval) # <--- used as.numeric() here because if tmp0 = NAs, then the code will fail.
      
      L <-  round(tmp[1],1)
      C <-  round(tmp[3],1)
      R <-  round(tmp[2],1)
      
      graphics::text(top_BrS - 0.12,lat_pos + .3, colnames(MBS_case_curr)[pos],cex=1 )
      graphics::text(top_BrS - 0.12,lat_pos + 1 / (2 * Jcause),C,cex = 1.5)
      graphics::text(top_BrS - 0.12,lat_pos - .2,paste(c(L,"   ",R),collapse = " "),
                     cex = 1.2)
      graphics::legend("topright",tmp0$label,bty = "n")
    }
    
    # x-axis for each cell:
    if (lat_pos>1){
      graphics::axis(1, seq(0,1,by = .1), lwd = 0, lwd.ticks = 0,#labels=rep("",length(seq(0,1,by=.2))),
                     pos = seq(.6, height +.6,by = 1)[lat_pos], cex.axis = 0.8,
                     lty = 2,col = "blue"
      )
    }
  }
  
  points_BrS_cell     <- function(lat_pos,pos,height,gap=0){ # pos for the measurement dimension, usually used as pos_vec[e].
    plotat <- get_plot_num(lat_pos,height) + gap
    
    graphics::points(c(fittedmean_case[pos],MBS_mean[,pos]),
                     plotat,
                     xlim=c(0,top_BrS),
                     ylim=c(0.5, height+0.5),
                     xaxt="n",xlab="positive rate",
                     ylab="",yaxt="n",
                     pch = c(2,20,20),
                     col = c("purple","dodgerblue2", "dodgerblue2"),
                     cex = c(1,2,2))
    
    
    graphics::points(c(theta_mean[pos],MBS_q2[,pos]),
                     plotat,
                     pch = c("+","|","|"),
                     col = c("purple",1,1),
                     cex = c(2,1,1))
    graphics::points(c(fittedmean_ctrl[pos],MBS_q1[,pos]),
                     plotat,
                     pch = c("*","|","|"),
                     col = c("purple",1,1),
                     cex = c(2,1,1))
    # connect case and control rates:
    graphics::segments(
      x0 = MBS_mean[1,pos],x1 = MBS_mean[2,pos],
      y0 = plotat[2],y1 = plotat[3],
      lty = 1,
      col = "dodgerblue2",
      lwd = 2
    )
    # case: rates
    graphics::segments(
      x0 = MBS_q1[1,pos],x1 = MBS_q2[1,pos],
      y0 = plotat[2], y1 = plotat[2],
      lty = 1
    )
    tmp.hpos <- ifelse(MBS_q2[1,pos]+0.15>0.95,MBS_q1[1,pos]-0.2,MBS_q2[1,pos]+0.15 )
    graphics::text(tmp.hpos, plotat[2], paste0(round(100*MBS_mean[1,pos],1),"%"),
                   srt=srtval,cex=cexval)
    # control:rates
    graphics::segments(
      x0 = MBS_q1[2,pos],x1 = MBS_q2[2,pos],
      y0 = plotat[3], y1 = plotat[3],
      lty = 1
    )
    tmp.hpos <- ifelse(MBS_q2[2,pos]+0.15>0.95,MBS_q1[2,pos]-0.2,MBS_q2[2,pos]+0.15 )
    graphics::text(tmp.hpos, plotat[3], paste0(round(100*MBS_mean[2,pos],1),"%"),
                   srt=srtval,cex=cexval)
    # poster means of TPR, FPR and fitted marginal rate:
    graphics::segments(
      x0 = theta_mean[pos],x1 = psi_mean[pos],
      y0 = plotat[1], y1 = plotat[1],
      lty = 4,col="gray"
    )
    
    if (!is.null(pos) && !is.na(pos)){#some pos can be NA: because certain cause has no measurements.
      # prior and posterior of TPR:
      if (prior_shape == "interval") {
        # prior of TPR:
        prior_plot_at <- lat_pos - .45 + gap
        tmp = stats::qbeta(c(0.025,0.975,0.25,0.75),alphaB[pos],betaB[pos])
        graphics::points(tmp,rep(prior_plot_at,4),pch = c("|","|","[","]"),col="gray")
        graphics::segments(tmp[1],prior_plot_at,
                           tmp[2],prior_plot_at,lty = 1,col="gray")
        graphics::segments(tmp[3],prior_plot_at,
                           tmp[4],prior_plot_at,lty = 1,col="gray",lwd=2)
        
        # posterior of TPR:
        post_plot_at <- lat_pos - .35 + gap
        tmp.post = as.matrix(theta_mat)[,pos]
        tmp  = stats::quantile(tmp.post, c(0.025,0.975,0.25,0.75))
        graphics::points(tmp,rep(post_plot_at,4),pch = c("|","|","[","]"),col = "purple")
        graphics::segments(tmp[1],post_plot_at,
                           tmp[2],post_plot_at,lty = 1)
        graphics::segments(tmp[3],post_plot_at,
                           tmp[4],post_plot_at,lty = 1,lwd=2)
      } else if (prior_shape == "boxplot") {
        tmp = stats::rbeta(10000,alphaB[pos],betaB[pos])
        graphics::boxplot(
          tmp,at = prior_plot_at, boxwex = 1 / 10 , col = "gray",
          add = TRUE,horizontal = TRUE,outline = FALSE,xaxt =
            "n"
        )
        tmp.post = as.matrix(theta_mat)[,pos]
        graphics::boxplot(
          tmp.post,at = post_plot_at,boxwex = 1 / 10,add = TRUE,
          horizontal = TRUE,outline = FALSE,xaxt = "n"
        )
      }
      
      #plot conditional odds ratio on the right:
      tmp0 <- get_COR(pos)
      tmp  <- as.numeric(tmp0$ORinterval)
      
      L <-  round(tmp[1],1)
      C <-  round(tmp[3],1)
      R <-  round(tmp[2],1)
      
      graphics::text(top_BrS - 0.12,lat_pos + .3+gap, colnames(MBS_case_curr)[pos],cex=1 )
      graphics::text(top_BrS - 0.12,lat_pos + 1 / (2 * Jcause)+gap,C,cex = 1.5)
      graphics::text(top_BrS - 0.12,lat_pos - .2+gap,paste(c(L,"   ",R),collapse = " "),
                     cex = 1.2)
      graphics::legend("topright",tmp0$label,bty = "n")
    }
    
    # x-axis for each cell:
    if (lat_pos>1){
      graphics::axis(1, seq(0,1,by = .1), lwd = 0, lwd.ticks = 0,#labels=rep("",length(seq(0,1,by=.2))),
                     pos = seq(.6,height +.6,by = 1)[lat_pos], cex.axis = 0.8,
                     lty = 2,col = "blue"
      )
    }
  }
  
  Jcause <- length(model_options$likelihood$cause_list)
  #
  # plotting:
  #
  #op <- graphics::par(mar=c(5.1,4.1,4.1,0))
  op <- graphics::par(mar=c(5.1,0,4.1,0))
  
  if (!is_length_all_one(pos_vec)){
    #stop("== Not implemented for combo latent status.==")
    warning("==[baker] Combo latent status implemented with measurements overlapping in BrS columns! ==\n")
  }
  
  first  <- TRUE
  cat("== Plotting BrS Slice: ", slice, ": ", unlist(names(data_nplcm$Mobs$MBS))[slice])
  for (e in 1:nrow(template_ord)){
    gap_seq <- 0
    if (!is.null(pos_vec[[e]]) && length(pos_vec[[e]])>1){
      gap_seq <- c(0,-0.1*(1:(length(pos_vec[[e]])-1)))
    }
    ct <- 0
    for (pos_curr in pos_vec[[e]]){
      if (!is.na(pos_curr)){
        ct <- ct +1
        if (first) {plot_BrS_cell(e,pos_curr, length(latent_seq),gap = gap_seq[ct]);first <- FALSE}
        if (!first) {points_BrS_cell(e,pos_curr, length(latent_seq),gap=gap_seq[ct])}
      }
    }
  }
  
  if (!is.null(bg_color) && !is.null(bg_color$BrS)){
    graphics::rect(graphics::par("usr")[1], graphics::par("usr")[3], graphics::par("usr")[2], graphics::par("usr")[4], col = 
                     bg_color$BrS)
    
    for (e in 1:nrow(template_ord)){
      gap_seq <- 0
      if (!is.null(pos_vec[[e]]) && length(pos_vec[[e]])>1){
        gap_seq <- c(0,-0.1*(1:(length(pos_vec[[e]])-1)))
      }
      
      ct <- 0
      for (pos_curr in pos_vec[[e]]){
        if (!is.na(pos_curr)){
          ct <- ct +1
          points_BrS_cell(e,pos_curr,length(latent_seq),gap=gap_seq[ct])
        }
      }
    }
  }
  
  #   #add graphics::axis labels on the left:
  #   graphics::axis(2,at = c(sapply(1:Jcause,get_plot_num,height=Jcause)),
  #        labels=rep(c("","case","ctrl"),Jcause),las=2)
  #   graphics::axis(2,at=(1:Jcause)-.45,labels=rep("",Jcause),las=2,cex.axis=.5)
  #   graphics::axis(2,at=(1:Jcause)-.35,labels=rep("",Jcause),las=2,cex.axis=.5)
  #   
  
  
  if (sum(template_ord)==0){
    warning(paste0("==[baker] Bronze-standard slice ", names(data_nplcm$Mobs$MBS)[slice], 
                   " has no measurements informative of the causes! Please check if measurements' columns correspond to causes.==\n"))  
    plotat <- c(sapply(seq_along(latent_seq),get_plot_num,length(latent_seq)))
    graphics::plot(rep(0,length(plotat)),
                   plotat,
                   xlim=c(0,top_BrS),
                   ylim=c(0.5, length(latent_seq)+0.5),
                   xaxt="n",xlab="positive rate",
                   ylab="",yaxt="n",
                   pch = c("","",""))
  }
  #add ticks from 0 to 1 for x-bar:
  graphics::axis(1,at = c(0,0.2,0.4,0.6,0.8,1),labels= as.character(c(0,0.2,0.4,0.6,0.8,1)),las=1)
  
  #add dashed lines to separate cells:
  if (length(latent_seq) > 1){
    graphics::abline(h=seq(1.5,length(latent_seq)-.5,by=1),lty=2,lwd=0.5,col="gray")
  }
  graphics::abline(v=1,lty=2,lwd=.5,col="gray")
  
  #add some texts:
  graphics::mtext(eval(paste0("BrS: ", names(data_nplcm$Mobs$MBS)[slice])),
                  line=1,cex=1.8)
}


#' get model fitted mean for conditional independence model
#' 
#' @inheritParams get_fitted_mean_nested
#' @return a list with model fitted means
get_fitted_mean_no_nested <- function(slice,res_nplcm,model_options,data_nplcm,
                                      clean_options){
  # order cause_list by posterior means:
  ord <- order_post_eti(res_nplcm,model_options)$ord
  pEti_mat_ord <- order_post_eti(res_nplcm,model_options)$pEti_mat_ord
  template_BrS <- lapply(clean_options$BrS_objects,"[[","template")
  names(template_BrS) <- lapply(clean_options$BrS_objects,"[[","nm_spec_test")
  template_ord <- template_BrS[[slice]][ord,,drop=FALSE]
  
  JBrS_curr <- ncol(data_nplcm$Mobs$MBS[[slice]])
  Jcause    <- length(model_options$likelihood$cause_list)
  
  thetaBS_nm <- paste0("^thetaBS_",slice)
  psiBS_nm   <- paste0("^psiBS_",slice)
  
  #posterior distribution of FPR:
  theta_mat    <- res_nplcm[,grep(thetaBS_nm,colnames(res_nplcm)),drop=FALSE]
  #posterior distribution of FPR:
  psi_mat    <- res_nplcm[,grep(psiBS_nm,colnames(res_nplcm)),drop=FALSE]
  
  ## model fitted postive rate for each pathogen
  fittedmean_case    <- rep(NA,JBrS_curr)
  names(fittedmean_case) <- colnames(data_nplcm$Mobs$MBS[[slice]])
  
  fitted_margin_case <- function(pEti_ord,theta,psi,template){
    mixture <-  pEti_ord
    tpr     <-  t(t(template)*theta)
    fpr     <- t(t(1-template)*psi)
    colSums(tpr*mixture + fpr*mixture)
  }
  
  res_case  <- colMeans(t(sapply(1:nrow(pEti_mat_ord),
                                 function(iter)
                                   fitted_margin_case(pEti_mat_ord[iter,], theta_mat[iter,],
                                                      psi_mat[iter,],template_ord))))
  res_case <- as.matrix(res_case)
  res_ctrl <- colMeans(as.matrix(psi_mat))
  
  make_list(res_case,res_ctrl)
}

#' get fitted mean for nested model with subclass mixing weights that are the same among cases
#' 
#' @param slice the slice of BrS data that are modeled
#' @param res_nplcm matrix of MCMC samples
#' @param model_options see [nplcm()]
#' @param data_nplcm see [nplcm()]
#' @param clean_options see [clean_perch_data()]
#' @return a matrix of no. of rows equal to retained MCMC samples, no. of columns
#' equal to the no. of measurement dimensions within a slice.
#' 
get_fitted_mean_nested <- function(slice,res_nplcm, model_options,
                                   data_nplcm,clean_options){
  # order cause_list by posterior means:
  ord <- order_post_eti(res_nplcm,model_options)$ord
  pEti_mat_ord <- order_post_eti(res_nplcm,model_options)$pEti_mat_ord
  template_BrS <- lapply(clean_options$BrS_objects,"[[","template")
  names(template_BrS) <- lapply(clean_options$BrS_objects,"[[","nm_spec_test")
  template_ord <- template_BrS[[slice]][ord,,drop=FALSE]
  
  JBrS_curr <- ncol(data_nplcm$Mobs$MBS[[slice]])
  Jcause    <- length(model_options$likelihood$cause_list)
  K_curr    <- model_options$likelihood$k_subclass[slice]  
  
  res_case <- matrix(NA,nrow = nrow(res_nplcm),ncol=JBrS_curr)
  res_ctrl <- matrix(0,nrow = nrow(res_nplcm),ncol=JBrS_curr)
  # formula: fittedmean_case[j] = \sum_e pEti[e]*\sum_k {Eta[e,k]*Theta[j,k]*templateBS[e,j]+Eta[e,k]*Psi[j,k]*(1-templateBS[e,j])}
  for (j in 1:JBrS_curr){
    # get ThetaBS[j,k]:
    ind_ThetaBS_tmp <- grep(paste0("^ThetaBS_",slice,"\\[",j,","),colnames(res_nplcm))
    if (length(ind_ThetaBS_tmp)!=K_curr){stop("== Check `ThetaBS` extraction from posterior samples! No. of subclasses not matched with specification.==")}
    ThetaBS_tmp <- res_nplcm[,ind_ThetaBS_tmp]
    # get PsiBS[j,k]:
    ind_PsiBS_tmp <- grep(paste0("^PsiBS_",slice,"\\[",j,","),colnames(res_nplcm))
    if (length(ind_PsiBS_tmp)!=K_curr){stop("== Check `PsiBS` extraction from posterior samples! No. of subclasses not matched with specification.==")}
    PsiBS_tmp <- res_nplcm[,ind_PsiBS_tmp]
    # get Eta[e,k]:
    ind_Eta_tmp <- grep(paste0("^Eta_",slice),colnames(res_nplcm))
    if (length(ind_Eta_tmp)!=K_curr){stop("== Check `Eta` extraction from posterior samples! No. of subclasses not matched with specification.==")}
    Eta_tmp <- res_nplcm[,ind_Eta_tmp]
    
    term_e <- matrix(0,nrow = nrow(res_nplcm),ncol=Jcause)
    
    for (e in 1:Jcause){
      # get templateBS[e,j]:
      indBS <- template_ord[e,j]
      # get pEti[e]:
      ind_pEti_tmp <- grep(paste0("^pEti\\[",ord[e],"\\]$"),colnames(res_nplcm))
      if (length(ind_pEti_tmp)!=1){stop("== Error in extracting etiology! ==")}
      pEti_tmp <- res_nplcm[,ind_pEti_tmp]
      
      # calculate by formula:
      for (k in 1:K_curr){
        if (indBS == 1){term_e[,e] <- term_e[,e] + Eta_tmp[,k]*ThetaBS_tmp[,k]}
        if (indBS == 0){term_e[,e] <- term_e[,e] + Eta_tmp[,k]*PsiBS_tmp[,k]}
      }
      term_e[,e] <- term_e[,e]*pEti_tmp
    }
    res_case[,j] <- rowSums(term_e)
    ## now handle controls:
    # get Lambda[k]:
    ind_Lambda_tmp <- grep(paste0("^Lambda_",slice),colnames(res_nplcm))
    if (length(ind_Lambda_tmp)!=K_curr){stop("== Check `Lambda` extraction from posterior samples! No. of subclasses not matched with specification.==")}
    Lambda_tmp <- res_nplcm[,ind_Lambda_tmp]
    
    for (k in 1:K_curr){
      res_ctrl[,j] <- res_ctrl[,j] + Lambda_tmp[,k]*PsiBS_tmp[,k]
    }
  }
  res_case <- colMeans(as.matrix(res_case))
  res_ctrl <- colMeans(as.matrix(res_ctrl))
  make_list(res_case,res_ctrl)
} 

#' get marginal TPR and FPR for no nested model
#' 
#' @param slice the slice of BrS data that are modeled
#' @param res_nplcm matrix of MCMC samples
#' @param model_options see [nplcm()]
#' @param data_nplcm see [nplcm()]
#' 
#' @return a matrix of no. of rows equal to retained MCMC samples, no. of columns
#' equal to the no. of measurement dimensions within a slice.
#' 
get_marginal_rates_no_nested <- function(slice, res_nplcm, model_options,data_nplcm){
  
  thetaBS_nm <- paste0("^thetaBS_",slice)
  psiBS_nm   <- paste0("^psiBS_",slice)
  
  res_tpr  <- res_nplcm[,grep(thetaBS_nm,colnames(res_nplcm)),drop=FALSE]
  
  #posterior distribution of FPR:
  res_fpr  <- res_nplcm[,grep(psiBS_nm,colnames(res_nplcm)),drop=FALSE]
  
  make_list(res_tpr,res_fpr)
}

#' get marginal TPR and FPR for nested model
#' 
#' @param slice the slice of BrS data that are modeled
#' @param res_nplcm matrix of MCMC samples
#' @param model_options see [nplcm()]
#' @param data_nplcm see [nplcm()]
#' 
#' @return a matrix of no. of rows equal to retained MCMC samples, no. of columns
#' equal to the no. of measurement dimensions within a slice.
#' 
get_marginal_rates_nested <- function(slice, res_nplcm, model_options,data_nplcm){
  JBrS_curr <- ncol(data_nplcm$Mobs$MBS[[slice]])
  Jcause    <- length(model_options$likelihood$cause_list)
  K_curr    <- model_options$likelihood$k_subclass[slice]  
  
  res_tpr <- matrix(0,nrow = nrow(res_nplcm),ncol=JBrS_curr)
  res_fpr <- matrix(0,nrow = nrow(res_nplcm),ncol=JBrS_curr)
  # formula: fittedmean_case[j] = \sum_e pEti[e]*\sum_k {Eta[e,k]*Theta[j,k]*templateBS[e,j]+Eta[e,k]*Psi[j,k]*(1-templateBS[e,j])}
  for (j in 1:JBrS_curr){
    # get ThetaBS[j,k]:
    ind_ThetaBS_tmp <- grep(paste0("^ThetaBS_",slice,"\\[",j,","),colnames(res_nplcm))
    if (length(ind_ThetaBS_tmp)!=K_curr){stop("==[baker] Check `ThetaBS` extraction from posterior samples! No. of subclasses not matched with specification.==")}
    ThetaBS_tmp <- res_nplcm[,ind_ThetaBS_tmp]
    # get PsiBS[j,k]:
    ind_PsiBS_tmp <- grep(paste0("^PsiBS_",slice,"\\[",j,","),colnames(res_nplcm))
    if (length(ind_PsiBS_tmp)!=K_curr){stop("==[baker] Check `PsiBS` extraction from posterior samples! No. of subclasses not matched with specification.==")}
    PsiBS_tmp <- res_nplcm[,ind_PsiBS_tmp]
    # get Eta[e,k]:
    ind_Eta_tmp <- grep(paste0("^Eta_",slice),colnames(res_nplcm))
    if (length(ind_Eta_tmp)!=K_curr){stop("==[baker] Check `Eta` extraction from posterior samples! No. of subclasses not matched with specification.==")}
    Eta_tmp <- res_nplcm[,ind_Eta_tmp]
    # get Lambda[k]:
    ind_Lambda_tmp <- grep(paste0("^Lambda_",slice),colnames(res_nplcm))
    if (length(ind_Lambda_tmp)!=K_curr){stop("==[baker] Check `Lambda` extraction from posterior samples! No. of subclasses not matched with specification.==")}
    Lambda_tmp <- res_nplcm[,ind_Lambda_tmp]
    
    # calculate by formula:
    for (k in 1:K_curr){
      res_tpr[,j] <- res_tpr[,j] + Eta_tmp[,k]*ThetaBS_tmp[,k]
      res_fpr[,j] <- res_fpr[,j] + Lambda_tmp[,k]*PsiBS_tmp[,k] # <-- may need to change Lambda into Eta if we mean FPR among cases.
    }
  }
  make_list(res_tpr,res_fpr)
}


### get marginal TPR and FPR for nested regression model (no template)
### 
### @param slice the slice of BrS data that are modeled
### @param res_nplcm matrix of MCMC samples
### @param model_options see [nplcm()]
### @param data_nplcm see [nplcm()]
### 
### @return a list; of dimension (number of subjects, dimension of the bronze-standard
### measurement slice, the number of MCMC iterations retained).
### 
###         
### 
##get_marginal_rates_nested_reg <- function(slice, res_nplcm, model_options,data_nplcm){
##  JBrS_curr <- ncol(data_nplcm$Mobs$MBS[[slice]])
##  Jcause    <- length(model_options$likelihood$cause_list)
##  K_curr    <- model_options$likelihood$k_subclass[slice]  
##  N_all     <- length(data_nplcm$Y)
##  res_tpr_all <- array(NA,c(N_all,JBrS_curr,nrow(res_nplcm)))
##  res_fpr_all <- array(NA,c(N_all,JBrS_curr,nrow(res_nplcm)))
##  for (i in 1:N_all){ # begin iteration over subjects:
##    res_tpr <- matrix(0,nrow = nrow(res_nplcm),ncol=JBrS_curr)
##    res_fpr <- matrix(0,nrow = nrow(res_nplcm),ncol=JBrS_curr)
##    is_case <- data_nplcm$Y[i]==1
##    # formula: fittedmean_case[j] = \sum_e pEti[e]*\sum_k {Eta[e,k]*Theta[j,k]*templateBS[e,j]+Eta[e,k]*Psi[j,k]*(1-templateBS[e,j])}
##    for (j in 1:JBrS_curr){
##      # get ThetaBS[j,k]:
##      ind_ThetaBS_tmp <- grep(paste0("^ThetaBS_",slice,"\\[",j,","),colnames(res_nplcm))
##      if (length(ind_ThetaBS_tmp)!=K_curr){stop("==[baker] Check `ThetaBS` extraction from posterior samples! No. of subclasses not matched with specification.==")}
##      ThetaBS_tmp <- res_nplcm[,ind_ThetaBS_tmp]
##      # get PsiBS[j,k]:
##      ind_PsiBS_tmp <- grep(paste0("^PsiBS_",slice,"\\[",j,","),colnames(res_nplcm))
##      if (length(ind_PsiBS_tmp)!=K_curr){stop("==[baker] Check `PsiBS` extraction from posterior samples! No. of subclasses not matched with specification.==")}
##      PsiBS_tmp <- res_nplcm[,ind_PsiBS_tmp]
##      if (is_case){
##        # get Eta[e,k]:
##        ind_Eta_tmp <- grep(paste0("^Eta_",slice,"\\[",i,","),colnames(res_nplcm))
##        if (length(ind_Eta_tmp)!=K_curr){stop("==[baker] Check `Eta` extraction from posterior samples! No. of subclasses not matched with specification.==")}
##        Eta_tmp <- res_nplcm[,ind_Eta_tmp]
##      } else{
##        # get Lambda[k]:
##        ind_Lambda_tmp <- grep(paste0("^Lambda_",slice,"\\[",i,","),colnames(res_nplcm))
##        if (length(ind_Lambda_tmp)!=K_curr){stop("==[baker] Check `Lambda` extraction from posterior samples! No. of subclasses not matched with specification.==")}
##        Lambda_tmp <- res_nplcm[,ind_Lambda_tmp]
##      }
##      # calculate by formula:
##      for (k in 1:K_curr){
##        if(is_case){multiplier <-  Eta_tmp[,k]}else{multiplier <- Lambda_tmp[,k]}
##        res_tpr[,j] <- res_tpr[,j] + multiplier*ThetaBS_tmp[,k]
##        res_fpr[,j] <- res_fpr[,j] + multiplier*PsiBS_tmp[,k]
##      }
##    }
##    res_tpr_all[i,,] <- res_tpr
##    res_fpr_all[i,,] <- res_fpr
##  }# end iteration over subjects.
##  make_list(res_tpr_all,res_fpr_all)
##}
##




#' Plot silver-standard (SS) panel
#' 
#' @param slice the index of measurement slice for SS.
#' @param data_nplcm See [nplcm()]
#' @param model_options See [nplcm()]
#' @param clean_options See [clean_perch_data()]
#' @param res_nplcm See [nplcm_read_folder()]
#' @param bugs.dat Data input for the model fitting.
#' @param bg_color A list with names "BrS", "SS", "pie" to specify background colors
#' @param select_latent a vector of character strings representing latent status. It is used for
#' just plotting a subset of latent status. For example, you can specify `select_latent = "HINF"`
#' to plot all latent status information relevant to `"HINF"`.
#' @param exact Default is `TRUE` to use `select_latent` as exact names of causes. If you want to
#' specify a name and plot all single or combo causes with that name, specify it to be `FALSE`.
#' @param top_SS Numerical value to specify the rightmost limit 
#' on the horizontal axis for the SS panel.
#' @param cexval Default is 1 - size of text of the SS percentages.
#' @param srtval Default is 0 - the direction of the text for the SS percentages.
#' @param prior_shape `interval` or `boxplot` - for how to represent
#' prior/posteriors of the TPR/FPRs of measurements.
#' 
#' @importFrom binom binom.confint
#' @family visualization functions
#' 
#' @return plotting function
#' 
plot_SS_panel <- function(slice,data_nplcm,model_options,
                          clean_options,bugs.dat,res_nplcm,
                          bg_color,
                          select_latent = NULL,
                          exact = TRUE,
                          top_SS = 1, 
                          cexval = 1,
                          srtval = 0,
                          prior_shape="interval"){
  template_SS    <- NULL
  check_combo_SS <- NULL
  if ("SS" %in% model_options$use_measurements){
    template_SS <- lapply(clean_options$SS_objects,"[[","template")
    names(template_SS) <- lapply(clean_options$SS_objects,"[[","nm_spec_test")
    check_combo_SS <- any(unlist(lapply(template_SS,rowSums))>1)
  }
  # order cause_list by posterior means:
  ord <- order_post_eti(res_nplcm,model_options)$ord
  pEti_mat_ord <- order_post_eti(res_nplcm,model_options)$pEti_mat_ord
  
  template_ord <- template_SS[[slice]][ord,,drop=FALSE]
  cause_list <- model_options$likelihood$cause_list
  cause_list_ord <- cause_list[ord]
  
  # focus on selected latent status:
  latent_seq <- get_latent_seq(cause_list,ord,select_latent,exact)$latent_seq
  template_ord <- template_ord[latent_seq,,drop = FALSE]
  
  
  thetaSS_nm <- paste0("thetaSS_",slice)
  alphaS_nm   <- paste0("alphaS_",slice)
  betaS_nm   <- paste0("betaS_",slice)
  # which model was fitted:
  parsed_model <- assign_model(model_options, data_nplcm)
  
  
  if (!any(unlist(parsed_model$regression))){
    
    if (parsed_model$SS_grp){stop("==[baker] Panel plot not available for 
                                  stratified SS TPRs. Please contact maintainer. ==\n")}
    #
    # plcm (just obtain TPR and FPR estimates):
    #
    theta_mat <- res_nplcm[,grep(thetaSS_nm,colnames(res_nplcm)),drop=FALSE]
    theta_mean <- colMeans(theta_mat)
    
    ## model fitted postive rate for each pathogen
    fittedmean_case    <- rep(NA,ncol(template_ord))
    names(fittedmean_case) <- colnames(template_ord)
    
    #fitted_margin_case(pEti_ord,theta,template)
    fitted_margin_case <- function(pEti_ord,theta,template){
      psi = 0
      mixture <-  pEti_ord
      #if (ncol(template)!=length(theta)){
      #  cat(theta,": ",length(theta),",ncol_tmp=",ncol(template),"\n")
      #}
      tpr     <-  t(t(template)*theta)
      fpr     <-  t(t(1-template)*psi)
      colSums(tpr*mixture + fpr*mixture)
    }
    
    fittedmean_case  <- colMeans(t(sapply(1:nrow(pEti_mat_ord),
                                          function(iter)
                                            fitted_margin_case(pEti_mat_ord[iter,latent_seq], 
                                                               theta_mat[iter,],
                                                               template_ord))))
    fittedmean_ctrl <- 0
    
    #plot_pos <- get_plot_pos(template_ord)  # 1 at the 5th means for the fifth smallest etiology, we should plot 1st dimension in this slice.
  }
  
  # get observed rates' summaries:
  Nd <- bugs.dat$Nd
  Nu <- bugs.dat$Nu
  
  MSS_curr <- data_nplcm$Mobs$MSS[[slice]]
  
  # positive rates and confidence intervals:
  #cases:
  MSS_case_curr <- MSS_curr[1:Nd,,drop=FALSE]
  count    <- as.integer(do.call(cbind,lapply(MSS_case_curr,sum,na.rm=TRUE))) #<-- added 'as.integer' to make pathogens appear by rows.
  NA_count <- apply(MSS_case_curr,2,function(v) sum(is.na(v)))
  tmp.case <- binom.confint(count,Nd-NA_count,conf.level = 0.95, methods = "ac")
  
  # case and control positive rate, lower and upper limit
  MSS_mean  <- rbind(round(tmp.case$mean,5))
  MSS_q1 <- rbind(tmp.case[,c("lower")])
  MSS_q2 <- rbind(tmp.case[,c("upper")])
  
  # prior parameters:
  alphaS         <- bugs.dat[[alphaS_nm]]
  betaS          <- bugs.dat[[betaS_nm]]
  
  pos_vec <- get_plot_pos(template_ord)
  
  plot_SS_cell <- function(lat_pos, pos, height,gap = 0){
    plotat <- get_plot_num(lat_pos,height) + gap
    graphics::plot(c(fittedmean_case[pos],MSS_mean[,pos]),
                   plotat[-3],
                   xlim=c(0,top_SS),
                   ylim=c(0.5, height+0.5),
                   xaxt="n",xlab="positive rate",
                   ylab="",yaxt="n",
                   pch = c(2,20),
                   col = c("purple", "dodgerblue2"),
                   cex = c(1,2))
    graphics::points(c(theta_mean[pos],MSS_q1[,pos],MSS_q2[,pos]), # <--- different than BrS here.
                     plotat[c(1,2,2)],
                     pch = c("+","|","|"),
                     col = c("purple",1,1),
                     cex = c(2,1,1))
    # label posterior mean of TPR:
    tmp.post <- as.matrix(theta_mat)[,pos]
    tmp.hpos <- stats::quantile(tmp.post,0.975) + 0.15
    graphics::text(tmp.hpos, lat_pos-0.35+gap, paste0(round(100*theta_mean[pos],1),"%"),
                   srt=srtval,cex=cexval,col="purple")
    
    # case: rates
    graphics::segments(
      x0 = MSS_q1[1,pos],x1 = MSS_q2[1,pos],
      y0 =plotat[2], y1 = plotat[2],
      lty = 1
    )
    tmp.hpos <- ifelse(MSS_q2[1,pos]+0.15>0.95,MSS_q1[1,pos]-0.2,MSS_q2[1,pos]+0.15 )
    graphics::text(tmp.hpos, plotat[2], paste0(round(100*MSS_mean[1,pos],1),"%"),
                   srt=srtval,cex=cexval)
    
    if (!is.null(pos) && !is.na(pos)){#some pos can be NA: because certain cause has no measurements.
      # prior and posterior of TPR:
      if (prior_shape == "interval") {
        # prior of TPR:
        prior_plot_at <- lat_pos - .45 +gap
        tmp = stats::qbeta(c(0.025,0.975,0.25,0.75),alphaS[pos],betaS[pos])
        graphics::points(tmp,rep(prior_plot_at,4),pch = c("|","|","[","]"),col="gray")
        graphics::segments(tmp[1],prior_plot_at,
                           tmp[2],prior_plot_at,lty = 1,col="gray")
        graphics::segments(tmp[3],prior_plot_at,
                           tmp[4],prior_plot_at,lty = 1,col="gray",lwd=2)
        
        # posterior of TPR:
        post_plot_at <- lat_pos - .35 +gap
        tmp.post = as.matrix(theta_mat)[,pos]
        tmp  = stats::quantile(tmp.post, c(0.025,0.975,0.25,0.75))
        graphics::points(tmp,rep(post_plot_at,4),pch = c("|","|","[","]"),col = "purple")
        graphics::segments(tmp[1],post_plot_at,
                           tmp[2],post_plot_at,lty = 1,col = "black")
        graphics::segments(tmp[3],post_plot_at,
                           tmp[4],post_plot_at,lty = 1,col = "black",lwd=2)
      } else if (prior_shape == "boxplot") {
        tmp = stats::rbeta(10000,alphaS[pos],betaS[pos])
        graphics::boxplot(
          tmp,at = prior_plot_at, boxwex = 1 / 10 , col = "gray",
          add = TRUE,horizontal = TRUE,outline = FALSE,xaxt =
            "n"
        )
        tmp.post = as.matrix(theta_mat)[,pos]
        graphics::boxplot(
          tmp.post,at = post_plot_at,boxwex = 1 / 10,add = TRUE,
          horizontal = TRUE,outline = FALSE,xaxt = "n"
        )
      }
      # print name of the dimension (e.g., pathogen):
      graphics::text(top_SS - 0.12,lat_pos + .3+gap, colnames(MSS_case_curr)[pos],cex=1 )
    }
  }
  
  points_SS_cell <- function(lat_pos, pos, height,gap=0){ # pos for the measurement dimension, usually used as pos_vec[e].
    plotat <- get_plot_num(lat_pos,height) +gap
    graphics::points(c(fittedmean_case[pos],MSS_mean[,pos]),
                     plotat[-3],
                     xlim=c(0,top_SS),
                     ylim=c(0.5, height+0.5),
                     xaxt="n",xlab="positive rate",
                     ylab="",yaxt="n",
                     pch = c(2,20),
                     col = c("purple", "dodgerblue2"),
                     cex = c(1,2))
    graphics::points(c(theta_mean[pos],MSS_q1[,pos],MSS_q2[,pos]), # <--- different than BrS here.
                     plotat[c(1,2,2)],
                     pch = c("+","|","|"),
                     col = c("purple",1,1),
                     cex = c(2,1,1))
    # label posterior mean of TPR:
    tmp.post <- as.matrix(theta_mat)[,pos]
    tmp.hpos <- stats::quantile(tmp.post,0.975) + 0.15
    graphics::text(tmp.hpos, lat_pos-0.35+gap, paste0(round(100*theta_mean[pos],1),"%"),
                   srt=srtval,cex=cexval,col="purple")
    
    # case: rates
    graphics::segments(
      x0 = MSS_q1[1,pos],x1 = MSS_q2[1,pos],
      y0 =plotat[2], y1 = plotat[2],
      lty = 1
    )
    tmp.hpos <- ifelse(MSS_q2[1,pos]+0.15>0.95,MSS_q1[1,pos]-0.2,MSS_q2[1,pos]+0.15 )
    graphics::text(tmp.hpos, plotat[2], paste0(round(100*MSS_mean[1,pos],1),"%"),
                   srt=srtval,cex=cexval)
    
    if (!is.null(pos) && !is.na(pos)){#some pos can be NA: because certain cause has no measurements.
      # prior and posterior of TPR:
      if (prior_shape == "interval") {
        # prior of TPR:
        prior_plot_at <- lat_pos - .45 +gap
        tmp = stats::qbeta(c(0.025,0.975,0.25,0.75),alphaS[pos],betaS[pos])
        graphics::points(tmp,rep(prior_plot_at,4),pch = c("|","|","[","]"),col="gray")
        graphics::segments(tmp[1],prior_plot_at,
                           tmp[2],prior_plot_at,lty = 1,col="gray")
        graphics::segments(tmp[3],prior_plot_at,
                           tmp[4],prior_plot_at,lty = 1,col="gray",lwd=2)
        
        # posterior of TPR:
        post_plot_at <- lat_pos - .35 +gap
        tmp.post = as.matrix(theta_mat)[,pos]
        tmp  = stats::quantile(tmp.post, c(0.025,0.975,0.25,0.75))
        graphics::points(tmp,rep(post_plot_at,4),pch = c("|","|","[","]"),col = "purple")
        graphics::segments(tmp[1],post_plot_at,
                           tmp[2],post_plot_at,lty = 1,col = "black")
        graphics::segments(tmp[3],post_plot_at,
                           tmp[4],post_plot_at,lty = 1,col = "black",lwd=2)
      } else if (prior_shape == "boxplot") {
        tmp = stats::rbeta(10000,alphaS[pos],betaS[pos])
        graphics::boxplot(
          tmp,at = prior_plot_at, boxwex = 1 / 10 , col = "gray",
          add = TRUE,horizontal = TRUE,outline = FALSE,xaxt =
            "n"
        )
        tmp.post = as.matrix(theta_mat)[,pos]
        graphics::boxplot(
          tmp.post,at = post_plot_at,boxwex = 1 / 10,add = TRUE,
          horizontal = TRUE,outline = FALSE,xaxt = "n"
        )
      }
      graphics::text(top_SS - 0.12,lat_pos + .3+gap, colnames(MSS_case_curr)[pos],cex=1 )
    }
    
    # x-axis for each cell:
    if (lat_pos>1){
      graphics::axis(1, seq(0,1,by = .1), lwd = 0, lwd.ticks = 0,#labels=rep("",length(seq(0,1,by=.2))),
                     pos = seq(.6,height +.6,by = 1)[lat_pos], cex.axis = 0.8,
                     lty = 2,col = "blue"
      )
    }
  }
  
  Jcause <- length(model_options$likelihood$cause_list)
  #
  # plotting:
  #
  op <- graphics::par(mar=c(5.1,0,4.1,0))
  
  if (check_combo_SS){
    #stop("== Not implemented for combo latent status.==")
    warning("==[baker] Combo latent status implemented with measurements overlapping in SS columns! ==\n")
  }
  
  first  <- TRUE
  cat("== Plotting SS Slice: ", slice, ": ", unlist(names(data_nplcm$Mobs$MSS))[slice])
  for (e in 1:nrow(template_ord)){
    gap_seq <- 0
    if (!is.null(pos_vec[[e]]) && length(pos_vec[[e]])>1){
      gap_seq <- c(0,-0.1*(1:(length(pos_vec[[e]])-1)))
    }
    ct <- 0
    for (pos_curr in pos_vec[[e]]){
      if (!is.na(pos_curr)){
        ct <- ct +1
        if (first) {plot_SS_cell(e,pos_curr,length(latent_seq),gap = gap_seq[ct]);first <- FALSE}
        if (!first) {points_SS_cell(e,pos_curr,length(latent_seq),gap=gap_seq[ct])}
      }
    }
  }
  
  if (!is.null(bg_color) && !is.null(bg_color$SS)){
    graphics::rect(graphics::par("usr")[1], graphics::par("usr")[3], graphics::par("usr")[2], graphics::par("usr")[4], col = 
                     bg_color$SS)
    
    for (e in 1:nrow(template_ord)){
      gap_seq <- 0
      if (!is.null(pos_vec[[e]]) && length(pos_vec[[e]])>1){
        gap_seq <- c(0,-0.1*(1:(length(pos_vec[[e]])-1)))
      }
      
      ct <- 0
      for (pos_curr in pos_vec[[e]]){
        if (!is.na(pos_curr)){
          ct <- ct +1
          points_SS_cell(e,pos_curr,length(latent_seq),gap=gap_seq[ct])
        }
      }
    }
  }
  
  if (sum(template_ord)==0){
    warning(paste0("==[baker] Silver-standard slice ", names(data_nplcm$Mobs$MSS)[slice], 
                   " has no measurements informative of the causes! Please check if measurements' columns correspond to causes.==\n"))  
    plotat <- c(sapply(seq_along(latent_seq),get_plot_num,length(latent_seq)))
    graphics::plot(rep(0,length(plotat)),
                   plotat,
                   xlim=c(0,top_SS),
                   ylim=c(0.5, length(latent_seq)+0.5),
                   xaxt="n",xlab="positive rate",
                   ylab="",yaxt="n",
                   pch = c("",""))
  }
  
  #add ticks from 0 to 1 for x-bar:
  graphics::axis(1,at = c(0,0.2,0.4,0.6,0.8,1),labels= c(0,0.2,0.4,0.6,0.8,1),las=1)
  
  #add dashed lines to separate cells:
  if (length(latent_seq) > 1){
    graphics::abline(h=seq(1.5,length(latent_seq)-.5,by=1),lty=2,lwd=0.5,col="gray")
  }
  
  #add some texts:
  graphics::mtext(eval(paste0("SS: ", names(data_nplcm$Mobs$MSS)[slice])),
                  line=1,cex=1.8)
}

#' Plot etiology (pie) panel
#' 
#' 
#' @param model_options See [nplcm()]
#' @param res_nplcm See [nplcm_read_folder()]
#' @param bugs.dat Data input for the model fitting.
#' @param bg_color A list with names "BrS", "SS", "pie" to specify background colors
#' @param select_latent a vector of character strings representing latent status. It is used for
#' just plotting a subset of latent status. For example, you can specify `select_latent = "HINF"`
#' @param exact Default is `TRUE` to use `select_latent` as exact names of causes. If you want to
#' specify a name and plot all single or combo causes with that name, specify it to be `FALSE`.
#' to plot all latent status information relevant to `"HINF"`.
#' @param top_pie Numerical value to specify the rightmost limit 
#' on the horizontal axis for the pie panel.
#' @param label_size the size of latent status labels on the right margin
#' @param ref_eti reference quantiles and means; a list: pEti_ref_q, pEti_ref_mean_ord
#' @param is_plot default to `TRUE` for plotting only; set to `FALSE` if to get summary.
#' 
#' @importFrom binom binom.confint
#' @family visualization functions
#' 
#' @return plotting function.
plot_pie_panel <- function(model_options,
                           res_nplcm,
                           bugs.dat,
                           bg_color,
                           select_latent = NULL,
                           exact = TRUE,
                           top_pie = 1,
                           label_size = 1,
                           ref_eti = NULL,is_plot=TRUE){
  
  res <- list()
  # order cause_list by posterior means:
  ord <- order_post_eti(res_nplcm,model_options)$ord
  pEti_mat_ord <- order_post_eti(res_nplcm,model_options)$pEti_mat_ord
  pEti_mean_ord <- colMeans(pEti_mat_ord)
  pEti_sd_ord <- apply(pEti_mat_ord,2,stats::sd)
  
  # quantiles for etiology: outer is 97.5% CI, inner is 50% CI
  pEti_q   <- apply(pEti_mat_ord,2,stats::quantile,probs=c(0.025,0.975,0.25,0.75))
  
  cause_list <- model_options$likelihood$cause_list
  cause_list_ord <- cause_list[ord]
  
  Nd <- bugs.dat$Nd
  Nu <- bugs.dat$Nu
  
  # focus on selected latent status:
  latent_seq <- get_latent_seq(cause_list,ord,select_latent,exact)$latent_seq
  
  original_num <- get_latent_seq(cause_list,ord,select_latent,exact)$original_num
  pEti_mat_ord   <- pEti_mat_ord[,latent_seq,drop=FALSE]
  pEti_mean_ord  <- pEti_mean_ord[latent_seq]
  pEti_sd_ord  <- pEti_sd_ord[latent_seq]
  pEti_q         <- pEti_q[,latent_seq,drop=FALSE]
  cause_list_ord <- cause_list_ord[latent_seq]
  
  
  Jcause <- length(model_options$likelihood$cause_list)
  alpha_ord <- bugs.dat$alphaEti[ord]
  
  if(is_plot){
    plot_pie_cell_first <- function(lat_pos,height,dotcolor="black",add=FALSE){
      # posterior mean of etiology:
      if (!add){
        graphics::plot(pEti_mean_ord[lat_pos],lat_pos,
                       yaxt="n",
                       xlim=c(0,top_pie),ylim=c(0.5,height+0.5),
                       col="purple",
                       ylab="",xlab="probability",
                       pch= 20,cex=2)
      }
      if (add){
        graphics::points(pEti_mean_ord[lat_pos],lat_pos,
                         yaxt="n",
                         xlim=c(0,top_pie),ylim=c(0.5,height+0.5),
                         col="purple",
                         ylab="",xlab="probability",
                         pch= 20,cex=2)
        
      }
    }
    
    points_pie_cell <- function(lat_pos,height,dotcolor="black"){
      if (lat_pos > 1){
        # posterior mean of etiology:
        graphics::points(pEti_mean_ord[lat_pos],lat_pos,
                         yaxt="n",
                         xlim=c(0,top_pie),ylim=c(lat_pos-0.5,lat_pos+0.5),
                         col="purple",
                         ylab="",xlab="probability",
                         pch= 20,cex=2)
      }
      # x-axis for each cell:
      if (lat_pos>1){
        graphics::axis(1, seq(0,1,by = .2), lwd = 0, lwd.ticks = 0,#labels=rep("",length(seq(0,1,by=.2))),
                       pos = seq(.625,height +.625,by = 1)[lat_pos], cex.axis = 0.8,lty =
                         2,col = "blue"
        )
      }
      
      graphics::points(pEti_q[,lat_pos],rep(lat_pos,4),pch=c("|","|","[","]"),cex=1,col="purple")
      
      pgrid = seq(0,1,by=0.01)
      graphics::segments(y0=lat_pos,x0=pEti_q[1,lat_pos],y1=lat_pos,x1=pEti_q[2,lat_pos],col=dotcolor)
      graphics::segments(y0=lat_pos,x0=pEti_q[3,lat_pos],y1=lat_pos,x1=pEti_q[4,lat_pos],col=dotcolor, lwd=2)
      graphics::text(.8,lat_pos,paste0("=",paste0(round(100*c(pEti_mean_ord),1)[lat_pos],"%")),srt=0,cex=2)
      if (!is.null(ref_eti)){
        match_id_ref <- match(cause_list_ord,colnames(ref_eti$pEti_ref_q))
        graphics::segments(y0=lat_pos+0.25,x0=ref_eti$pEti_ref_q[1,match_id_ref[lat_pos]],
                           y1=lat_pos+0.25,x1=ref_eti$pEti_ref_q[2,match_id_ref[lat_pos]],col=dotcolor,lwd=1)
        graphics::points(ref_eti$pEti_ref_mean_ord[match_id_ref[lat_pos]],lat_pos+0.25,col="orange",cex=2)
        graphics::text(.8,lat_pos+0.25,paste0("=",paste0(round(100*c(ref_eti$pEti_ref_mean_ord),1)[match_id_ref[lat_pos]],"%")),srt=0,cex=2,col="orange")
      }
      graphics::text(.65,lat_pos,bquote(hat(pi)[.(ord[lat_pos])]),srt=0,cex=2)
      #prior density:
      tmp.density = stats::dbeta(pgrid,alpha_ord[lat_pos],sum(alpha_ord[-lat_pos]))
      graphics::points(pgrid,tmp.density/(3*max(tmp.density))+lat_pos-0.45,type="l",col="gray",lwd=4,lty=2)
      #posterior density:
      tmp.post.density = stats::density(pEti_mat_ord[,lat_pos],from=0,to=1)
      tmp.x = tmp.post.density$x
      tmp.y = tmp.post.density$y
      graphics::points(tmp.x,tmp.y/(3*max(tmp.y))+lat_pos-0.45,lwd=4,type="l",col="purple")
      
    }
    
    #
    # plot etiology information:
    #
    
    op <- graphics::par(mar=c(5.1,0,4.1,10))
    
    cat("== Plotting pies == \n")
    
    first <- TRUE
    for (e in seq_along(latent_seq)){
      if (first){plot_pie_cell_first(e,length(latent_seq))}
      points_pie_cell(e,length(latent_seq))
      first <- FALSE
    }
    
    if (!is.null(bg_color) && !is.null(bg_color$pie)){
      graphics::rect(graphics::par("usr")[1], graphics::par("usr")[3], graphics::par("usr")[2], graphics::par("usr")[4], col = 
                       bg_color$pie)
      
      first <- TRUE
      for (e in seq_along(latent_seq)){
        if (first){plot_pie_cell_first(e,length(latent_seq),add=TRUE)}
        points_pie_cell(e,length(latent_seq))
        first <- FALSE
      }
    }
    # cause names on the right edge:
    graphics::axis(4,at=1:length(latent_seq),labels=paste(paste(cause_list_ord,original_num,sep=" ("),")",sep=""),
                   las=2,cex.axis=label_size)
    # cell bottom axis:
    if (length(latent_seq)>1){
      graphics::abline(h=seq(1.5,length(latent_seq)-.5,by=1),lty=2,lwd=0.5,col="gray")
    }
    
    graphics::mtext(expression(underline(hat(pi))),line=1,cex=1.8)
    graphics::legend("topright",c("prior","posterior"),lty=c(2,1),col=c("gray","purple"),
                     lwd = 4,horiz=TRUE,cex=1,bty="n")
    
  }
  ord_now <- ord[latent_seq] 
  # return to original non-ordered results:
  res0 <- make_list(cause_list_ord[order(ord_now)],pEti_mean_ord[order(ord_now)],
                    pEti_sd_ord[order(ord_now)],pEti_q[c(1,2),order(ord_now),drop=FALSE])
  res <- t(do.call("rbind",res0[-1]))
  rownames(res) <- res0[[1]]
  colnames(res) <- c("post.mean","post.sd","CrI_025","CrI_0975")
  
  if (!is_plot){return(res)}
}
