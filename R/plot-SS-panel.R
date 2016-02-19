#' Plot silver-standard (SS) panel
#' 
#' 
#' @param slice the index of measurement slice for SS.
#' @param data_nplcm See \code{\link{nplcm}}
#' @param model_options See \code{\link{nplcm}}
#' @param clean_options See \code{\link{clean_perch_data}}
#' @param res_nplcm See \code{\link{nplcm_read_folder}}
#' @param bugs.dat Data input for the model fitting.
#' @param bg_color A list with names "BrS", "SS", "pie" to specify background colors
#' @param select_latent a vector of character strings representing latent status. It is used for
#' just plotting a subset of latent status. For example, you can specify \code{select_latent = "HINF"}
#' to plot all latent status information relevant to \code{"HINF"}.
#' @param exact Default is \code{TRUE} to use \code{select_latent} as exact names of causes. If you want to
#' specify a name and plot all single or combo causes with that name, specify it to be \code{FALSE}.
#' @param top_SS Numerical value to specify the rightmost limit 
#' on the horizontal axis for the SS panel.
#' @param cexval Default is 1 - size of text of the SS percentages.
#' @param srtval Default is 0 - the direction of the text for the SS percentages.
#' @param prior_shape \code{interval} or \code{boxplot} - for how to represent
#' prior/posteriors of the TPR/FPRs of measurements.
#' 
#' @importFrom binom binom.confint
#' @family visualization functions
#' @export


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
    
    if (parsed_model$SS_grp){stop("== Panel plot not available for stratified SS TPRs. Please contact maintainer. ==")}
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
    plot(c(fittedmean_case[pos],MSS_mean[,pos]),
         plotat[-3],
         xlim=c(0,top_SS),
         ylim=c(0.5, height+0.5),
         xaxt="n",xlab="positive rate",
         ylab="",yaxt="n",
         pch = c(2,20),
         col = c("purple", "dodgerblue2"),
         cex = c(1,2))
    points(c(theta_mean[pos],MSS_q1[,pos],MSS_q2[,pos]), # <--- different than BrS here.
           plotat[c(1,2,2)],
           pch = c("+","|","|"),
           col = c("purple",1,1),
           cex = c(2,1,1))
    # label posterior mean of TPR:
    tmp.post <- as.matrix(theta_mat)[,pos]
    tmp.hpos <- quantile(tmp.post,0.975) + 0.15
    text(tmp.hpos, lat_pos-0.35+gap, paste0(round(100*theta_mean[pos],1),"%"),
         srt=srtval,cex=cexval,col="purple")
    
    # case: rates
    segments(
      x0 = MSS_q1[1,pos],x1 = MSS_q2[1,pos],
      y0 =plotat[2], y1 = plotat[2],
      lty = 1
    )
    tmp.hpos <- ifelse(MSS_q2[1,pos]+0.15>0.95,MSS_q1[1,pos]-0.2,MSS_q2[1,pos]+0.15 )
    text(tmp.hpos, plotat[2], paste0(round(100*MSS_mean[1,pos],1),"%"),
         srt=srtval,cex=cexval)
    
    if (!is.null(pos) && !is.na(pos)){#some pos can be NA: because certain cause has no measurements.
      # prior and posterior of TPR:
      if (prior_shape == "interval") {
        # prior of TPR:
        prior_plot_at <- lat_pos - .45 +gap
        tmp = qbeta(c(0.025,0.975,0.25,0.75),alphaS[pos],betaS[pos])
        points(tmp,rep(prior_plot_at,4),pch = c("|","|","[","]"),col="gray")
        segments(tmp[1],prior_plot_at,
                 tmp[2],prior_plot_at,lty = 1,col="gray")
        segments(tmp[3],prior_plot_at,
                 tmp[4],prior_plot_at,lty = 1,col="gray",lwd=2)
        
        # posterior of TPR:
        post_plot_at <- lat_pos - .35 +gap
        tmp.post = as.matrix(theta_mat)[,pos]
        tmp  = quantile(tmp.post, c(0.025,0.975,0.25,0.75))
        points(tmp,rep(post_plot_at,4),pch = c("|","|","[","]"),col = "purple")
        segments(tmp[1],post_plot_at,
                 tmp[2],post_plot_at,lty = 1,col = "black")
        segments(tmp[3],post_plot_at,
                 tmp[4],post_plot_at,lty = 1,col = "black",lwd=2)
      } else if (prior_shape == "boxplot") {
        tmp = rbeta(10000,alphaS[pos],betaS[pos])
        boxplot(
          tmp,at = prior_plot_at, boxwex = 1 / 10 , col = "gray",
          add = TRUE,horizontal = TRUE,outline = FALSE,xaxt =
            "n"
        )
        tmp.post = as.matrix(theta_mat)[,pos]
        boxplot(
          tmp.post,at = post_plot_at,boxwex = 1 / 10,add = TRUE,
          horizontal = TRUE,outline = FALSE,xaxt = "n"
        )
      }
      # print name of the dimension (e.g., pathogen):
      text(top_SS - 0.12,lat_pos + .3+gap, colnames(MSS_case_curr)[pos],cex=1 )
    }
  }
  
  points_SS_cell <- function(lat_pos, pos, height,gap=0){ # pos for the measurement dimension, usually used as pos_vec[e].
    plotat <- get_plot_num(lat_pos,height) +gap
    points(c(fittedmean_case[pos],MSS_mean[,pos]),
           plotat[-3],
           xlim=c(0,top_SS),
           ylim=c(0.5, height+0.5),
           xaxt="n",xlab="positive rate",
           ylab="",yaxt="n",
           pch = c(2,20),
           col = c("purple", "dodgerblue2"),
           cex = c(1,2))
    points(c(theta_mean[pos],MSS_q1[,pos],MSS_q2[,pos]), # <--- different than BrS here.
           plotat[c(1,2,2)],
           pch = c("+","|","|"),
           col = c("purple",1,1),
           cex = c(2,1,1))
    # label posterior mean of TPR:
    tmp.post <- as.matrix(theta_mat)[,pos]
    tmp.hpos <- quantile(tmp.post,0.975) + 0.15
    text(tmp.hpos, lat_pos-0.35+gap, paste0(round(100*theta_mean[pos],1),"%"),
         srt=srtval,cex=cexval,col="purple")
    
    # case: rates
    segments(
      x0 = MSS_q1[1,pos],x1 = MSS_q2[1,pos],
      y0 =plotat[2], y1 = plotat[2],
      lty = 1
    )
    tmp.hpos <- ifelse(MSS_q2[1,pos]+0.15>0.95,MSS_q1[1,pos]-0.2,MSS_q2[1,pos]+0.15 )
    text(tmp.hpos, plotat[2], paste0(round(100*MSS_mean[1,pos],1),"%"),
         srt=srtval,cex=cexval)
    
    if (!is.null(pos) && !is.na(pos)){#some pos can be NA: because certain cause has no measurements.
      # prior and posterior of TPR:
      if (prior_shape == "interval") {
        # prior of TPR:
        prior_plot_at <- lat_pos - .45 +gap
        tmp = qbeta(c(0.025,0.975,0.25,0.75),alphaS[pos],betaS[pos])
        points(tmp,rep(prior_plot_at,4),pch = c("|","|","[","]"),col="gray")
        segments(tmp[1],prior_plot_at,
                 tmp[2],prior_plot_at,lty = 1,col="gray")
        segments(tmp[3],prior_plot_at,
                 tmp[4],prior_plot_at,lty = 1,col="gray",lwd=2)
        
        # posterior of TPR:
        post_plot_at <- lat_pos - .35 +gap
        tmp.post = as.matrix(theta_mat)[,pos]
        tmp  = quantile(tmp.post, c(0.025,0.975,0.25,0.75))
        points(tmp,rep(post_plot_at,4),pch = c("|","|","[","]"),col = "purple")
        segments(tmp[1],post_plot_at,
                 tmp[2],post_plot_at,lty = 1,col = "black")
        segments(tmp[3],post_plot_at,
                 tmp[4],post_plot_at,lty = 1,col = "black",lwd=2)
      } else if (prior_shape == "boxplot") {
        tmp = rbeta(10000,alphaS[pos],betaS[pos])
        boxplot(
          tmp,at = prior_plot_at, boxwex = 1 / 10 , col = "gray",
          add = TRUE,horizontal = TRUE,outline = FALSE,xaxt =
            "n"
        )
        tmp.post = as.matrix(theta_mat)[,pos]
        boxplot(
          tmp.post,at = post_plot_at,boxwex = 1 / 10,add = TRUE,
          horizontal = TRUE,outline = FALSE,xaxt = "n"
        )
      }
      text(top_SS - 0.12,lat_pos + .3+gap, colnames(MSS_case_curr)[pos],cex=1 )
    }
  }
  
  Jcause <- length(model_options$likelihood$cause_list)
  #
  # plotting:
  #
  op <- par(mar=c(5.1,0,4.1,0))
  
  if (!is_length_all_one(pos_vec)){
    #stop("== Not implemented for combo latent status.==")
    warning("== Combo latent status implemented with measurements overlapping in SS columns! ==")
  }
  
  first  <- TRUE
  cat("\n == Plotting SS Slice: ", slice, ": ", unlist(names(data_nplcm$Mobs$MSS))[slice])
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
    rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = 
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
    warning(paste0("== Silver-standard slice ", names(data_nplcm$Mobs$MSS)[slice], " has no measurements informative of the causes! Please check if measurements' columns correspond to causes.=="))  
    plotat <- c(sapply(seq_along(latent_seq),get_plot_num,length(latent_seq)))
    plot(rep(0,length(plotat)),
           plotat,
           xlim=c(0,top_SS),
           ylim=c(0.5, length(latent_seq)+0.5),
           xaxt="n",xlab="positive rate",
           ylab="",yaxt="n",
           pch = c("",""))
  }
  
  #add ticks from 0 to 1 for x-bar:
  axis(1,at = c(0,0.2,0.4,0.6,0.8,1),labels= c(0,0.2,0.4,0.6,0.8,1),las=1)
  
  #add dashed lines to separate cells:
  if (length(latent_seq) > 1){
    abline(h=seq(1.5,length(latent_seq)-.5,by=1),lty=2,lwd=0.5,col="gray")
  }
  
  #add some texts:
  mtext(eval(paste0("SS: ", names(data_nplcm$Mobs$MSS)[slice])),
        line=1,cex=1.8)
}

