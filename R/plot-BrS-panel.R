#' Plot bronze-standard (BrS) panel
#' 
#' 
#' @param slice the index of measurement slice for BrS.
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
#' @param top_BrS Numerical value to specify the rightmost limit 
#' on the horizontal axis for the BrS panel.
#' @param cexval Default is 1 - size of text of the BrS percentages.
#' @param srtval Default is 0 - the direction of the text for the BrS percentages.
#' @param prior_shape \code{interval} or \code{boxplot} - for how to represent
#' prior/posteriors of the TPR/FPRs of measurements.
#' @param silent Default is TRUE to not print any warning messages; FALSE otherwise.
#' 
#' @family visualization functions
#' @export

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
  cat("\n == Plotting BrS Slice: ", slice, ": ", unlist(names(data_nplcm$Mobs$MBS))[slice])
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
#' @examples 
#' \dontrun{
#' result_folder <- c("C:/2015_09_17_01KEN")
#' out           <- nplcm_read_folder(result_folder)
#' data_nplcm    <- list(Mobs  = out$Mobs, Y = out$Y)
#' slice         <- 1
#' # fitted positive rates for pathogens separately among cases and controls:
#' get_fitted_mean_no_nested(slice,out$res_nplcm,out$model_options,data_nplcm,
#'                       out$clean_options)
#' # names of pathogens:
#' colnames(out$Mobs$MBS[[slice]])
#' }
#' @return a list with model fitted means
#' @export
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
#' @param model_options see \code{\link{nplcm}}
#' @param data_nplcm see \code{\link{nplcm}}
#' @param clean_options see \code{\link{clean_perch_data}}
#' @return a matrix of no. of rows equal to retained MCMC samples, no. of columns
#' equal to the no. of measurement dimensions within a slice.
#' 
#' @examples 
#' \dontrun{
#' result_folder <- c("C:/2015_09_17_01KEN_nplcm")
#' out           <- nplcm_read_folder(result_folder)
#' data_nplcm    <- list(Mobs  = out$Mobs, Y = out$Y)
#' slice         <- 1
#' # fitted positive rates for pathogens separately among cases and controls:
#' get_fitted_mean_nested(slice,out$res_nplcm,out$model_options,data_nplcm,
#'                       out$clean_options)
#' # names of pathogens:
#' colnames(out$Mobs$MBS[[slice]])
#' }
#' 
#' @export
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
#' @param model_options see \code{\link{nplcm}}
#' @param data_nplcm see \code{\link{nplcm}}
#' 
#' @return a matrix of no. of rows equal to retained MCMC samples, no. of columns
#' equal to the no. of measurement dimensions within a slice.
#' 
#' @export
#' 
#' 
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
#' @param model_options see \code{\link{nplcm}}
#' @param data_nplcm see \code{\link{nplcm}}
#' 
#' @return a matrix of no. of rows equal to retained MCMC samples, no. of columns
#' equal to the no. of measurement dimensions within a slice.
#' 
#' @export
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



###' get marginal TPR and FPR for nested regression model (no template)
###' 
###' @param slice the slice of BrS data that are modeled
###' @param res_nplcm matrix of MCMC samples
###' @param model_options see \code{\link{nplcm}}
###' @param data_nplcm see \code{\link{nplcm}}
###' 
###' @return a list; of dimension (number of subjects, dimension of the bronze-standard
###' measurement slice, the number of MCMC iterations retained).
###' 
###' @export
###' 
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

