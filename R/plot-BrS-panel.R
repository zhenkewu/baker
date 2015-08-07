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
#' @param top_BrS Numerical value to specify the rightmost limit 
#' on the horizontal axis for the BrS panel.
#' @param cexval Default is 1 - size of text of the BrS percentages.
#' @param srtval Default is 0 - the direction of the text for the BrS percentages.
#' @param prior_shape \code{interval} or \code{boxplot} - for how to represent
#' prior/posteriors of the TPR/FPRs of measurements.
#' @param silent Default is TRUE to not print any warning messages; FALSE otherwise.
#' @importFrom binom binom.confint
#' 
#' @export

plot_BrS_panel <- function(slice,data_nplcm,model_options,
                           clean_options,bugs.dat,res_nplcm,
                           bg_color,
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
  
  thetaBS_nm <- paste0("thetaBS_",slice)
  psiBS_nm   <- paste0("psiBS_",slice)
  alphaB_nm   <- paste0("alphaB_",slice)
  betaB_nm   <- paste0("betaB_",slice)
  # which model was fitted:
  parsed_model <- assign_model(model_options, data_nplcm)
  
  if (!parsed_model$nested & !any(unlist(parsed_model$regression))){
    #
    # plcm (just obtain TPR and FPR estimates):
    #
    theta_mat <- res_nplcm[,grep(thetaBS_nm,colnames(res_nplcm)),drop=FALSE]
    theta_mean <- colMeans(theta_mat)
    
    #posterior distribution of FPR:
    psi_mat  <- res_nplcm[,grep(psiBS_nm,colnames(res_nplcm)),drop=FALSE]
    psi_mean   <- colMeans(psi_mat)
    
    ## model fitted postive rate for each pathogen
    fittedmean_case    <- rep(NA,ncol(template_ord))
    names(fittedmean_case) <- colnames(template_ord)
    
    fitted_margin_case <- function(pEti_ord,theta,psi,template){
      mixture <-  pEti_ord
      tpr     <-  t(t(template)*theta)
      fpr     <- t(t(1-template)*psi)
      colSums(tpr*mixture + fpr*mixture)
    }
    
    fittedmean_case  <- colMeans(t(sapply(1:nrow(pEti_mat_ord),
                                          function(iter)
                                            fitted_margin_case(pEti_mat_ord[iter,], theta_mat[iter,],
                                                               psi_mat[iter,],template_ord))))
    fittedmean_control <- psi_mean
    
    #plot_pos <- get_plot_pos(template_ord)  # 1 at the 5th means for the fifth smallest etiology, we should plot 1st dimension in this slice.
  }
  
  # get observed rates' summaries:
  Nd <- bugs.dat$Nd
  Nu <- bugs.dat$Nu
  
  MBS_curr <- data_nplcm$Mobs$MBS[[slice]]
  
  # positive rates and confidence intervals:
  #cases:
  MBS_case_curr <- MBS_curr[1:Nd,,drop=FALSE]
  count    <- do.call(cbind,lapply(MBS_case_curr,table))["1",]
  NA_count <- apply(MBS_case_curr,2,function(v) sum(is.na(v)))
  tmp.case <- binom.confint(count,Nd-NA_count,conf.level = 0.95, methods = "ac")
  
  #controls:
  MBS_ctrl_curr <- MBS_curr[-(1:Nd),,drop=FALSE]
  count    <- do.call(cbind,lapply(MBS_ctrl_curr,table))["1",]
  NA_count <- apply(MBS_ctrl_curr,2,function(v) sum(is.na(v)))
  tmp.ctrl <- binom.confint(count, Nu-NA_count, conf.level = 0.95, methods = "ac")
  
  # case and control positive rate, lower and upper limit
  MBS_mean  <- rbind(round(tmp.case$mean,5),round(tmp.ctrl$mean,5))
  MBS_q1 <- rbind(tmp.case[,c("lower")],tmp.ctrl[,c("lower")])
  MBS_q2 <- rbind(tmp.case[,c("upper")],tmp.ctrl[,c("upper")])
  
  # prior parameters:
  alphaB         <- bugs.dat[[alphaB_nm]]
  betaB          <- bugs.dat[[betaB_nm]]

  pos_vec <- get_plot_pos(template_ord)
  
  get_COR <- function(pos){
    if (is.null(pos) || is.na(pos)){return(NULL)}
    y <- c(rep(1,Nd), rep(0,Nu))
    brs.data <- as.data.frame(rbind(MBS_case_curr,MBS_ctrl_curr))
    dat.reg  <- as.data.frame(cbind(y,brs.data))
    fit      <- glm(y~.,data=dat.reg,family=binomial,na.action="na.omit")
    
    if (sum(is.na(coef(fit)))==0 & sum(diag(vcov(fit))>100)==0){
      res0 = cbind(exp(suppressMessages(confint(fit))),exp(fit$coef))[-1,]
      res  = list(ORinterval = matrix(res0,nrow = ncol(brs.data),ncol=3)[pos,],
                  label = "conditional OR")
    }else{
      if (!silent){
       print("Conditional OR not calculatable. Switch to mariginal OR.")
      }
      res0 = matrix(NA,nrow=1,ncol=3)
      res0 = data.frame(res0)
      tb <- table(dat.reg$y,brs.data[,pos])
      fit_tmp  <- glm(y~brs.data[,pos],family=binomial,na.action = "na.omit")
      if (length(vcov(fit_tmp))>1 && vcov(fit_tmp)[2,2]<100 && ncol(tb)==2){
        #print(l)
        res0 <- cbind(exp(suppressMessages(confint(fit_tmp))),exp(fit_tmp$coef))[-1,]
      }
      res = list(ORinterval=res0,label="marginal OR")
    }
    res
  }
  
  plot_BrS_cell  <- function(lat_pos, pos, height,gap = 0){
    plotat <- get_plot_num(lat_pos,height) + gap
    
    plot(c(fittedmean_case[pos],MBS_mean[,pos]),
           plotat,
           xlim=c(0,top_BrS),
           ylim=c(0.5, height+0.5),
           xaxt="n",xlab="positive rate",
           ylab="",yaxt="n",
           pch = c(2,20,20),
           col = c("purple","dodgerblue2", "dodgerblue2"),
           cex = c(1,2,2))
    
    
    points(c(theta_mean[pos],MBS_q2[,pos]),
           plotat,
           pch = c("+","|","|"),
           col = c("purple",1,1),
           cex = c(2,1,1))
    points(c(psi_mean[pos],MBS_q1[,pos]),
           plotat,
           pch = c("*","|","|"),
           col = c("purple",1,1),
           cex = c(2,1,1))
    # connect case and control rates:
    segments(
      x0 = MBS_mean[1,pos],x1 = MBS_mean[2,pos],
      y0 = plotat[2],y1 = plotat[3],
      lty = 1,
      col = "dodgerblue2",
      lwd = 2
    )
    # case: rates
    segments(
      x0 = MBS_q1[1,pos],x1 = MBS_q2[1,pos],
      y0 = plotat[2], y1 = plotat[2],
      lty = 1
    )
    tmp.hpos <- ifelse(MBS_q2[1,pos]+0.15>0.95,MBS_q1[1,pos]-0.2,MBS_q2[1,pos]+0.15 )
    text(tmp.hpos, plotat[2], paste0(round(100*MBS_mean[1,pos],1),"%"),
         srt=srtval,cex=cexval)
    # control:rates
    segments(
      x0 = MBS_q1[2,pos],x1 = MBS_q2[2,pos],
      y0 = plotat[3], y1 = plotat[3],
      lty = 1
    )
    tmp.hpos <- ifelse(MBS_q2[2,pos]+0.15>0.95,MBS_q1[2,pos]-0.2,MBS_q2[2,pos]+0.15 )
    text(tmp.hpos, plotat[3], paste0(round(100*MBS_mean[2,pos],1),"%"),
         srt=srtval,cex=cexval)
    # poster means of TPR, FPR and fitted marginal rate:
    segments(
      x0 = theta_mean[pos],x1 = psi_mean[pos],
      y0 = plotat[1], y1 = plotat[1],
      lty = 4,col="gray"
    )
    
    if (!is.null(pos) && !is.na(pos)){#some pos can be NA: because certain cause has no measurements.
      # prior and posterior of TPR:
      if (prior_shape == "interval") {
        # prior of TPR:
        prior_plot_at <- lat_pos - .45 + gap
        tmp = qbeta(c(0.025,0.975,0.25,0.75),alphaB[pos],betaB[pos])
        points(tmp,rep(prior_plot_at,4),pch = c("|","|","[","]"),col="gray")
        segments(tmp[1],prior_plot_at,
                 tmp[2],prior_plot_at,lty = 1,col="gray")
        segments(tmp[3],prior_plot_at,
                 tmp[4],prior_plot_at,lty = 1,col="gray",lwd=2)
        
        # posterior of TPR:
        post_plot_at <- lat_pos - .35 + gap
        tmp.post = as.matrix(theta_mat)[,pos]
        tmp  = quantile(tmp.post, c(0.025,0.975,0.25,0.75))
        points(tmp,rep(post_plot_at,4),pch = c("|","|","[","]"),col = "purple")
        segments(tmp[1],post_plot_at,
                 tmp[2],post_plot_at,lty = 1)
        segments(tmp[3],post_plot_at,
                 tmp[4],post_plot_at,lty = 1,lwd=2)
      } else if (prior_shape == "boxplot") {
        tmp = rbeta(10000,alphaB[pos],betaB[pos])
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
      
      #plot conditional odds ratio on the right:
      tmp0 <- get_COR(pos)
      tmp  <- tmp0$ORinterval
      
      L <-  round(tmp[1],1)
      C <-  round(tmp[3],1)
      R <-  round(tmp[2],1)
      
      text(top_BrS - 0.12,lat_pos + .3, colnames(MBS_case_curr)[pos],cex=1 )
      text(top_BrS - 0.12,lat_pos + 1 / (2 * Jcause),C,cex = 1.5)
      text(top_BrS - 0.12,lat_pos - .2,paste(c(L,"   ",R),collapse = " "),
           cex = 1.2)
      legend("topright",tmp0$label,bty = "n")
    }
    
    # x-axis for each cell:
    if (lat_pos>1){
      axis(1, seq(0,1,by = .05), lwd = 0, lwd.ticks = 0,#labels=rep("",length(seq(0,1,by=.2))),
           pos = seq(.625,Jcause +.625,by = 1)[lat_pos], cex.axis = 0.8,
           lty = 2,col = "blue"
      )
    }
  }
  
  points_BrS_cell     <- function(lat_pos,pos,height,gap=0){ # pos for the measurement dimension, usually used as pos_vec[e].
    plotat <- get_plot_num(lat_pos,height) + gap
    
    points(c(fittedmean_case[pos],MBS_mean[,pos]),
           plotat,
           xlim=c(0,top_BrS),
           ylim=c(0.5, height+0.5),
           xaxt="n",xlab="positive rate",
           ylab="",yaxt="n",
           pch = c(2,20,20),
           col = c("purple","dodgerblue2", "dodgerblue2"),
           cex = c(1,2,2))
    
    
    points(c(theta_mean[pos],MBS_q2[,pos]),
           plotat,
           pch = c("+","|","|"),
           col = c("purple",1,1),
           cex = c(2,1,1))
    points(c(psi_mean[pos],MBS_q1[,pos]),
           plotat,
           pch = c("*","|","|"),
           col = c("purple",1,1),
           cex = c(2,1,1))
    # connect case and control rates:
    segments(
      x0 = MBS_mean[1,pos],x1 = MBS_mean[2,pos],
      y0 = plotat[2],y1 = plotat[3],
      lty = 1,
      col = "dodgerblue2",
      lwd = 2
    )
    # case: rates
    segments(
      x0 = MBS_q1[1,pos],x1 = MBS_q2[1,pos],
      y0 = plotat[2], y1 = plotat[2],
      lty = 1
    )
    tmp.hpos <- ifelse(MBS_q2[1,pos]+0.15>0.95,MBS_q1[1,pos]-0.2,MBS_q2[1,pos]+0.15 )
    text(tmp.hpos, plotat[2], paste0(round(100*MBS_mean[1,pos],1),"%"),
         srt=srtval,cex=cexval)
    # control:rates
    segments(
      x0 = MBS_q1[2,pos],x1 = MBS_q2[2,pos],
      y0 = plotat[3], y1 = plotat[3],
      lty = 1
    )
    tmp.hpos <- ifelse(MBS_q2[2,pos]+0.15>0.95,MBS_q1[2,pos]-0.2,MBS_q2[2,pos]+0.15 )
    text(tmp.hpos, plotat[3], paste0(round(100*MBS_mean[2,pos],1),"%"),
         srt=srtval,cex=cexval)
    # poster means of TPR, FPR and fitted marginal rate:
    segments(
      x0 = theta_mean[pos],x1 = psi_mean[pos],
      y0 = plotat[1], y1 = plotat[1],
      lty = 4,col="gray"
    )
    
    if (!is.null(pos) && !is.na(pos)){#some pos can be NA: because certain cause has no measurements.
      # prior and posterior of TPR:
      if (prior_shape == "interval") {
        # prior of TPR:
        prior_plot_at <- lat_pos - .45 + gap
        tmp = qbeta(c(0.025,0.975,0.25,0.75),alphaB[pos],betaB[pos])
        points(tmp,rep(prior_plot_at,4),pch = c("|","|","[","]"),col="gray")
        segments(tmp[1],prior_plot_at,
                 tmp[2],prior_plot_at,lty = 1,col="gray")
        segments(tmp[3],prior_plot_at,
                 tmp[4],prior_plot_at,lty = 1,col="gray",lwd=2)
        
        # posterior of TPR:
        post_plot_at <- lat_pos - .35 + gap
        tmp.post = as.matrix(theta_mat)[,pos]
        tmp  = quantile(tmp.post, c(0.025,0.975,0.25,0.75))
        points(tmp,rep(post_plot_at,4),pch = c("|","|","[","]"),col = "purple")
        segments(tmp[1],post_plot_at,
                 tmp[2],post_plot_at,lty = 1)
        segments(tmp[3],post_plot_at,
                 tmp[4],post_plot_at,lty = 1,lwd=2)
      } else if (prior_shape == "boxplot") {
        tmp = rbeta(10000,alphaB[pos],betaB[pos])
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
      
      #plot conditional odds ratio on the right:
      tmp0 <- get_COR(pos)
      tmp  <- tmp0$ORinterval
      
      L <-  round(tmp[1],1)
      C <-  round(tmp[3],1)
      R <-  round(tmp[2],1)
      
      text(top_BrS - 0.12,lat_pos + .3+gap, colnames(MBS_case_curr)[pos],cex=1 )
      text(top_BrS - 0.12,lat_pos + 1 / (2 * Jcause)+gap,C,cex = 1.5)
      text(top_BrS - 0.12,lat_pos - .2+gap,paste(c(L,"   ",R),collapse = " "),
           cex = 1.2)
      legend("topright",tmp0$label,bty = "n")
    }
    
    # x-axis for each cell:
    if (lat_pos>1){
      axis(1, seq(0,1,by = .05), lwd = 0, lwd.ticks = 0,#labels=rep("",length(seq(0,1,by=.2))),
           pos = seq(.625,Jcause +.625,by = 1)[lat_pos], cex.axis = 0.8,
           lty = 2,col = "blue"
      )
    }
  }
  
  Jcause <- length(model_options$likelihood$cause_list)
  #
  # plotting:
  #
  #op <- par(mar=c(5.1,4.1,4.1,0))
  op <- par(mar=c(5.1,0,4.1,0))
  
  if (!is_length_all_one(pos_vec)){
    #stop("== Not implemented for combo latent status.==")
    warning("== Combo latent status implemented with measurements overlapping in BrS columns! ==")
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
        if (first) {plot_BrS_cell(e,pos_curr,Jcause,gap = gap_seq[ct]);first <- FALSE}
        if (!first) {points_BrS_cell(e,pos_curr,Jcause,gap=gap_seq[ct])}
      }
    }
  }
  
  if (!is.null(bg_color) && !is.null(bg_color$BrS)){
    rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = 
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
          points_BrS_cell(e,pos_curr,Jcause,gap=gap_seq[ct])
        }
      }
    }
  }
  
#   #add axis labels on the left:
#   axis(2,at = c(sapply(1:Jcause,get_plot_num,height=Jcause)),
#        labels=rep(c("","case","ctrl"),Jcause),las=2)
#   axis(2,at=(1:Jcause)-.45,labels=rep("",Jcause),las=2,cex.axis=.5)
#   axis(2,at=(1:Jcause)-.35,labels=rep("",Jcause),las=2,cex.axis=.5)
#   
  #add ticks from 0 to 1 for x-bar:
  axis(1,at = c(0,0.2,0.4,0.6,0.8,1),labels= c(0,0.2,0.4,0.6,0.8,1),las=1)
  
  #add dashed lines to separate cells:
  abline(h=seq(1.5,Jcause-.5,by=1),lty=2,lwd=0.5,col="gray")
  abline(v=1,lty=2,lwd=.5,col="gray")
  
  #add some texts:
  mtext(eval(paste0("BrS: ", names(data_nplcm$Mobs$MBS)[slice])),
                   line=1,cex=1.8)
}

