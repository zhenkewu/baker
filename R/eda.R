#' get an individual's data from the output of [clean_perch_data()]
#'
#' @param data_nplcm data for fitting nplcm; See [nplcm]
#' @param ID patient id: `patid`.
#'
#' @return a list with the inquired patient's data
#'
#' @examples 
#' \dontrun{
#' show_individual(data_nplcm,"G01431")
#' }
#' 
#' @family exploratory data analysis functions
#' @export
show_individual <- function(data_nplcm,ID) {
  index <- which(data_nplcm$X$patid == ID)
  Mobs  <- list()
  for (i in seq_along(data_nplcm$Mobs)) {
    Mobs[[i]] <- lapply(data_nplcm$Mobs[[i]],function(df)
      df[index,])
  }
  list(Mobs = Mobs,
       X    = data_nplcm$X[index,])
}

#' Visualize pairwise log odds ratios (LOR) for data that are available in
#' both cases and controls
#'
#' @details `plot_logORmat` visualizes a matrix of pairwise log odds ratios (LOR)
#'  for cases (upper) and controls (lower). LOR is at the top of the cell. 
#'  Below it, its standard error is in smaller type, using the same color as the LOR. 
#'  Then the estimate is divided by its standard error. We put the actual value when
#'   the Z-statistics has an absolute value greater than $2$; a plus (red) or minus (blue)
#'  if between $1$ and $2$; blank otherwise. 
#'
#' @param data_nplcm See [assign_model()].
#' @param pathogen_display The pathogen vector in desired order for display.
#' It can be of larger length than that of `pathogen_BrS`.
#' @param BrS_slice Default is 1 - the set of BrS data to visualize.
#' @param logOR_rounding Rounding number of the log odds ratio. Default is 2.
#' 
#' @return Figure of LOR matrix and relevant s.e. and significance information.
#' 
#' @family exploratory data analysis functions
#' @export

plot_logORmat = function(data_nplcm,
                         pathogen_display,
                         BrS_slice = 1,
                         logOR_rounding = 2){
  Y <- data_nplcm$Y
  cat("== Visualizing pairwise log odds ratios for bronze-standard data set: ", BrS_slice, ": ",names(data_nplcm$Mobs$MBS[BrS_slice]) ,". ==\n")
  MBS.case <- as.matrix(data_nplcm$Mobs$MBS[[BrS_slice]][Y==1,,drop=FALSE])
  MBS.ctrl <- as.matrix(data_nplcm$Mobs$MBS[[BrS_slice]][Y==0,,drop=FALSE])
  pathogen_BrS <- colnames(MBS.case)
  
  if(length(pathogen_BrS)==1){stop("== Cannot do log odds ratio plot with only one measurement! ==\n")}
  
  if (is.null(data_nplcm$Mobs$MBS) || is.na(data_nplcm$Mobs["MBS"])){
    stop("==No bronze-standard data!==")
  }
  
  # this plotting function will change graphical paramter "mar"; use
  # on.exit function to reset the graphical parameters after we have executed 
  # the function:
  default_par <- list(mar = graphics::par("mar"))
  on.exit(graphics::par(default_par))
  
  J   <- ncol(MBS.case)
  Nd  <-  sum(Y)
  Nu  <-  length(Y)-Nd
  
  logORmat    <- matrix(NA,nrow=J,ncol=J)
  logORmat.se <- matrix(NA,nrow=J,ncol=J)
  
  # reorder the data columns by pathogen_display; the original
  # order of the columns follows pathogen_BrS:
  index_display <- my_reorder(pathogen_display,pathogen_BrS)
  pathogen_name <- pathogen_BrS[index_display]
  MBS.case      <- MBS.case[,index_display]
  MBS.ctrl      <- MBS.ctrl[,index_display]
  MBS           <-  as.matrix(rbind(MBS.case,MBS.ctrl))
  
  for (j2 in 1:(J-1)){ #case (j2,j1); ctrl (j1,j2).
    for (j1 in (j2+1):J){
      
      ## cases: (upper)
      x = MBS.case[,j2]
      y = MBS.case[,j1]
      
      fit = stats::glm(y~x,family = stats::binomial(link="logit"))
      
      if ("x" %in% rownames(summary(fit)$coef)){
        logORmat[j2,j1] = round(summary(fit)$coef["x",1],3)
        logORmat.se[j2,j1] = round(summary(fit)$coef["x",2],3)
      }
      ##controls: (lower)
      x = MBS.ctrl[,j2]
      y = MBS.ctrl[,j1]
      
      fit = stats::glm(y~x,family = stats::binomial(link="logit"))
      
      if ("x" %in% rownames(summary(fit)$coef)){
        logORmat[j1,j2] = round(summary(fit)$coef["x",1],3)
        logORmat.se[j1,j2] = round(summary(fit)$coef["x",2],3)
      }
    }
  }
  
  #cell.num = logORmat/logORmat.se
  tmp       = logORmat
  tmp[abs(logORmat.se)>10]=NA
  cell.num.logOR = tmp
  
  tmp2 = logORmat.se
  tmp2[abs(logORmat.se)>10]=NA
  cell.num.prec = 1/tmp2^2
  
  #   cor = cell.num.logOR
  #   cor.se = 1/sqrt(cell.num.prec)
  
  circle.cor = function(cor, cor.se, axes = FALSE, xlab = "",ylab = "",
                        asp = 1,title="",...) {
    n = nrow(cor)
    # size of the numbers in the boxes:
    cex_main= min(2,20/n)
    cex_se  = min(1.5,15/n)
    
    graphics::par(mar = c(0, 1, 5, 0), bg = "white",xpd=TRUE)
    graphics::plot(c(0, n + 0.8), c(0, n + 0.8), axes = axes, xlab = "",
         ylab = "", asp = 1, type = "n")
    ##add grid
    graphics::segments(rep(0.5, n + 1), 0.5 + 0:n, rep(n + 0.5, n + 1),
             0.5 + 0:n, col = "gray")
    graphics::segments(0.5 + 0:n, rep(0.5, n + 1), 0.5 + 0:n, rep(n + 0.5,
                                                        n), col = "gray")
    cor.txt<- round(t(cor)[,n:1],logOR_rounding)
    cor.se.txt <-round(t(cor.se)[,n:1],logOR_rounding)
    cor.txt3<- round(t(cor)[,n:1],3)
    cor.se.txt3 <-round(t(cor.se)[,n:1],3)
    
    graphics::text(1,J+0.3,"logOR",cex=cex_main/2)
    graphics::text(1,J,"s.e.",cex=cex_se/2)
    graphics::text(1,J-0.3,"std.logOR",cex=cex_se/2)
    
    for (i in 1:n){
      for (j in 1:n){
        graphics::text(i,j+0.3,cor.txt[i,j],col=ifelse(cor.txt[i,j]>0,"red","blue"),cex=cex_main)
        graphics::text(i,j,cor.se.txt[i,j],col=ifelse(cor.txt[i,j]>0,"red","blue")
             ,cex=cex_se)
        abs.std.logOR <- abs(cor.txt3[i,j]/cor.se.txt3[i,j])
        if (!is.na(abs.std.logOR) && abs.std.logOR>1){
          if (abs.std.logOR>2){
            graphics::text(i,j-0.3,round(cor.txt3[i,j]/cor.se.txt3[i,j],1),
                 col=ifelse(cor.txt[i,j]>0,"red","blue"),cex=cex_se)
          }else{
            graphics::text(i,j-0.3,ifelse(cor.txt[i,j]>0,"+","-"),
                 col=ifelse(cor.txt[i,j]>0,"red","blue"),cex=cex_se)
          }
        }
        
      }
    }
    # diagonal line:
    graphics::segments(0.5+1,.5+n-1,.5+n,0.5,col="black",lty=3,lwd=3)
    graphics::mtext(title,3,cex=1,line=1)
  }
  
  # put texts in the boxes:
  circle.cor(cell.num.logOR,1/sqrt(cell.num.prec))
  
  # put pathogen names on rows and columns:
  for (s in rev(1:length(pathogen_name))){
    #   graphics::text(-2,graphics::par("usr")[1]+0.03*(length(pathogen_name)-s)*diff(graphics::par("usr")[1:2])+5,
    #        paste0(s,":",pathogen_name[s]),las=2,
    #        cex=1)
    graphics::text(0.25,J-s+1,paste0(pathogen_name[s],":(",s,")"),cex=min(1.5,20/J),srt=45,adj=1)
    graphics::text(s,J+0.7,paste0("(",s,"):",pathogen_name[s]),cex=min(1.5,20/J),srt=45,adj=0)
  }
  # labels for cases and controls:
  graphics::text(J+1,J/2,"cases",cex=2,srt=-90)
  graphics::text(J/2,0,"controls",cex=2)
}




#' silver-standard data summary
#' 
#' @param SS_dat a data frame of silver-standard data. It can usually 
#' be obtained by `data_nplcm$Mobs$MSS[[1]]`, meaning the first SS measurement
#' slice.
#' 
#' @param Y a vector of case control status: 1 for case; 0 for control.
#' 
#' @return a vector of number of positives
#' 
#' @examples 
#' \dontrun{
#' summarize_SS(data_nplcm$Mobs$MSS[[1]], data_nplcm$Y)
#' }
#' 
#' @family exploratory data analysis functions
#' @export
summarize_SS <- function(SS_dat, Y){
  # get observed rates' summaries:
  Nd <- sum(Y==1)
  Nu <- sum(Y==0)
  N  <- Nd+Nu
  
  MSS_curr <- SS_dat
  MSS_nm <- colnames(SS_dat)
  # positive rates and confidence intervals:
  #cases:
  MSS_case_curr <- MSS_curr[1:Nd,,drop=FALSE]
  count    <- as.integer(do.call(cbind,lapply(MSS_case_curr,sum,na.rm=TRUE))) #<-- added 'as.integer' to make pathogens appear by rows.
  NA_count <- apply(MSS_case_curr,2,function(v) sum(is.na(v)))
  tmp.case <- binom.confint(count,Nd-NA_count,conf.level = 0.95, methods = "ac")
  
  # case and control positive rate, lower and upper limit
  MSS_mean  <- rbind(round(tmp.case$mean,5))
  MSS_q1    <- rbind(tmp.case[,c("lower")])
  MSS_q2    <- rbind(tmp.case[,c("upper")])
  
  make_list(count, NA_count, Nd, MSS_mean, MSS_q1, MSS_q2, MSS_nm)
}



#' summarize bronze-standard data
#' 
#' @param BrS_dat bronze-standard data, which is usually `data_nplcm$Mobs$MBS[[1]]`
#' 
#' @param Y A vector of case/control status: 1 for case; 0 for control
#' 
#' @return a list of summaries for BrS data
#' @examples 
#' \dontrun{
#' summarize_BrS(data_nplcm$Mobs$MBS[[1]], data_nplcm$Y)
#' }
#' 
#' @family exploratory data analysis functions
#' @export

summarize_BrS <- function(BrS_dat,Y){
  # get observed rates' summaries:
  Nd <- sum(Y==1)
  Nu <- sum(Y==0)
  N  <- Nd+Nu
  
  MBS_curr <- BrS_dat
  MBS_nm <- colnames(BrS_dat)
  # positive rates and confidence intervals:
  #cases:
  
  MBS_case_curr <- MBS_curr[1:Nd,,drop=FALSE]
  count    <- as.integer(do.call(cbind,lapply(MBS_case_curr,sum,na.rm=TRUE)))
  NA_count <- apply(MBS_case_curr,2,function(v) sum(is.na(v)))
  tmp.case <- binom.confint(count,Nd-NA_count,conf.level = 0.95, methods = "ac")
  
  #controls:
  MBS_ctrl_curr <- MBS_curr[-(1:Nd),,drop=FALSE]
  count    <- as.integer(do.call(cbind,lapply(MBS_ctrl_curr,sum,na.rm=TRUE)))
  NA_count <- apply(MBS_ctrl_curr,2,function(v) sum(is.na(v)))
  tmp.ctrl <- binom.confint(count, Nu-NA_count, conf.level = 0.95, methods = "ac")
  
  # case and control positive rate, lower and upper limit
  MBS_mean  <- rbind(round(tmp.case$mean,5),round(tmp.ctrl$mean,5));rownames(MBS_mean) = c("case","control")
  MBS_q1 <- rbind(tmp.case[,c("lower")],tmp.ctrl[,c("lower")]); rownames(MBS_q1) = c("case","control")
  MBS_q2 <- rbind(tmp.case[,c("upper")],tmp.ctrl[,c("upper")]); rownames(MBS_q2) = c("case","control")

  make_list(Nd, MBS_mean, MBS_q1, MBS_q2, MBS_nm)
}


#' get top patterns from a slice of bronze-standard measurement
#' 
#' @param BrS_dat bronze-standard data, which is usually `data_nplcm$Mobs$MBS[[1]]`
#' 
#' @param Y A vector of case/control status: 1 for case; 0 for control
#' @param case_status 1 for case; 0 for controls
#' @param n_pat the number of top patterns one wants to show
#' @param exclude_missing DEFAULT is TRUE for excluding any individual with missing measurements.
#' 
#' @examples 
#' \dontrun{
#' 
#' res <- get_top_pattern(data_nplcm$Mobs$MBS[[1]],data_nplcm$Y,1,30,FALSE)
#' }
#' 
#' @return a list of results: `obs_pat` - observed rates; 
#' `pattern_names`; `exist_other` - if
#' actual no. of patterns is larger than `n_pat`; `N`- No. of individuals
#' with `Y = case_status`.
#' 
#' @family exploratory data analysis functions
#' @export
get_top_pattern <- function(BrS_dat,Y,case_status,n_pat,exclude_missing = TRUE){
  # getting data:
  
  # get observed data:
  curr_observed <- BrS_dat
  # length of a pattern (e.g., 10001 means length is 5)
  len_pat       <- ncol(curr_observed)
  curr_Y        <- Y # case control status.
  
  # subsetting into case or control based on input `case_status`:
  observed  <- curr_observed[curr_Y==case_status,,drop=FALSE]
  N <- sum(curr_Y==case_status)
  
  # convert numeric vector into a character string: e.g., c(1,0,0,1,1) into "10011"; for faster pattern matching:
  collapse_byrow <- function(mat){
    NA2dot(apply(mat,1,paste,collapse = "" ))
  }
  
  observed_pat  <- collapse_byrow(observed)
  
  # counting patterns:
  pat               <- sort(table(observed_pat),decreasing=TRUE) # observed patten.
  n_pat_used        <- min(n_pat,length(pat)) # actually used pattern number.
  pat_high_frac     <- pat[1:n_pat_used]/length(observed_pat)
  pattern_names     <- names(pat_high_frac)
  obs_pat           <- pat_high_frac
  exist_other       <- (length(pat)-n_pat_used) >0
  if (exist_other){
    pattern_names  <- c(names(pat_high_frac),"other")
    obs_pat        <- c(pat_high_frac,1-sum(pat_high_frac))
  }
  
  
  if (exclude_missing){
    ind_missing       <- grep("\\.",names(pat)) # pick out patterns with missing measurements.
    n_missing         <- sum(pat[ind_missing])  # the total number of individuals with missing measurements.
    pat_high_frac     <- pat[1:n_pat_used]/(length(observed_pat)-n_missing) # divide by no. of individuals with complete measurements.
    pat_high_frac_no_missing     <- pat_high_frac
    if (length(ind_missing)){
      pat_high_frac_no_missing     <- pat_high_frac[-ind_missing] # delete patterns with missingness.
    }
    pat_high_name_no_missing     <- names(pat_high_frac_no_missing) # get names to display in the plot.
    
    n_pat_used_no_missing <- length(pat_high_name_no_missing) # the length of patterns without missingness.
    
    exist_other <- (length(pat)-length(ind_missing))>0
    pattern_names        <- c(pat_high_name_no_missing)
    obs_pat              <- pat_high_frac_no_missing
    if (exist_other){
      pattern_names  <- c(pat_high_name_no_missing,"other")
      obs_pat        <- c(pat_high_frac_no_missing,1-sum(pat_high_frac_no_missing))
    }
  }
  names(obs_pat) <- pattern_names
  make_list(obs_pat,pattern_names,exist_other,N)
}

#' visualize trend of pathogen observation rate for NPPCR data (both cases and controls)
#' 
#' @details This function shows observed
#' positive rate for continuous covariates,e.g., enrollment date 
#' in PERCH application. Smoothing is done by penalized splines implemented by 
#' `mgcv` package. The penalized spline smoothing term is constructed by 
#' [mgcv::smooth.construct.ps.smooth.spec()]
#' 
#' @param data_nplcm Data set produced by [clean_perch_data()]
#' @param patho the index of pathogen
#' @param slice the slice of BrS data for visualization; default is 1.
#' @param slice_SS the slice of SS data to add onto BrS plots; default is 1, usually
#' representing blood culture measurements.
#' 
#' @import mgcv
#' 
#' @return A figure with smoothed positive rate and confidence bands for cases
#' and controls, respectively. The right margin shows marginal positive rates.
#' 
#' 
#' @family exploratory data analysis functions
#' 
#' @export

visualize_season <- function(data_nplcm, patho, slice = 1,slice_SS = 1){
  ord_all <- order(data_nplcm$X$ENRLDATE)
  
  Y <- data_nplcm$Y[ord_all]
  X <- data_nplcm$X[ord_all,]
  Nd <- sum(Y)
  Nu <- sum(1-Y)
  curr_MBS <- data_nplcm$Mobs$MBS[[slice]][ord_all,]
  # dataframe with enrollment dates in the last column: Rdate
  curr_dat <- cbind(curr_MBS,Rdate= X$ENRLDATE)
  
  smooth_season <- function(){
    # prepare data by centering and transformation:
    std_date      <- dm_Rdate_FPR(curr_dat$Rdate,Y,effect = "fixed")
    # case fit and predict:
    pred_date1    <- seq(min(std_date[Y==1]),max(std_date[Y==1]),length=1000)
    datcase_path  <- data.frame(M = curr_MBS[Y==1,patho],Rdate = std_date[Y==1])
    fit_case      <- gam(M~s(Rdate,bs="ps",k=10,m=c(2,1)),  
                         data = datcase_path,family=stats::binomial(logit),method="REML")
    
    fitted_case  <- mgcv::predict.gam(fit_case,type="response")[1:Nd]
    pred_case    <- mgcv::predict.gam(fit_case,data.frame(Rdate=pred_date1), type = "link", se.fit = TRUE)
    
    # control fit and predict:
    pred_date0   <- seq(min(std_date[Y==0]),max(std_date[Y==0]),length=1000)
    datctrl_path <- data.frame(M = curr_MBS[Y==0,patho],Rdate = std_date[Y==0])
    fit_ctrl     <- gam(M~s(Rdate,bs="ps",k=10,m=c(2,1)), 
                        data = datctrl_path,family=stats::binomial(logit),method="REML")
    
    fitted_ctrl  <- mgcv::predict.gam(fit_ctrl,type="response")[1:Nu]
    pred_ctrl    <- mgcv::predict.gam(fit_ctrl,data.frame(Rdate=pred_date0), type = "link", se.fit = TRUE)
    
    make_list(fitted_case,fitted_ctrl,fit_case,fit_ctrl,pred_case,pred_ctrl,
              pred_date1,pred_date0,std_date)
  }
  
  # patho specifies which pathogen to smooth:
  response.case <- curr_MBS[Y==1,patho]
  response.ctrl <- curr_MBS[Y==0,patho]
  # calculate the smoothed curve of pathogen detection over seasons:
  out <- smooth_season()  
  fitted_case <- out$fitted_case
  fitted_ctrl <- out$fitted_ctrl
  fit_case   <- out$fit_case
  fit_ctrl   <- out$fit_ctrl
  pred_case  <- out$pred_case
  pred_ctrl  <- out$pred_ctrl
  pred_date1 <- out$pred_date1
  pred_date0 <- out$pred_date0
  std_date   <- out$std_date
  
  # because pred_date1/0 are standardized, transform back to original scale:
  pred.date.case.plot <- pred_date1*stats::sd(curr_dat$Rdate[Y==0])+mean(curr_dat$Rdate[Y==0])
  pred.date.ctrl.plot <- pred_date0*stats::sd(curr_dat$Rdate[Y==0])+mean(curr_dat$Rdate[Y==0])
  
  #
  # plotting raw data:
  #
  # some date transformations:
  X$date_plot  <- as.Date(X$ENRLDATE)
  X$date_month_centered <- as.Date(cut(X$date_plot,breaks="2 months"))+30
  X$date_month <- as.Date(cut(X$date_plot,breaks="2 months"))
  
  color2 <- grDevices::rgb(190, 190, 190, alpha=200, maxColorValue=255)
  color1 <- grDevices::rgb(216,191,216, alpha=200, maxColorValue=255)
  #
  #cases:
  #
  
  dat_case <- cbind(X[Y==1,],fitted_case)
  # fitted curve:
  graphics::plot(fitted_case ~ date_plot, data=dat_case,
       type = "l",lwd=4,
       xlab = "",xaxt = "n", axes=F,xlim=c(min(X$date_plot),max(X$date_plot)+80),
       ylab = colnames(curr_MBS)[patho],ylim=c(-.2,1.4))
  
  
  last_interval <- max(X$date_month)
  lubridate::month(last_interval) <- lubridate::month(last_interval) +2
  graphics::axis(1, c(X$date_month,last_interval), format(c(X$date_month,last_interval), "%Y %b"), 
       cex.axis = .7)
  graphics::axis(2,at = seq(0,1,by=0.2),labels=seq(0,1,by=0.2))
  
  
  #graphics::points(upr~as.Date(pred.date.case.plot),lty=2,type="l",lwd=2)
  #graphics::points(lwr~as.Date(pred.date.case.plot),lty=2,type="l",lwd=2)
  
  #rug plot:
  graphics::points(dat_case$date_plot,c(-0.1,1.15)[response.case+1],pch="|")
  
  #
  # controls:
  #
  dat_ctrl <- cbind(X[Y==0,],fitted_ctrl)
  # get value in linear predictor scale:
  upr <- pred_ctrl$fit + (1.96 * pred_ctrl$se.fit)
  lwr <- pred_ctrl$fit - (1.96 * pred_ctrl$se.fit)
  # transform to the right scale:
  upr <- fit_ctrl$family$linkinv(upr)
  lwr <- fit_ctrl$family$linkinv(lwr)
  graphics::polygon(c(as.Date(pred.date.ctrl.plot), rev(as.Date(pred.date.ctrl.plot))),
          c(lwr,rev(upr)),col=color2,border=NA)
  #plot control actual data:
  graphics::points(fitted_ctrl ~ date_plot,data=dat_ctrl,
         pch=2,cex=2,col="dodgerblue2",lwd=5,type="l",lty=1)
  
  ma <- function(x,n=60){stats::filter(x,rep(1/n,n), sides=2)}
  
  dat_ctrl$runmean <- ma(response.ctrl)
  graphics::points(runmean ~ date_plot,data=dat_ctrl,lty=2,pch=1,cex=0.5,type="o",col="dodgerblue2")
  
  
  # raw moving-window prevalences:
  # interval_raw <- aggregate(response.ctrl~dat_ctrl$date_month_centered, FUN=mean)
  # graphics::points(interval_raw[,2]~as.Date(interval_raw[,1]),
  #        lty=2,pch=1,cex=0.7,col="dodgerblue2",type="o",
  #        lwd=1)
  
  #### plot some case data here to prevent overlapping:
  # get value in linear predictor scale:
  upr <- pred_case$fit + (1.96 * pred_case$se.fit)
  lwr <- pred_case$fit - (1.96 * pred_case$se.fit)
  # transform back to the right scale:
  upr <- fit_case$family$linkinv(upr)
  lwr <- fit_case$family$linkinv(lwr)
  graphics::polygon(c(as.Date(pred.date.case.plot), rev(as.Date(pred.date.case.plot))),
          c(lwr,rev(upr)),col=color1,border=NA)
  
  graphics::points(fitted_case ~ date_plot, data=dat_case,
         type = "l",lwd=4)
  
  # raw moving-window prevalences:
  #interval_raw <- aggregate(response.case~dat_case$date_month_centered, FUN=mean)
  #graphics::points(interval_raw[,2]~as.Date(interval_raw[,1]),
  #       lty=1,pch=20,cex=0.7,type="o")
  
  dat_case$runmean <- ma(response.case)
  graphics::points(runmean ~ date_plot,data=dat_case,lty=1,pch=20,cex=0.5)
  
  
  #graphics::points(upr~as.Date(pred.date.ctrl.plot),lty=2,type="l",col="dodgerblue2",lwd=2)
  #graphics::points(lwr~as.Date(pred.date.ctrl.plot),lty=2,type="l",col="dodgerblue2",lwd=2)
  # rug plot:
  graphics::points(dat_ctrl$date_plot,c(-0.2,1.05)[response.ctrl+1],pch="|", col="dodgerblue2")
  
  
  #
  # plot overall marginal mean:
  #
  delta            <- 40
  case.overall.loc <- max(X$date_plot)+10
  ctrl.overall.loc <- case.overall.loc+delta
  
  graphics::text(case.overall.loc,mean(response.case,na.rm=TRUE)+0.3,
       paste0(round(mean(response.case)*100,1),"%"),col="black",pch=20,srt=90,cex=2)
  graphics::text(ctrl.overall.loc,mean(response.ctrl,na.rm=TRUE)+0.3,
       paste0(round(mean(response.ctrl)*100,1),"%"),col="dodgerblue2",pch=20,srt=90,cex=2)
  
  ncase = sum(!is.na(response.case))
  nctrl = sum(!is.na(response.ctrl))
  tmp.case = binom::binom.confint(sum(response.case,na.rm=TRUE), ncase, conf.level = 0.95, methods = "ac")
  tmp.ctrl = binom::binom.confint(sum(response.ctrl,na.rm=TRUE), nctrl, conf.level = 0.95, methods = "ac")
  Bcomp = rbind(round(tmp.case$mean,5),round(tmp.ctrl$mean,5))
  Bcompq1 = rbind(tmp.case[,c("lower")],tmp.ctrl[,c("lower")])
  Bcompq2 = rbind(tmp.case[,c("upper")],tmp.ctrl[,c("upper")])
  
  graphics::points(case.overall.loc,mean(response.case),col="black",pch=20,cex=2)
  graphics::points(ctrl.overall.loc,mean(response.ctrl),col="dodgerblue2",pch=20,cex=2)
  
  graphics::segments(case.overall.loc,mean(response.case),ctrl.overall.loc,mean(response.ctrl))
  graphics::segments(case.overall.loc,Bcompq1[1,],case.overall.loc,Bcompq2[1,])
  graphics::segments(ctrl.overall.loc,Bcompq1[2,],ctrl.overall.loc,Bcompq2[2,],col="dodgerblue2")
  
  ## plot silver-standard data
  
  if (!is.null(data_nplcm$Mobs$MSS) && 
      colnames(curr_MBS)[patho]%in%colnames(data_nplcm$Mobs$MSS[[slice_SS]])){
    
    patho_SS <- which(colnames(data_nplcm$Mobs$MSS[[slice_SS]])==colnames(curr_MBS)[patho])
    
    curr_MSS <- (data_nplcm$Mobs$MSS[[slice_SS]][ord_all,,drop=FALSE])[data_nplcm$Y[ord_all]==1,,drop=FALSE]
    ind_BrS_not_missing <- which(!is.na(curr_MSS[,patho_SS]))
    graphics::points(dat_case$date_plot[ind_BrS_not_missing],
           c(-0.5,1.3)[curr_MSS[ind_BrS_not_missing,patho_SS]+1],pch="|",col="red",lwd=2)
    graphics::mtext(
         paste0(names(data_nplcm$Mobs$MSS)[slice_SS],"(+)%:",round(100*mean(curr_MSS[ind_BrS_not_missing,patho_SS]),2),"%"),
         cex=0.6,col=2,at=1.3,side=2,las=1)
  }
  
  
}


