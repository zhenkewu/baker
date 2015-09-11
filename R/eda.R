#' Visualize pairwise log odds ratios (LOR) for data that are available in
#' both cases and controls
#'
#' @details \code{plot_logORmat} visualizes a matrix of values (log odds ratios here). 
#' Cases' values are on upper right and controls on lower
#' left. Log odds ratio (LOR) is at the top of the cell.  Below it, its standard
#' error is in smaller type, using the same color as the LOR.  Then the
#' estimate is divided by its standard error. If it is less than 1 in
#' absolute value, we put a plus (red) or minus (blue) when
#' the Z-stat is between 1-2 in absolute value and put the actual
#' value when the Z is greater than 2.
#'
#' @param data_nplcm See \code{\link{assign_model}}.
#' @param pathogen_display The pathogen vector in desired order for display.
#' It can be of larger length than that of \code{pathogen_BrS}.
#' @param BrS_slice Default is 1 - the set of BrS data to visualize.
#' @param logOR_rounding Rounding number of the log odds ratio. Default is 2.
#' 
#' @return Figure of LOR matrix and relavent s.e. and significance information.
#' @export

plot_logORmat = function(data_nplcm,
                         pathogen_display,
                         BrS_slice = 1,
                         logOR_rounding = 2){
  Y <- data_nplcm$Y
  cat("== Visualizing pairwise log odds ratios for bronze-standard data set: ", BrS_slice, ": ",names(data_nplcm$Mobs$MBS[BrS_slice]) ,". ==")
  MBS.case <- as.matrix(data_nplcm$Mobs$MBS[[BrS_slice]][Y==1,,drop=FALSE])
  MBS.ctrl <- as.matrix(data_nplcm$Mobs$MBS[[BrS_slice]][Y==0,,drop=FALSE])
  pathogen_BrS <- colnames(MBS.case)
  
  if(length(pathogen_BrS)==1){stop("== Cannot do log odds ratio plot with only one measurement! ==")}
  
  if (is.null(data_nplcm$Mobs$MBS) || is.na(data_nplcm$Mobs["MBS"])){
    stop("==No bronze-standard data!==")
  }
  
  # this plotting function will change graphical paramter "mar"; use
  # on.exit function to reset the graphical parameters after we have executed 
  # the function:
  default_par <- list(mar = par("mar"))
  on.exit(par(default_par))
  
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
      
      fit = glm(y~x,family = binomial(link="logit"))
      
      if ("x" %in% rownames(summary(fit)$coef)){
        logORmat[j2,j1] = round(summary(fit)$coef["x",1],3)
        logORmat.se[j2,j1] = round(summary(fit)$coef["x",2],3)
      }
      ##controls: (lower)
      x = MBS.ctrl[,j2]
      y = MBS.ctrl[,j1]
      
      fit = glm(y~x,family = binomial(link="logit"))
      
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
    
    par(mar = c(0, 0, 5, 0), bg = "white",xpd=TRUE)
    plot(c(0, n + 0.8), c(0, n + 0.8), axes = axes, xlab = "",
         ylab = "", asp = 1, type = "n")
    ##add grid
    segments(rep(0.5, n + 1), 0.5 + 0:n, rep(n + 0.5, n + 1),
             0.5 + 0:n, col = "gray")
    segments(0.5 + 0:n, rep(0.5, n + 1), 0.5 + 0:n, rep(n + 0.5,
                                                        n), col = "gray")
    cor.txt<- round(t(cor)[,n:1],logOR_rounding)
    cor.se.txt <-round(t(cor.se)[,n:1],logOR_rounding)
    cor.txt3<- round(t(cor)[,n:1],3)
    cor.se.txt3 <-round(t(cor.se)[,n:1],3)
    
    text(1,J+0.3,"logOR",cex=cex_main/2)
    text(1,J,"s.e.",cex=cex_se/2)
    text(1,J-0.3,"std.logOR",cex=cex_se/2)
    
    for (i in 1:n){
      for (j in 1:n){
        text(i,j+0.3,cor.txt[i,j],col=ifelse(cor.txt[i,j]>0,"red","blue"),cex=cex_main)
        text(i,j,cor.se.txt[i,j],col=ifelse(cor.txt[i,j]>0,"red","blue")
             ,cex=cex_se)
        abs.std.logOR <- abs(cor.txt3[i,j]/cor.se.txt3[i,j])
        if (!is.na(abs.std.logOR) && abs.std.logOR>1){
          if (abs.std.logOR>2){
            text(i,j-0.3,round(cor.txt3[i,j]/cor.se.txt3[i,j],1),
                 col=ifelse(cor.txt[i,j]>0,"red","blue"),cex=cex_se)
          }else{
            text(i,j-0.3,ifelse(cor.txt[i,j]>0,"+","-"),
                 col=ifelse(cor.txt[i,j]>0,"red","blue"),cex=cex_se)
          }
        }
        
      }
    }
    # diagonal line:
    segments(0.5+1,.5+n-1,.5+n,0.5,col="black",lty=3,lwd=3)
    mtext(title,3,cex=1,line=1)
  }
  
  # put texts in the boxes:
  circle.cor(cell.num.logOR,1/sqrt(cell.num.prec))
  
  # put pathogen names on rows and columns:
  for (s in rev(1:length(pathogen_name))){
    #   text(-2,par("usr")[1]+0.03*(length(pathogen_name)-s)*diff(par("usr")[1:2])+5,
    #        paste0(s,":",pathogen_name[s]),las=2,
    #        cex=1)
    text(-0,J-s+1,paste0(pathogen_name[s],":(",s,")"),cex=min(1.5,20/J),adj=1)
    text(s,J+0.7,paste0("(",s,"):",pathogen_name[s]),cex=min(1.5,20/J),srt=45,adj=0)
  }
  # labels for cases and controls:
  text(J+1,J/2,"cases",cex=2,srt=-90)
  text(J/2,0,"controls",cex=2)
}




#' silver-standard data summary
#' 
#' @param SS_dat a data frame of silver-standard data. It can usually 
#' be obtained by \code{data_nplcm$Mobs$MSS[[1]]}, meaning the first SS measurement
#' slice.
#' 
#' @param Y a vector of case control status: 1 for case; 0 for control.
#' 
#' @return a vector of number of positives
#' 
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

#summarize_SS(data_nplcm$Mobs$MSS[[1]], data_nplcm$Y)

#' summarize bronze-standard data
#' 
#' @param BrS_dat bronze-stanadrd data, which is usually \code{data_nplcm$Mobs$MBS[[1]]}
#' 
#' @param Y A vector of case/control status: 1 for case; 0 for control
#' 
#' @return a list of summaries for BrS data
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
#' @param BrS_dat bronze-stanadrd data, which is usually \code{data_nplcm$Mobs$MBS[[1]]}
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
#' @return a list of results: \code{obs_pat} - observed rates; 
#' \code{pattern_names}; \code{exist_other} - if
#' actual no. of patterns is larger than \code{n_pat}; \code{N}- No. of individuals
#' with \code{Y = case_status}.
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


