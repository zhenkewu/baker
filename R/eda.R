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
#' @param logOR_rounding Rounding number of the log odds ratio. Default is 2.
#' @param brs_slice Default is 1 - the set of BrS data to visualize.
#' 
#' @return Figure of LOR matrix and relavent s.e. and significance information.
#' @export

plot_logORmat = function(data_nplcm,
                         pathogen_display,
                         logOR_rounding = 2,
                         brs_slice = 1){
  Y <- data_nplcm$Y
  cat("== Visualizing pairwise log odds ratios for bronze-standard data set: ", brs_slice,". ==")
  MBS.case <- as.matrix(data_nplcm$Mobs$MBS[[brs_slice]][Y==1,,drop=FALSE])
  MBS.ctrl <- as.matrix(data_nplcm$Mobs$MBS[[brs_slice]][Y==0,,drop=FALSE])
  pathogen_BrS <- data_nplcm$Mname$Mname_BrS
  
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




















#         n1       <- nrow(MBS.case)
#         n0       <- nrow(MBS.ctrl)
#         J        <- ncol(MBS.case)

# #1. comparing the total no. of positives for cases and controls
# if (eda_options$total_positives == TRUE){
#   ct1 = rowSums(MBS.case)
#   ct0 = rowSums(MBS.ctrl)
#   tb1 = rep(NA,J+1);names(tb1)=c(0,1:J)
#   tb0 = tb1
#   for (j in 1:(J+1)){
#     tb1[j] = sum(ct1==j-1);tb0[j]=sum(ct0==j-1)
#   }
#
#   barplot(rbind(tb1/n1,tb0/n0),beside=TRUE,
#           #ylim=c(0,max(c(tb1/n1,tb0/n0))+.1),
#           legend.text = c("case", "control"),
#           xlab="no. of positives",
#           ylab="sample frequency",
#           main=paste(eda_options$X_names,eda_options$X_values,
#                      sep="="))
# }




