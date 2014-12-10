#' Matrix of pairwise log odds ratios (LOR)
#'
#' The matrix of values is with cases on upper right and controls on lower
#' left. Log odds ratio is at the top of the cell.  Below it, its standard
#' error is in smaller type, using the same color as the LOR.  Then the
#' estimate is divided by its standard error. If it is less than 1 in
#' absolute value, put just a plus (red) or minus (blue) when
#' the Z-stat is between 1-2 in absolute value and put the actual
#' value when the Z is greater than 2.
#'
#' @param MBS.case Case bronze-standard measurements.
#' @param MBS.ctrl Control bronze-standard measurements.
#' @param pathogen_display The pathogen vector in desired order for display.
#' It can be of larger length than that of \code{pathogen_BrS}.
#' @param pathogen_BrS The string of pathogen names corresponding to each
#'  row/column in the plot. All pathogens have BrS measures.
#'  @param cex_main The size of LOR and Z-stat.
#'  @param cex_se The size of standard error of LOR.
#' @importFrom RColorBrewer brewer.pal
#' @return Figure of LOR matrix and relavent s.e. and significance information.
#' @export
#'
logORmat = function(MBS.case,MBS.ctrl,pathogen_display,pathogen_BrS,
                           cex_main=2,cex_se=1){
  # dat = datobs
  J <- ncol(MBS.case)
  Y   <- c(rep(1,nrow(MBS.case)),rep(0,nrow(MBS.ctrl)))
  Nd  <-  sum(Y)
  MBS <-  as.matrix(rbind(MBS.case,MBS.ctrl))
  Nu <-  length(Y)-Nd

  logORmat    <- matrix(NA,nrow=J,ncol=J)
  logORmat.se <- matrix(NA,nrow=J,ncol=J)

  for (j2 in 1:(J-1)){ #case (j2,j1); ctrl (j1,j2).
      for (j1 in (j2+1):J){

          ## cases:
          x = MBS.case[,j2]
          y = MBS.case[,j1]

          fit = glm(y~x,family = binomial(link="logit"))

          if ("x" %in% rownames(summary(fit)$coef)){
            logORmat[j2,j1] = round(summary(fit)$coef["x",1],3)
            logORmat.se[j2,j1] = round(summary(fit)$coef["x",2],3)
          }
          ##controls:
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
  cell.num.std.logOR = tmp

  tmp2 = logORmat.se
  tmp2[abs(logORmat.se)>10]=NA
  cell.num.prec = 1/tmp2^2

  colors <-rev(brewer.pal(10,"PuOr"))
  pal <- colorRampPalette(colors)

  val2col<-function(z, zlim, col = heat.colors(12), breaks){
    if(!missing(breaks)){
      if(length(breaks) != (length(col)+1)){stop("must have one more break than colour")}
    }
    if(missing(breaks) & !missing(zlim)){
      zlim[2] <- zlim[2]+c(zlim[2]-zlim[1])*(1E-3)#adds a bit to the range in both directions
      zlim[1] <- zlim[1]-c(zlim[2]-zlim[1])*(1E-3)
      breaks <- seq(zlim[1], zlim[2], length.out=(length(col)+1))
    }
    if(missing(breaks) & missing(zlim)){
      zlim <- range(z, na.rm=TRUE)
      zlim[2] <- zlim[2]+c(zlim[2]-zlim[1])*(1E-3)#adds a bit to the range in both directions
      zlim[1] <- zlim[1]-c(zlim[2]-zlim[1])*(1E-3)
      breaks <- seq(zlim[1], zlim[2], length.out=(length(col)+1))
    }
    CUT <- cut(z, breaks=breaks)
    colorlevels <- col[match(CUT, levels(CUT))] # assign colors to heights for each point
    return(colorlevels)
  }


  circle.cor = function(cor, cor.se, axes = FALSE, xlab = "",ylab = "",
                        asp = 1,title="",...) {
    n = nrow(cor)
    par(mar = c(0, 0, 5, 0), bg = "white")
    plot(c(0, n + 0.8), c(0, n + 0.8), axes = axes, xlab = "",
         ylab = "", asp = 1, type = "n")
    ##add grid
    segments(rep(0.5, n + 1), 0.5 + 0:n, rep(n + 0.5, n + 1),
             0.5 + 0:n, col = "gray")
    segments(0.5 + 0:n, rep(0.5, n + 1), 0.5 + 0:n, rep(n + 0.5,
                                                        n), col = "gray")
    ##define circles' background color.
    ##black for positive correlation coefficient and white for negative
    bg = matrix(0,nrow=n,ncol=n)

    for (i in 1:n){
      for (j in 1:n){
        if (!is.na(cor[i,j])){
          bg[i,j] = val2col(z=cor[i,j],c(-10,10),col=pal(64))
        }
      }
    }
    #bg = pal(64)
    #bg[cor > 0] = "black"
    #bg[cor <= 0] = "white"
    ##plot n*n circles using vector language, suggested by Yihui Xie

    #symbols(rep(1:n, each = n), rep(n:1, n), add = TRUE, inches = F,
    #        circles = as.vector(sqrt(abs(cor.se)/max(cor.se[!is.na(cor.se)]))/2),
    #        bg = as.vector(bg))

    text(rep(0, n), 1:n, n:1, col = "black")
    text(1:n, rep(n + 1), 1:n, col = "black")

    cor.txt<- round(t(cor)[,n:1],1)
    cor.se.txt <-round(t(cor.se)[,n:1],1)
    for (i in 1:n){
      for (j in 1:n){
        text(i,j+0.3,cor.txt[i,j],col=ifelse(cor.txt[i,j]>0,"red","blue"),cex=cex_main)
        text(i,j,cor.se.txt[i,j],col=ifelse(cor.txt[i,j]>0,"red","blue")
             ,cex=cex_se)
        abs.std.logOR <- abs(cor.txt[i,j]/sqrt(abs(cor.se.txt[i,j])))
        if (!is.na(abs.std.logOR) && abs.std.logOR>1){
          if (abs.std.logOR>2){
            text(i,j-0.3,round(cor.txt[i,j]/cor.se.txt[i,j],1),
                    col=ifelse(cor.txt[i,j]>0,"red","blue"),cex=cex_main)
          }else{
            text(i,j-0.3,ifelse(cor.txt[i,j]>0,"+","-"),
                 col=ifelse(cor.txt[i,j]>0,"red","blue"),cex=cex_main)
          }
        }

      }
    }
    segments(0.5,.5+n,.5+n,0.5,col="black",lty=3,lwd=3)
    mtext(title,3,cex=1,line=1)
  }

  circle.cor(cell.num.std.logOR,1/sqrt(cell.num.prec))

  pathogen_name <- pathogen_BrS[reorder(pathogen_display,pathogen_BrS)]
  for (s in rev(1:length(pathogen_name))){
    mtext(paste0(s,":",pathogen_name[s]),side=2,las=2,
          at=par("usr")[1]+0.03*(length(pathogen_name)-s)*diff(par("usr")[1:2])+5,
          cex=1,line=0)
  }
  mtext("cases",4,cex=2,line=6)
  mtext("controls",1,cex=2,line=-1)
  #mtext(paste(eda_options$X_names,eda_options$X_values,sep="="),3,cex=2,line=-0.5)

}

















