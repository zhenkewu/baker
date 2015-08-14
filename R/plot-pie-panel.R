#' Plot etiology (pie) panel
#' 
#' 
#' @param model_options See \code{\link{nplcm}}
#' @param res_nplcm See \code{\link{nplcm_read_folder}}
#' @param bugs.dat Data input for the model fitting.
#' @param bg_color A list with names "BrS", "SS", "pie" to specify background colors
#' @param top_pie Numerical value to specify the rightmost limit 
#' on the horizontal axis for the pie panel.
#' @param label_size the size of latent status labels on the right margin
#' 
#' @importFrom binom binom.confint
#' 
#' @export

plot_pie_panel <- function(model_options,
                           res_nplcm,
                           bugs.dat,
                           bg_color,
                           top_pie = 1,
                           label_size = 1 ){

  # order cause_list by posterior means:
  ord <- order_post_eti(res_nplcm,model_options)$ord
  pEti_mat_ord <- order_post_eti(res_nplcm,model_options)$pEti_mat_ord
  pEti_mean_ord <- colMeans(pEti_mat_ord)
  
  # quantiles for etiology: outer is 97.5% CI, inner is 50% CI
  pEti_q   <- apply(pEti_mat_ord,2,quantile,probs=c(0.025,0.975,0.25,0.75))
  
  cause_list <- model_options$likelihood$cause_list
  cause_list_ord <- cause_list[ord]
  
  Nd <- bugs.dat$Nd
  Nu <- bugs.dat$Nu
  
  Jcause <- length(model_options$likelihood$cause_list)
  alpha_ord <- bugs.dat$alpha[ord]
  
  plot_pie_cell_first <- function(lat_pos,height,dotcolor="black",add=FALSE){
    # posterior mean of etiology:
    if (!add){
    plot(pEti_mean_ord[lat_pos],lat_pos,
         yaxt="n",
         xlim=c(0,top_pie),ylim=c(0.5,Jcause+0.5),
         col="purple",
         ylab="",xlab="probability",
         pch= 20,cex=2)
    }
    if (add){
      points(pEti_mean_ord[lat_pos],lat_pos,
           yaxt="n",
           xlim=c(0,top_pie),ylim=c(0.5,Jcause+0.5),
           col="purple",
           ylab="",xlab="probability",
           pch= 20,cex=2)
      
    }
  }
  
  points_pie_cell <- function(lat_pos,height,dotcolor="black"){
    if (lat_pos > 1){
      # posterior mean of etiology:
      points(pEti_mean_ord[lat_pos],lat_pos,
           yaxt="n",
           xlim=c(0,top_pie),ylim=c(lat_pos-0.5,lat_pos+0.5),
           col="purple",
           ylab="",xlab="probability",
           pch= 20,cex=2)
    }
    # x-axis for each cell:
    if (lat_pos>1){
      axis(1, seq(0,1,by = .2), lwd = 0, lwd.ticks = 0,#labels=rep("",length(seq(0,1,by=.2))),
           pos = seq(.625,Jcause +.625,by = 1)[lat_pos], cex.axis = 0.8,lty =
             2,col = "blue"
      )
    }
    
    points(pEti_q[,lat_pos],rep(lat_pos,4),pch=c("|","|","[","]"),cex=1,col="purple")
    
    pgrid = seq(0,1,by=0.01)
    segments(y0=lat_pos,x0=pEti_q[1,lat_pos],y1=lat_pos,x1=pEti_q[2,lat_pos],col=dotcolor)
    segments(y0=lat_pos,x0=pEti_q[3,lat_pos],y1=lat_pos,x1=pEti_q[4,lat_pos],col=dotcolor, lwd=2)
    text(.8,lat_pos,paste0("=",paste0(round(100*c(pEti_mean_ord),1)[lat_pos],"%")),srt=0,cex=2)
    text(.65,lat_pos,bquote(hat(pi)[.(ord[lat_pos])]),srt=0,cex=2)
    #prior density:
    tmp.density = dbeta(pgrid,alpha_ord[lat_pos],sum(alpha_ord[-lat_pos]))
    points(pgrid,tmp.density/(3*max(tmp.density))+lat_pos-0.45,type="l",col="gray",lwd=4,lty=2)
    #posterior density:
    tmp.post.density = density(pEti_mat_ord[,lat_pos],from=0,to=1)
    tmp.x = tmp.post.density$x
    tmp.y = tmp.post.density$y
    points(tmp.x,tmp.y/(3*max(tmp.y))+lat_pos-0.45,lwd=4,type="l",col="purple")
    
  }
  
  #
  # plot etiology information:
  #
  op <- par(mar=c(5.1,0,4.1,10))
  
  cat("\n == Plotting pies == ")
  
  first <- TRUE
  for (e in 1:Jcause){
    if (first){plot_pie_cell_first(e,Jcause)}
    points_pie_cell(e,Jcause)
    first <- FALSE
  }
  
  if (!is.null(bg_color) && !is.null(bg_color$pie)){
    rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = 
           bg_color$pie)
    
    first <- TRUE
    for (e in 1:Jcause){
      if (first){plot_pie_cell_first(e,Jcause,add=TRUE)}
      points_pie_cell(e,Jcause)
      first <- FALSE
    }
  }
  # cause names on the right edge:
  axis(4,at=1:Jcause,labels=paste(paste(cause_list_ord,ord,sep=" ("),")",sep=""),
       las=2,cex.axis=label_size)
  # cell bottom axis:
  abline(h=seq(1.5,Jcause-.5,by=1),lty=2,lwd=0.5,col="gray")
  
  mtext(expression(underline(hat(pi))),line=1,cex=1.8)
  legend("topright",c("prior","posterior"),lty=c(2,1),col=c("gray","purple"),
         lwd = 4,horiz=TRUE,cex=1,bty="n")
  
  
}
