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
#' 
#' @importFrom binom binom.confint
#' @family visualization functions
#' @export

plot_pie_panel <- function(model_options,
                           res_nplcm,
                           bugs.dat,
                           bg_color,
                           select_latent = NULL,
                           exact = TRUE,
                           top_pie = 1,
                           label_size = 1,
                           ref_eti = NULL){
  
  # order cause_list by posterior means:
  ord <- order_post_eti(res_nplcm,model_options)$ord
  pEti_mat_ord <- order_post_eti(res_nplcm,model_options)$pEti_mat_ord
  pEti_mean_ord <- colMeans(pEti_mat_ord)
  
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
  pEti_q         <- pEti_q[,latent_seq,drop=FALSE]
  cause_list_ord <- cause_list_ord[latent_seq]
  
  
  Jcause <- length(model_options$likelihood$cause_list)
  alpha_ord <- bugs.dat$alphaEti[ord]
  
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
  
  cat("\n == Plotting pies == \n")
  
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
