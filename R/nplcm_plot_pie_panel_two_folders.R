#' Plot etiology (pie) panel
#' 
#' Now only works for singleton etiologies.
#' 
#' @param model_options_list See \code{\link{nplcm}}
#' @param res_nplcm_list A list of posterior samples
#' @param bugs.dat_list Data input for the model fitting.
#' @param index_for_order The index of the list to be used for 
#' ordering the pathogens. Because each model output might give
#' different ordering of pathogens (based on posterior mean).
#' @param dir_nplcm_list A list of file paths to the folders containing 
#' posterior samples
#' @param top_pie Numerical value to specify the rightmost limit 
#' on the horizontal axis for the pie panel.
#' 
#' @importFrom binom binom.confint
#' 
#' @export

nplcm_plot_pie_panel_two_folders <- function(model_options_list,
                                             res_nplcm_list,
                                             bugs.dat_list,
                                             index_for_order,
                                             dir_nplcm_list,
                                             top_pie = 1){
  #
  # now only deal with singleton etiologies:
  # 
  
  # total no. of causes:
  Jcause     <- length(model_options_list[[index_for_order]]$cause_list)
  # extract and process some data and posterior samples:
  SubVarName <- rep(NA,Jcause)
  for (j in 1:Jcause){
    SubVarName[j] = paste("pEti","[",j,"]",sep="")
  }
  
  # get etiology fraction MCMC samples:
  pEti_mat   <- res_nplcm_list[[index_for_order]][,SubVarName]
  pEti_mean  <- colMeans(pEti_mat)
  # order the causes by posterior mean:
  ord <- order(pEti_mean)
  
  prior_flag <- TRUE
  dotcolor_vec <- c("black","dodgerblue2")
  
  seq_folder_iter <- c(index_for_order,seq_along(model_options_list)[-index_for_order])
  for (ind_plot in 1:length(model_options_list)){  
    
    ind_folder <- seq_folder_iter[ind_plot]
    pEti_mat   <- res_nplcm_list[[ind_folder]][,SubVarName]
    pEti_mean  <- colMeans(pEti_mat)
    pEti_mean_ord <- pEti_mean[ord]
    pEti_mat_ord  <- pEti_mat[,ord]
    
    
    model_options <- model_options_list[[ind_folder]]
    
    
    # quantiles for etiology: outer is 97.5% CI, inner is 50% CI
    pEti_q1   <- apply(pEti_mat,2,quantile,probs=0.025)[ord]
    pEti_q2   <- apply(pEti_mat,2,quantile,probs=0.975)[ord]
    pEti_innerq1   <- apply(pEti_mat,2,quantile,probs=0.25)[ord]
    pEti_innerq2   <- apply(pEti_mat,2,quantile,probs=0.75)[ord]
    
    # complete list of pathogens:
    if (is.null(model_options$SSonly) || model_options$SSonly==FALSE){
      pathogen_list     <- model_options$pathogen_BrS_list
    } else{
      pathogen_list     <- c(model_options$pathogen_BrS_list,
                             model_options$pathogen_SSonly_list)
    }
    
    Jfull               <- length(pathogen_list)
    
    if (Jfull != Jcause){
      stop("== The number of causes is different from the total number of pathogens.
           The multiple-cause visualization, including NoA, is being developed.
           Please contact developer. Thanks. ==")
    }
    JBrS                <- length(model_options$pathogen_BrS_list)
    pathogens_ord       <- pathogen_list[ord]
    
    bugs.dat <- bugs.dat_list[[ind_folder]]
    Nd <- bugs.dat$Nd
    Nu <- bugs.dat$Nu
    
    #
    # plot etiology information:
    #
    dotcolor <- dotcolor_vec[ind_folder]
    
    op <- par(mar=c(5.1,0,4.1,10))
    if(prior_flag){
      plot(c(pEti_mean_ord),1:(Jfull)+(ind_plot-1)*0.2-0.25,
           yaxt="n",#xaxt="n",
           xlim=c(0,top_pie),ylim=c(0.5,Jfull+0.5),col=dotcolor,
           ylab="",xlab="probability",
           pch=c(20),cex=2)
    }else{
      points(c(pEti_mean_ord),1:(Jfull)+(ind_plot-1)*0.2-.25,
             yaxt="n",#xaxt="n",
             xlim=c(0,top_pie),ylim=c(0.5,Jfull+0.5),col=dotcolor,
             ylab="",xlab="probability",
             pch=c(20),cex=2)
    }
    
    axis(4,at=1:Jfull,labels=paste(paste(pathogens_ord,ord,sep=" ("),")",sep=""),las=2,cex.axis=1.5)
    abline(h=seq(1.5,Jfull-.5,by=1),lty=2,lwd=0.5,col="blue")
    #draw axis within plot:
    for (s in 1:(Jfull-1)){
      axis(1, seq(0,1,by=.2), lwd=0,lwd.ticks=1,#labels=rep("",length(seq(0,1,by=.2))),
           pos = seq(1.5,Jfull-.5,by=1)[s], cex.axis = 0.8,lty=2,col="blue")
    }
    points(c(pEti_q1),1:(Jfull)+(ind_plot-1)*0.2-0.25,pch="|",cex=1,col=dotcolor)
    points(c(pEti_q2),1:(Jfull)+(ind_plot-1)*0.2-0.25,pch="|",cex=1,col=dotcolor)
    points(c(pEti_innerq1),1:(Jfull)+(ind_plot-1)*0.2-0.25,pch="[",cex=1,col=dotcolor)
    points(c(pEti_innerq2),1:(Jfull)+(ind_plot-1)*0.2-0.25,pch="]",cex=1,col=dotcolor)
    
    mtext(expression(underline(hat(pi))),line=1,cex=1.8)
    
    #       legend("topright",c("prior","posterior"),lty=c(2,1),col=c("gray","black"),
    #              lwd = 4,horiz=TRUE,cex=1,bty="n")
    #       pgrid = seq(0,1,by=0.01)
    #       
    #       alpha <- bugs.dat$alpha#eti_prior_set(model_options)
    #       
    for (s in 1:(Jfull)){
      segments(y0=s+(ind_plot-1)*0.2-0.25,x0=c(pEti_q1)[s],y1=s+(ind_plot-1)*0.2-0.25,x1=c(pEti_q2)[s],col=dotcolor)
      segments(y0=s+(ind_plot-1)*0.2-0.25,x0=c(pEti_innerq1)[s],y1=s+(ind_plot-1)*0.2-0.25,x1=c(pEti_innerq2)[s],col = dotcolor,lwd=2)
      text(.9,s+(ind_plot-1)*0.3-0.25,paste0("=",paste0(round(100*c(pEti_mean_ord),1)[s],"%")),srt=0,cex=2,col=dotcolor)
      text(.75,s+(ind_plot-1)*0.3-0.25,bquote(hat(pi)[.(ord[s])]),srt=0,cex=2,col=dotcolor)
      #         if (prior_flag){
      #             tmp.density = dbeta(pgrid,alpha[ord[s]],sum(alpha[-ord[s]]))
      #             points(pgrid,tmp.density/(3*max(tmp.density))+s-0.45,type="l",col="gray",lwd=4,lty=2)
      #         }
      #         ##posterior density
      #         tmp.post.density = density(pEti_mat_ord[,s],from=0,to=1)
      #         tmp.x = tmp.post.density$x
      #         tmp.y = tmp.post.density$y
      #         points(tmp.x,tmp.y/(3*max(tmp.y))+s-0.45,col=dotcolor,lwd=4,type="l")
      #         
    }

    
    prior_flag <- ifelse(prior_flag,FALSE,TRUE)
    }
  
  legend("topleft",# inset=c(0.1,0.12),
         legend=c(dir_nplcm_list[[1]],dir_nplcm_list[[2]])[seq_folder_iter], 
           lty=c(1,1),pch=c(20,20), 
         title="Folders",col=dotcolor_vec[seq_folder_iter],bty='n')
  
}
  