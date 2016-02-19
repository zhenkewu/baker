#' Plot grouped (two groups) posterior and cascade distributions (next top-three) 
#' within each group 
#'
#' DN: 1. the pathogens should have been ordered as in the perch_clean_data()
#' function. This function picks out pathogens not based on pathogen category lookup
#' table, but based on the order of pathogens that enter the analysis.
#'
#' @param DIR_NPLCM The file path to the result folder
#' @param dir_taxo File path to taxonomy information (.csv). The default is \code{NULL}.
#' If specified, it will overide \code{patho_taxo_dir} in \code{clean_options} read 
#' from \code{DIR_NPLCM}.
#' @param ksFrac A number between 0 and 1, which is the fraction of samples used to
#' calculate kernel density
#' @param levellabel The contour line to be drawn in the final plot. Default is
#' \code{5}, which represents the contour of the \code{95} percent support region.
#' @return A figure with group etiology.
#'
#' @family visualization functions
#' @export

plot_group_etiology <- function(DIR_NPLCM,dir_taxo=NULL,ksFrac = 1,levellabel = 5){

  # read NPLCM outputs:
  out           <- nplcm_read_folder(DIR_NPLCM)
  # organize ouputs:
  Mobs          <- out$Mobs
  Y             <- out$Y
  model_options <- out$model_options
  clean_options <- out$clean_options
  res_nplcm     <- out$res_nplcm
  bugs.dat      <- out$bugs.dat
  rm(out)
  
  Jcause <- length(model_options$likelihood$cause_list)
  
  template_BrS <- template_SS <- NULL
  check_combo_BrS <- check_combo_SS <- NULL
  
  if ("BrS" %in% model_options$use_measurements){
    template_BrS <- lapply(clean_options$BrS_objects,"[[","template")
    check_combo_BrS <- any(unlist(lapply(template_BrS,rowSums))>1)
    }
  if ("SS" %in% model_options$use_measurements){
    template_SS  <- lapply(clean_options$SS_objects,"[[","template")
    check_combo_SS <- any(unlist(lapply(template_SS,rowSums))>1)
    }
  if (is.null(template_BrS) && is.null(template_SS)){stop("== No BrS or SS used in the fit. ==")}
  
  
  if ((!is.null(check_combo_BrS) && check_combo_BrS) | 
      (!is.null(check_combo_SS) && check_combo_SS) ){
    stop("== Grouped visualziation not implemented for combo latent status. ==")
  }
  
  if (is.null(dir_taxo)){
   taxo_cause_list <- assign_taxo_cause_list(model_options$likelihood$cause_list,
                                            clean_options$patho_taxo_dir)
  } else{
    taxo_cause_list <- assign_taxo_cause_list(model_options$likelihood$cause_list,
                                              dir_taxo)
  }
  
  bacteria_index <- which(taxo_cause_list=="B")
  virus_index <- which(taxo_cause_list=="V")
  
  
  
  old_par <- par(no.readonly=TRUE)
  on.exit(par(old_par))
  
  
  # start plotting: ------------------------------------------------
  
  layout(matrix(c(1,1,2,3), 2, 2, byrow = TRUE),
           widths=c(4,4), heights=c(5,3))
  #layout.show(n=3)
  cexval <- 1
  srtval <- 0
  
  # extract and process some data and posterior samples:
  SubVarName <- rep(NA,Jcause)
  for (j in 1:Jcause){
    SubVarName[j] = paste("pEti","[",j,"]",sep="")
  }
  
  #get etiology fraction MCMC samples:
  pEti_mat   <- res_nplcm[,SubVarName]
  pEti_mean  <- colMeans(pEti_mat)
  pEti_mean0 <- pEti_mean
  
  # order the pathogens by posterior mean:
  ord <- order(pEti_mean)
  cause_list_ord <- model_options$likelihood$cause_list[ord]

  pEti_mean_ord <- pEti_mean[ord]
  pEti_mat_ord  <- pEti_mat[,ord]
  

  
  ##############################################
  ## plot 1: pr (virus is a cause)
  ##############################################
  #op<-par()
  par(mai=c(0.5,3,1.09333,3))
  pr.v <- rowSums(pEti_mat[,ord[taxo_cause_list=="V"]])
  
  plot(density(pr.v),type="l",xlim=c(0,1),xlab="",bty="n",
       ylab="",yaxt="n",lwd=5,main="",cex.lab=2,cex.axis=2)
  
  alphaE <- bugs.dat$alpha

  
  #tmp.rnd = rbeta(10000,sum(alphaE[-bacteria_index]),
  #                      sum(alphaE[bacteria_index]))
  
  tmp.x  = seq(0.001,0.999,by=0.001)
  tmp.y  = dbeta(tmp.x,sum(alphaE[-bacteria_index]),sum(alphaE[bacteria_index]))
  points(tmp.x,tmp.y,type="l",lty=2,col="grey",lwd=5)
  mtext("Probability of Viral Cause",side=3,line=0,cex=2)
  #op
  
 if (length(bacteria_index)>=3){
      #################################################
      ## plot 2: pr(B1, B2, others | cause is a Bacteria)
      ################################################
      par(mar = c(6.1,0,1,4.1),mai=c(.5,1,.5,.5),xpd=TRUE)
      
      pEti0 = pEti_mat[,ord[which(ord %in% bacteria_index)]]
      
      pEti1  = data.frame(others=rowSums(pEti0[,-c(ncol(pEti0)-1,ncol(pEti0)),drop=FALSE]),
                          B2 = pEti0[,ncol(pEti0)-1],
                          B1 = pEti0[,ncol(pEti0)] )
      pEti = pEti1/rowSums(pEti1)
      
      names(pEti) = c("B-rest",
                        rev(cause_list_ord[which(ord %in% bacteria_index )])[2],
                        rev(cause_list_ord[which(ord %in% bacteria_index )])[1])
      
      #ksFrac = 1
      #levellabel = 5#c(5,20,50)
      postMean = colMeans(pEti)
      
      robCompositions::ternaryDiag(
        as.data.frame(matrix(postMean,nrow=1,ncol=3)),
        pch = "",
        col = "black",
        cex=3,name=rep(NA,3),
        lwd=3)
      
      ind.temp = floor(seq(1,nrow(pEti),len=floor(ksFrac*nrow(pEti))))
      temp = robCompositions::isomLR(pEti)[ind.temp,]
      
      # a 2 dim grid with in [0,1]*[0,1]:
      coord1 = seq(0.001,0.999,by=0.01)
      coord2 = sqrt(0.75)*coord1
      coord.grid = expand.grid(coord1=coord1,coord2=coord2)
      
      # transform Euclidean distance to probabilities, there will
      # be several rows with negatives because the (coord1, coord2) do not correspond
      # to a probabiliy vector.
      prob.grid.temp  = data.frame(z=coord.grid$coord2/sqrt(0.75),
                                   y=coord.grid$coord1-1/sqrt(3)*coord.grid$coord2,
                                   x=1-coord.grid$coord2/sqrt(0.75)-coord.grid$coord1+
                                     1/sqrt(3)*coord.grid$coord2)
      prob.grid = rev(prob.grid.temp)
      
      ilr.grid    = suppressWarnings(robCompositions::isomLR(prob.grid))
      
      fhat1       = ks::kde(x=temp,H=ks::Hpi(temp),compute.cont=T)
      fhat2       = ks::kde(x=temp,H=ks::Hpi(temp),eval.points=ilr.grid)
      
      tri.z       = matrix(fhat2$estimate,nrow=length(coord1),ncol=length(coord2))
      
      f = function(x,y){
        (y<=sqrt(3)*x & -y/sqrt(3)>=(x-1))
      }
      
      for (i in 1:length(coord1)){
        for (j in 1:length(coord2)){
          tri.z[i,j]<-ifelse(f(coord1[i],coord2[j]),tri.z[i,j],NA)
        }
      }
      
      
      num.levels      = levellabel#c(5,10,20,50)
      cont.levels     = paste(num.levels,"%",sep="")
      label.levels     = paste(100-num.levels,"%",sep="")
      
      points(1/2*(postMean[3]+2*postMean[2]),sqrt(3)/2*postMean[3],
             cex=2,col="red",pch=15)
      contour(coord1,coord2,tri.z,
              levels =(fhat1$cont)[cont.levels],
              labels = label.levels,add=TRUE,lwd=4,col=c("blue"),lty=1)
      mcex=2
      name = names(pEti)
      mtext(name[1], side = 1, line = 1, at = -0.1, cex = mcex)
      mtext(name[2], side = 1, line = 1, at = 1.1, cex = mcex)
      mtext(name[3], side = 3,  line= -1, at=0.5 ,cex = mcex)
  }
  #################################################
  ## plot 3: pr(V1, V2, others | cause is a virus)
  ################################################
  #op<-par()
  
  if (length(virus_index)>=3){
     if (length(bacteria_index)<3){frame()}
      par(mar = c(6.1,4.1,1,2.1),mai=c(.5,.5,.5,1))
      pEti0 = pEti_mat[,ord[which(ord %in% virus_index)]]
      pEti1  = data.frame(V2 = pEti0[,ncol(pEti0)-1],
                          others=rowSums(pEti0[,-c(ncol(pEti0)-1,ncol(pEti0)),drop=FALSE]),
                          V1 = pEti0[,ncol(pEti0)])
      pEti = pEti1/rowSums(pEti1)
      names(pEti) = c(rev(cause_list_ord[which(ord %in% virus_index)])[2],
                      "V-rest",
                      rev(cause_list_ord[which(ord %in% virus_index)])[1])
      
      #ksFrac = 1
      #levellabel = c(5,20,50)
      postMean = colMeans(pEti)
      robCompositions::ternaryDiag(
        as.data.frame(matrix(postMean,nrow=1,ncol=3)),
        pch = "",
        col = "black",
        cex=3,name=rep(NA,3),
        lwd=3)
      
      ind.temp = floor(seq(1,nrow(pEti),len=floor(ksFrac*nrow(pEti))))
      #temp = compositions::ilr(pEti)[ind.temp,]
      temp = robCompositions::isomLR(pEti)[ind.temp,]
      
      coord1 = seq(0.001,0.999,by=0.01)
      coord2 = sqrt(0.75)*coord1
      
      coord.grid = expand.grid(coord1=coord1,coord2=coord2)
      prob.grid.temp  = data.frame(z=coord.grid$coord2/sqrt(0.75),
                                   y=coord.grid$coord1-1/sqrt(3)*coord.grid$coord2,
                                   x=1-coord.grid$coord2/sqrt(0.75)-coord.grid$coord1+
                                     1/sqrt(3)*coord.grid$coord2)
      prob.grid = rev(prob.grid.temp)
      
      ilr.grid    = suppressWarnings(robCompositions::isomLR(prob.grid))
      
      fhat1       = ks::kde(x=temp,H=ks::Hpi(temp),compute.cont=T)
      fhat2       = ks::kde(x=temp,H=ks::Hpi(temp),eval.points=ilr.grid)
      
      tri.z       = matrix(fhat2$estimate,nrow=length(coord1),ncol=length(coord2))
      
      f = function(x,y){
        (y<=sqrt(3)*x & -y/sqrt(3)>=(x-1))
      }
      
      for (i in 1:length(coord1)){
        for (j in 1:length(coord2)){
          tri.z[i,j]<-ifelse(f(coord1[i],coord2[j]),tri.z[i,j],NA)
        }
      }
      
      num.levels      = levellabel
      cont.levels     = paste(num.levels,"%",sep="")
      label.levels     = paste(100-num.levels,"%",sep="")
      
      points(1/2*(postMean[3]+2*postMean[2]),sqrt(3)/2*postMean[3],
             cex=2,col="red",pch=15)
      contour(coord1,coord2,tri.z,
              levels =(fhat1$cont)[cont.levels],
              labels = label.levels,add=TRUE,lwd=4,col=c("blue"),lty=1)
      mcex=2
      name = names(pEti)
      mtext(name[1], side = 1, line = 1, at = -0.1, cex = mcex)
      mtext(name[2], side = 1, line = 1, at = 1.1, cex = mcex)
      mtext(name[3], side = 3,  line= -1, at=0.5 ,cex = mcex)
      #op
    }
  }



#' Plot grouped (two groups) posterior and cascade distributions (next top-three) 
#' within each group
#'
#' @details By default the function automatically orders pathogens by posterior 
#' means (increase clockwise)
#' 
#' @param selected a vector of character strings of length 3: the names of etiologies
#' to visualize.
#' @param DIR_NPLCM The file path to the result folder
#' @param ksFrac A number between 0 and 1, which is the fraction of samples used to
#' calculate kernel density
#' @param levellabel The contour line to be drawn in the final plot. Default is
#' \code{5}, which represents the contour of the \code{95} percent support region.
#' @return A figure with group etiology.
#' @family visualization functions
#' @export

plot_selected_etiology <- function(selected, DIR_NPLCM,ksFrac = 1,levellabel = 5){
  
  
  # read NPLCM outputs:
  out           <- nplcm_read_folder(DIR_NPLCM)
  # organize ouputs:
  Mobs          <- out$Mobs
  Y             <- out$Y
  model_options <- out$model_options
  clean_options <- out$clean_options
  res_nplcm     <- out$res_nplcm
  bugs.dat      <- out$bugs.dat
  rm(out)
  
  Jcause <- length(model_options$likelihood$cause_list)
  
  template_BrS <- template_SS <- NULL
  check_combo_BrS <- check_combo_SS <- NULL
  
  if ("BrS" %in% model_options$use_measurements){
    template_BrS <- lapply(clean_options$BrS_objects,"[[","template")
    check_combo_BrS <- any(unlist(lapply(template_BrS,rowSums))>1)
  }
  if ("SS" %in% model_options$use_measurements){
    template_SS  <- lapply(clean_options$SS_objects,"[[","template")
    check_combo_SS <- any(unlist(lapply(template_SS,rowSums))>1)
  }
  if (is.null(template_BrS) && is.null(template_SS)){stop("== No BrS or SS used in the fit. ==")}
  
  
  if ((!is.null(check_combo_BrS) && check_combo_BrS) | 
      (!is.null(check_combo_SS) && check_combo_SS) ){
    stop("== Grouped visualziation not implemented for combo latent status. ==")
  }
  
  match_index <- match(selected, model_options$likelihood$cause_list)
  if (sum(is.na(match_index))>0){stop("== There is a selected name not in the fitted model (cause list). ==")}
  
  vis_index   <- match_index
  
  old_par <- par(no.readonly=TRUE)
  on.exit(par(old_par))
  
  
  # start plotting: ------------------------------------------------
  layout(matrix(c(1), 1, 1, byrow = TRUE),
         widths=c(6), heights=c(6))
  #layout.show(n=3)
  cexval <- 1
  srtval <- 0
  
  # extract and process some data and posterior samples:
  SubVarName <- rep(NA,Jcause)
  for (j in 1:Jcause){
    SubVarName[j] = paste("pEti","[",j,"]",sep="")
  }
  
  #get etiology fraction MCMC samples:
  pEti_mat   <- res_nplcm[,SubVarName]
  pEti_mean  <- colMeans(pEti_mat)
  pEti_mean0 <- pEti_mean
  
  # order the pathogens by posterior mean:
  ord <- order(pEti_mean)
  cause_list_ord <- model_options$likelihood$cause_list[ord]
  
  pEti_mean_ord <- pEti_mean[ord]
  pEti_mat_ord  <- pEti_mat[,ord]
  
  
  
  #################################################
  ## plot any three etiologies:
  ################################################
  par(mar = c(6.1,0,1,4.1),mai=c(.5,1,.5,.5),oma=c(5,5,5,5),xpd=TRUE)
  
  
  pEti0 = pEti_mat[,ord[which(ord %in% vis_index)]]
  
  pEti1  = data.frame(others=rowSums(pEti0[,-c(ncol(pEti0)-1,ncol(pEti0)),drop=FALSE]),
                      B2 = pEti0[,ncol(pEti0)-1],
                      B1 = pEti0[,ncol(pEti0)] )
  pEti = pEti1/rowSums(pEti1)
  
  names(pEti) = c(rev(cause_list_ord[which(ord %in% vis_index )])[3],
                  rev(cause_list_ord[which(ord %in% vis_index )])[2],
                  rev(cause_list_ord[which(ord %in% vis_index )])[1])
  
  #ksFrac = 1
  #levellabel = 5#c(5,20,50)
  postMean = colMeans(pEti)
  
  robCompositions::ternaryDiag(
    as.data.frame(matrix(postMean,nrow=1,ncol=3)),
    pch = "",
    col = "black",
    cex=3,name=rep(NA,3),
    lwd=3)
  
  ind.temp = floor(seq(1,nrow(pEti),len=floor(ksFrac*nrow(pEti))))
  temp = robCompositions::isomLR(pEti)[ind.temp,]
  
  # a 2 dim grid with in [0,1]*[0,1]:
  coord1 = seq(0.001,0.999,by=0.01)
  coord2 = sqrt(0.75)*coord1
  coord.grid = expand.grid(coord1=coord1,coord2=coord2)
  
  # transform Euclidean distance to probabilities, there will
  # be several rows with negatives because the (coord1, coord2) do not correspond
  # to a probabiliy vector.
  prob.grid.temp  = data.frame(z=coord.grid$coord2/sqrt(0.75),
                               y=coord.grid$coord1-1/sqrt(3)*coord.grid$coord2,
                               x=1-coord.grid$coord2/sqrt(0.75)-coord.grid$coord1+
                                 1/sqrt(3)*coord.grid$coord2)
  prob.grid = rev(prob.grid.temp)
  
  ilr.grid    = suppressWarnings(robCompositions::isomLR(prob.grid))
  
  fhat1       = ks::kde(x=temp,H=ks::Hpi(temp),compute.cont=T)
  fhat2       = ks::kde(x=temp,H=ks::Hpi(temp),eval.points=ilr.grid)
  
  tri.z       = matrix(fhat2$estimate,nrow=length(coord1),ncol=length(coord2))
  
  f = function(x,y){
    (y<=sqrt(3)*x & -y/sqrt(3)>=(x-1))
  }
  
  for (i in 1:length(coord1)){
    for (j in 1:length(coord2)){
      tri.z[i,j]<-ifelse(f(coord1[i],coord2[j]),tri.z[i,j],NA)
    }
  }
  
  
  num.levels      = levellabel#c(5,10,20,50)
  cont.levels     = paste(num.levels,"%",sep="")
  label.levels     = paste(100-num.levels,"%",sep="")
  
  points(1/2*(postMean[3]+2*postMean[2]),sqrt(3)/2*postMean[3],
         cex=2,col="red",pch=15)
  contour(coord1,coord2,tri.z,
          levels =(fhat1$cont)[cont.levels],
          labels = label.levels,add=TRUE,lwd=4,col=c("blue"),lty=1)
  mcex=2
  name = names(pEti)
  mtext(name[1], side = 1, line = 1, at = -0.1, cex = mcex)
  mtext(name[2], side = 1, line = 1, at = 1.1, cex = mcex)
  mtext(name[3], side = 3,  line= -1, at=0.5 ,cex = mcex)
  
}
