#' Plot viral etiology posterior and distribution within bacterial and viral
#' etiologies
#'
#' DN: 1. the pathogens should have been ordered as in the perch_clean_data()
#' function. This function picks out pathogens not based on pathogen category lookup
#' table, but based on the order of pathogens that enter the analysis.
#'
#' @param DIR_NPLCM The file path to the result folder
#' @param ksFrac A number between 0 and 1, which is the fraction of samples used to
#' calculate kernel density
#' @param levellabel The contour line to be drawn in the final plot. Default is
#' \code{5}, which represents the contour of the \code{95} percent support region.
#' @importFrom coda read.coda
#' @importFrom robCompositions ternaryDiag
#' @import ks
#' @return A figure with group etiology.
#'
#' @export

nplcm_plot_group_etiology <- function(DIR_NPLCM,
                                      ksFrac = 1,
                                      levellabel = 5){
          #read in data from the folder directory:
          bugs.dat <- dget(paste(DIR_NPLCM,"data.txt",sep="/"))
          for (bugs.variable.name in names(bugs.dat)) {
            if (!is.null(dim(bugs.dat[[bugs.variable.name]]))) {
              dim(bugs.dat[[bugs.variable.name]]) <- rev(dim(bugs.dat[[bugs.variable.name]]))
              bugs.dat[[bugs.variable.name]] <- aperm(bugs.dat[[bugs.variable.name]])
            }
            assign(bugs.variable.name, bugs.dat[[bugs.variable.name]])
          }

          model_options  <- dget(paste(DIR_NPLCM,"model_options.txt",sep="/"))
          #compatibility checking:
          if (length(model_options$M_use)!=length(model_options$TPR_prior)){
            stop("The number of measurement source(s) is different from
                                 the number of TPR prior option!
                                 Make them equal, and match with order!")
          }


          # start plotting: ------------------------------------------------
            layout(matrix(c(1,1,2,3), 2, 2, byrow = TRUE),
                   widths=c(4,4), heights=c(5,3))
            #layout.show(n=3)
            cexval <- 1
            srtval <- 0

          #reading nplcm outputs:
          res_nplcm <- read.coda(paste(DIR_NPLCM,"coda1.txt",sep="/"),
                                 paste(DIR_NPLCM,"codaIndex.txt",sep="/"),
                                 quiet=TRUE)

          Jfull <- length(grep("pEti",colnames(res_nplcm)))

          # extract and process some data and posterior samples:
          SubVarName <- rep(NA,Jfull)
          for (j in 1:Jfull){
            SubVarName[j] = paste("pEti","[",j,"]",sep="")
          }


          #get etiology fraction MCMC samples:
          pEti_mat   <- res_nplcm[,SubVarName]
          pEti_mean  <- colMeans(pEti_mat)
          pEti_mean0 <- pEti_mean

          # order the pathogens by posterior mean:
          ord <- order(pEti_mean)
          if (is.null(model_options$SSonly) || model_options$SSonly==FALSE){
            pathogens_ord <- model_options$pathogen_list[ord]
          } else{
            pathogens_ord <- c(model_options$pathogen_list,
                               model_options$pathogen_SSonly_list)[ord]
          }
          pEti_mean_ord <- pEti_mean[ord]
          pEti_mat_ord  <- pEti_mat[,ord]

          JSS <- bugs.dat$JSS

            ##############################################
            ## plot 1: pr (virus is a cause)
            ##############################################
            op<-par()
            par(mai=c(0.5,3,1.09333,3))
            pr.v <- rowSums(pEti_mat[,ord[which(ord>JSS & ord <= bugs.dat$Jfull)]])

            plot(density(pr.v),type="l",xlim=c(0,1),xlab="",bty="n",
                 ylab="",yaxt="n",lwd=5,main="",cex.lab=2,cex.axis=2)

            alphaE <- bugs.dat$alpha
            if (is.null(model_options$SSonly)||model_options$SSonly==FALSE){
                  bacteria_index <- 1:JSS
            } else {
                  bacteria_index <- c(1:JSS,(1:bugs.dat$JSS_only)+bugs.dat$Jfull)
            }

            #tmp.rnd = rbeta(10000,sum(alphaE[-bacteria_index]),
            #                      sum(alphaE[bacteria_index]))

            tmp.x  = seq(0.001,0.999,by=0.001)
            tmp.y  = dbeta(tmp.x,sum(alphaE[-bacteria_index]),sum(alphaE[bacteria_index]))
            points(tmp.x,tmp.y,type="l",lty=2,col="grey",lwd=5)
            mtext("Probability of Viral Cause",side=3,line=0,cex=2)
            op


            #################################################
            ## plot 2: pr(B1, B2, others | cause is a Bacteria)
            ################################################
            par(mar = c(6.1,0,1,4.1),mai=c(.5,1,.5,.5),xpd=TRUE)

            if (is.null(model_options$SSonly)||model_options$SSonly==FALSE){

                  pEti0 = pEti_mat[,ord[which(ord < JSS+1)]]
            }else{
                  pEti0 = pEti_mat[,ord[which(ord < JSS+1 | ord > bugs.dat$Jfull )]]
            }

            pEti1  = data.frame(others=rowSums(pEti0[,-c(ncol(pEti0)-1,ncol(pEti0))]),
                                B2 = pEti0[,ncol(pEti0)-1],
                                B1 = pEti0[,ncol(pEti0)] )
            pEti = pEti1/rowSums(pEti1)

            if (is.null(model_options$SSonly)||model_options$SSonly==FALSE){

              names(pEti) = c("B-rest",
                              rev(pathogens_ord[which(ord<JSS+1)])[2],
                              rev(pathogens_ord[which(ord<JSS+1)])[1])
            }else{
              names(pEti) = c("B-rest",
                              rev(pathogens_ord[which(ord < JSS+1 | ord > bugs.dat$Jfull )])[2],
                              rev(pathogens_ord[which(ord < JSS+1 | ord > bugs.dat$Jfull )])[1])
            }

            #ksFrac = 1
            #levellabel = 5#c(5,20,50)
            postMean = colMeans(pEti)

            ternaryDiag(
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

            fhat1       = kde(x=temp,H=Hpi(temp),compute.cont=T)
            fhat2       = kde(x=temp,H=Hpi(temp),eval.points=ilr.grid)

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

            #################################################
            ## plot 3: pr(V1, V2, others | cause is a virus)
            ################################################
            op<-par()
            par(mar = c(6.1,4.1,1,2.1),mai=c(.5,.5,.5,1))
            pEti0 = pEti_mat[,ord[which(ord > JSS & ord <= bugs.dat$Jfull)]]
            pEti1  = data.frame(V2 = pEti0[,ncol(pEti0)-1],
                                others=rowSums(pEti0[,-c(ncol(pEti0)-1,ncol(pEti0))]),
                                V1 = pEti0[,ncol(pEti0)])
            pEti = pEti1/rowSums(pEti1)
            names(pEti) = c(rev(pathogens_ord[which(ord > JSS & ord <= bugs.dat$Jfull)])[2],
                            "V-rest",
                            rev(pathogens_ord[which(ord > JSS & ord <= bugs.dat$Jfull)])[1])

            #ksFrac = 1
            #levellabel = c(5,20,50)
            postMean = colMeans(pEti)
            ternaryDiag(
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

            fhat1       = kde(x=temp,H=Hpi(temp),compute.cont=T)
            fhat2       = kde(x=temp,H=Hpi(temp),eval.points=ilr.grid)

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
            op

}
