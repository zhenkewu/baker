#' Calculate marginal log odds ratios
#' 
#' @inheritParams simulate_nplcm
#' 
#' @examples
#'rm(list=ls())
#'library(cobs)
#'K.true  <- 2   # no. of latent subclasses in actual simulation. If eta = c(1,0), effectively, it is K.true=1
#'J       <- 5   # no. of pathogens.
#'N       <- 500 # no. of cases/controls.
#'
#'col_seq_cause <-  c("#DB9D85","#A2B367","#47BEA2","#70B3DA","#CD99D8")#colorspace::rainbow_hcl(5, start = 30, end = 300)
#'
#'subclass_mix_seq <- seq(0,1,by=0.05)
#'res      <- array(NA,c(J,J,length(subclass_mix_seq)))
#'res_cond <- array(NA,c(J,J,length(subclass_mix_seq),J))
#'
#'it <- layout(matrix(1:J^2,nrow=J,ncol=J,byrow=TRUE),
#'             heights = rep(3,J),
#'             widths  = rep(3,J)) 
#'
#'par(oma=c(8,10,8,3));  
#'
#'adj_seq <- c(0.15,0.5,0.85) # for roman numerals:
#'cex1       <- 2
#'cex_label1 <- 1
#'cex2       <- 2
#'cex_label2 <- 2
#'cex_margin_marks <- 2
#'
#'for (scn in c(1,2,3)){
#'  for (iter in seq_along(subclass_mix_seq)){
#'    curr_mix <- subclass_mix_seq[iter]
#'    lambda <- c(curr_mix,1-curr_mix)
#'    eta    <- c(curr_mix,1-curr_mix) 
#'    # if it is c(1,0),then it is conditional independence model, and
#'    # only the first column of parameters in PsiBS, ThetaBS matter!
#'    
#'    seed_start <- 20150923
#'    
#'    # set fixed simulation sequence:
#'    set.seed(seed_start)
#'    
#'    if (scn == 3){
#'      ThetaBS_withNA <- cbind(c(0.95,0.9,0.1,0.5,0.5),
#'                              c(0.95,0.1,0.9,0.5,0.5))
#'      PsiBS_withNA   <- cbind(c(0.4,0.4,0.05,0.2,0.2),
#'                              c(0.05,0.05,0.4,0.05,0.05))
#'    }
#'    
#'    if (scn == 2){
#'      ThetaBS_withNA <- cbind(c(0.95,0.5,0.5,0.5,0.5),
#'                              c(0.95,0.5,0.5,0.5,0.5))
#'      PsiBS_withNA   <- cbind(c(0.4,0.4,0.05,0.2,0.2),
#'                              c(0.05,0.05,0.4,0.05,0.05))
#'    }
#'    
#'    if (scn == 1){
#'      ThetaBS_withNA <- cbind(c(0.95,0.5,0.5,0.5,0.5),
#'                              c(0.95,0.5,0.5,0.5,0.5))
#'      PsiBS_withNA   <- cbind(c(0.3,0.3,0.15,0.2,0.2),
#'                              c(0.15,0.15,0.3,0.05,0.05))
#'    }
#'    
#'    # the following paramter names are set using names in the 'baker' package:
#'    set_parameter0 <- list(
#'      cause_list      = c(LETTERS[1:J]),
#'      etiology        = c(0.5,0.2,0.15,0.1,0.05), #same length as cause_list
#'      #etiology        = rep(0.2,J), #same length as cause_list
#'      pathogen_BrS    = LETTERS[1:J],
#'      meas_nm         = list(MBS = c("MBS1")),
#'      Lambda          = lambda,              #ctrl mix
#'      Eta             = t(replicate(J,eta)), #case mix, row number equal to Jcause.
#'      PsiBS           = PsiBS_withNA,
#'      ThetaBS         = ThetaBS_withNA,
#'      Nu      =     N, # control size.
#'      Nd      =     N  # case size.
#'    )
#'    
#'    res[,,iter] <- round(compute_logOR_single_cause(set_parameter0),2)
#'    
#'    for (pick in 1:J){
#'      set_parameter <- set_parameter0
#'      set_parameter$ThetaBS <- set_parameter0$PsiBS
#'      set_parameter$ThetaBS[pick,] <- set_parameter0$ThetaBS[pick,]
#'      set_parameter$etiology <- rep(0,J); set_parameter$etiology[pick] <- 1
#'      res_cond[,,iter,pick] <- round(compute_logOR_single_cause(set_parameter),2)
#'    }
#'  }
#'  
#'  ind <- sapply(c(0,0.5,1),function(x) which(subclass_mix_seq==x))
#'  logOR_lim <- c(-2.15,2.15)
#'  col_seq <- c("dodgerblue2","orange")
#'  logOR_seq <- log(c(0.25,0.5,1,2,4))
#'  for (j in 1:J){
#'    for (l in 1:J){
#'      
#'      par(mar=c(0,0,0,0)); 
#'      if (j==J){
#'        par(mar=c(0,0,0,0))
#'      }
#'      if (l%%J==0){
#'        par(mar=c(0,0,0,1)) 
#'      }
#'      if (l%%J==1){
#'        par(mar=c(0,1,0,0))
#'      }
#'      if (!(j==l)){
#'        plot(res[j,l,],type="l",xlab="",ylab="",
#'             ylim=logOR_lim, lwd=5,
#'             xaxt="n",
#'             yaxt="n",
#'             col=col_seq[1+(l>j)],
#'             #lty=c(2,1)[1+(l>j)],
#'             lty=1,
#'             bty="n"
#'        )
#'        box(col="lightgray")
#'        abline(h=0,col="lightgray",lwd=3,lty=3)
#'        
#'        if (j<l){
#'          matplot(res_cond[j,l,,],type="l",add=TRUE,pch=LETTERS[1:J],lwd=2,lty=2,
#'                  col=col_seq_cause)
#'        }
#'        lab_ord <- c(j,l); if (j>l){lab_ord <- rev(lab_ord)}
#'        mtext(paste0("(",set_parameter$pathogen_BrS[lab_ord[1]],",", 
#'                     set_parameter$pathogen_BrS[lab_ord[2]],")"), 
#'              side=3, adj=0.1,line=-2)
#'        
#'        if (l%%J==1){
#'          axis(2,at = logOR_seq, 
#'               labels = round(exp(logOR_seq),1),
#'               las=2,cex.axis=cex1)
#'        }
#'        
#'        if (l%%J==0){
#'          axis(4,at = logOR_seq, 
#'               labels = round(exp(logOR_seq),1),
#'               las=2,cex.axis=cex1)
#'        }
#'        
#'        if (j==J){
#'          axis(1,at=seq_along(subclass_mix_seq)[ind],labels=rep("",length(ind)),cex.axis = cex1,las=1)
#'          axis(1,at=seq_along(subclass_mix_seq)[ind]+c(1,rep(0,length(ind)-2),-1),labels=subclass_mix_seq[ind],cex.axis = cex1,las=1,tick=FALSE)
#'        }
#'        if (j==1){
#'          axis(3,at=seq_along(subclass_mix_seq)[ind],labels=rep("",length(ind)),cex.axis = cex1,las=1)
#'          axis(3,at=seq_along(subclass_mix_seq)[ind]+c(1,rep(0,length(ind)-2),-1),labels=subclass_mix_seq[ind],cex.axis = cex1,las=1,tick=FALSE)
#'        }
#'        if (j==5 & l==1){
#'          mtext(expression(atop("Odds Ratio","(log-scale)")), side = 2, line = 4, 
#'                cex=cex_label1, las=2)
#'        }
#'        if (j==5){
#'          mtext(expression(lambda[o]),side=1,line=4,cex=cex_label1)
#'        }
#'      }else{
#'        
#'        plot(1, type="n", axes=F, xlab="", ylab="", bty="n",
#'             xlim=c(0,1),ylim=c(0,1))
#'        
#'        
#'        if (j==3){
#'          text(labels=expression(CASES%up%""),x=.7,
#'               y=0.55,srt=-49,col=col_seq[2],cex=1.8,adj=0.5,font=4)
#'          text(labels=expression(CONTROLS%down%""),x=.42,
#'               y=0.38,srt=-49,col=col_seq[1],cex=1.8,adj=0.5,font=4)
#'        }
#'        if (j!=1 & j!=J){
#'          dg <- par("usr") 
#'          segments(dg[1],dg[4],dg[2],dg[3], col='lightgray',lwd=3)
#'        }
#'        if (j==J){
#'          legend("top",LETTERS[1:J],lty=2,col=col_seq_cause,cex = 1.5,lwd=2,
#'                 bty="n",horiz=FALSE)
#'        }
#'      }
#'    }
#'  }
#'}
#'
#'
#'
#' @return A figure showing pairwise odds ratios for cases (upper right, solid lines) 
#' and controls (lower left, broken lines) as the first subclass weight increases 
#' from 0 to 1. Pairwise independence is represented by the dotted horizontal lines 
#' for reference.
#' 
#' @export
#'
compute_logOR_single_cause <- function(set_parameter){
  eti <- set_parameter$etiology
  theta <- set_parameter$ThetaBS
  psi <- set_parameter$PsiBS
  lambda <- set_parameter$Lambda
  eta <- set_parameter$Eta
  K <- max(sum(lambda!=0), sum(eta[1,]!=0))
  k_ind_case <- 1:K
  k_ind_ctrl <- 1:K
  if (K==1){
    k_ind_case <- which(eta[1,]!=0)  
    k_ind_ctrl <- which(lambda!=0)  
  }
  J <- length(eti)
  MAT <- matrix(NA,nrow=J, ncol = J)
  for (j in 1:(J-1)){
    for (l in (j+1):J){
      A <- 0; B <- 0; C <- 0; D <- 0
      for (c in 1:J){
        ind_j <- 0+(c==j)
        ind_l <- 0+(c==l)
        A_subclass_temp <- 0
        B_subclass_temp <- 0
        C_subclass_temp <- 0
        D_subclass_temp <- 0
        for (k in k_ind_case){
          A_subclass_temp <- A_subclass_temp + eta[j,k]*theta[j,k]^ind_j*psi[j,k]^(1-ind_j)*(theta[l,k])^(ind_l)*(psi[l,k])^(1-ind_l)
          B_subclass_temp <- B_subclass_temp + eta[j,k]*(1-theta[j,k])^ind_j*(1-psi[j,k])^(1-ind_j)*(theta[l,k])^(ind_l)*(psi[l,k])^(1-ind_l)
          C_subclass_temp <- C_subclass_temp + eta[j,k]*(1-theta[j,k])^ind_j*(1-psi[j,k])^(1-ind_j)*(1-theta[l,k])^(ind_l)*(1-psi[l,k])^(1-ind_l)
          D_subclass_temp <- D_subclass_temp + eta[j,k]*theta[j,k]^ind_j*psi[j,k]^(1-ind_j)*(1-theta[l,k])^(ind_l)*(1-psi[l,k])^(1-ind_l)
        }
        A <- A + eti[c]*A_subclass_temp  
        B <- B + eti[c]*B_subclass_temp  
        C <- C + eti[c]*C_subclass_temp  
        D <- D + eti[c]*D_subclass_temp  
      }
      
      MAT[j,l] <- log(A)-log(B)+log(C)-log(D)
    }
  }
  
  for (l in 1:(J-1)){
    for (j in (l+1):J){
      A <- 0
      B <- 0
      C <- 0
      D <- 0
      for (k in k_ind_ctrl){
        A <- A + lambda[k]*psi[j,k]*psi[l,k]
        B <- B + lambda[k]*(1-psi[j,k])*(1-psi[l,k])
        C <- C + lambda[k]*psi[j,k]*(1-psi[l,k])
        D <- D + lambda[k]*(1-psi[j,k])*psi[l,k]
      }
      MAT[j,l] <- log(A)+log(B)-log(C)-log(D)
    }  
  }
  MAT
}
