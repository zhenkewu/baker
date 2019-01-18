## ----eval=TRUE,include=FALSE---------------------------------------------
#setwd("~/Desktop/baker_research")

## ----eval=T,message=FALSE------------------------------------------------
rm(list=ls())
library(baker)
library(lubridate)
library(rjags)
library(R2jags)

## ----eval=T--------------------------------------------------------------
J = 6                         # number of causes
N.SITE = 7                    # number of sites
cause_list = c(LETTERS[1:J])  # cause list
K = 2                         # number of subclasses
lambda = c(1,0)               # subclass weights for control group
eta = c(1,0)                  # subclass weights for case group
N = 300                       # number of subjects for control/case group


# etiology for all sites
etiology_allsites = list(c(0.5,0.2,0.15,0.05,0.05,0.05),
                     c(0.2,0.5,0.15,0.05,0.05,0.05),
                     c(0.2,0.15,0.5,0.05,0.05,0.05),
                     c(0.2,0.15,0.05,0.5,0.05,0.05),
                     c(0.2,0.15,0.05,0.05,0.5,0.05),
                     c(0.2,0.15,0.05,0.05,0.05,0.5),
                     c(0.05,0.2,0.15,0.5,0.05,0.05))
# FPR for MBS1 for all sites
PsiBS_MBS1_allsites = list(c(0.3,0.3,0.15,0.2,0.2,0.2),   
                           c(0.5,0.5,0.15,0.2,0.2,0.2),
                           c(0.2,0.2,0.15,0.2,0.2,0.2),
                           c(0.2,0.2,0.15,0.2,0.2,0.2),
                           c(0.1,0.1,0.15,0.2,0.2,0.2),
                           c(0.1,0.1,0.15,0.2,0.2,0.2),
                           c(0.1,0.1,0.15,0.2,0.2,0.2))

set.seed(20181018)
data_nplcm_list = lapply(1:N.SITE,function(siteID){
  set_parameter <- list(
    cause_list      = cause_list,
    etiology        = etiology_allsites[[siteID]], 
    pathogen_BrS    = LETTERS[1:J],
    SS = T,
    pathogen_SS     = LETTERS[1:2],
    meas_nm         = list(MBS = c("MBS1"),MSS=c("MSS1")),
    Lambda          = lambda,               # contral subclass weight for BrS   
    Eta             = t(replicate(J,eta)),  # case subclass weight for BrS
    PsiBS           = cbind(PsiBS_MBS1_allsites[[siteID]],PsiBS_MBS1_allsites[[siteID]]), # FPR
    PsiSS           = cbind(rep(0,J),rep(0,J)),
    ThetaBS         = cbind(c(0.95,0.9,0.85,0.9,0.9,0.9),     # TPR for MBS1
                            c(0.95,0.9,0.85,0.9,0.9,0.9)),
    ThetaSS         = cbind(c(0.25,0.10,0.15,0.05,0.15,0.15), # TPS for MSS1
                            c(0.25,0.10,0.15,0.05,0.15,0.15)),
    Nu      =     N, 
    Nd      =     N 
  )
  simu_out   <- simulate_nplcm(set_parameter)
  data_nplcm <- simu_out$data_nplcm
  data_nplcm$X = data.frame(SITE=rep(siteID,(set_parameter$Nu+set_parameter$Nu)))      # set X as a data frame containing SITE
  return(data_nplcm)
})

# put cases on top of controls
data_nplcm_order = unlist(lapply(1:N.SITE,function(i) ((i-1)*2*N+1):((i-1)*2*N+N)))
data_nplcm_order = c(data_nplcm_order, setdiff(1:(N*2*N.SITE),data_nplcm_order))
data_nplcm = list(Mobs=list(MBS = list(MBS1 = Reduce(rbind, lapply(data_nplcm_list, function(l) l$Mobs$MBS$MBS1))[data_nplcm_order,]),
                            MSS = list(MSS1 = Reduce(rbind, lapply(data_nplcm_list, function(l) l$Mobs$MSS$MSS1))[data_nplcm_order,]),
                            MGS=NULL),
                  Y = Reduce(c, lapply(data_nplcm_list, function(l) l$Y))[data_nplcm_order],
                  X = data.frame(SITE=Reduce(rbind, lapply(data_nplcm_list, function(l) l$X))[data_nplcm_order,]))

## ----eval=T--------------------------------------------------------------
BrS_object_1 <- make_meas_object(LETTERS[1:J],"MBS","1","BrS",cause_list)
SS_object_1 <- make_meas_object(LETTERS[1:2],"MSS","1","SS",cause_list)

model_options <- list(likelihood = list(cause_list = cause_list,              # <---- fitted causes.
                                        k_subclass = c(1),                    # <---- no. of subclasses.
                                        Eti_formula = ~ -1 + as.factor(SITE), # <---- etiology regression formula; only for cases.
                                        FPR_formula = list(
                                          MBS1 =  ~ -1+as.factor(SITE))),     
                      use_measurements = c("BrS","SS"),                       # <---- which measurements to use to inform etiology
                      prior = list(
                        Eti_prior = t(sapply(1:N.SITE, function(i) overall_uniform(1, cause_list))), # <--- etiology prior.
                        TPR_prior = list(
                          BrS = list(info  = "informative",
                                     input = "match_range",
                                     val   = list(
                                       MBS1 = list(up = list(rep(0.99,length(BrS_object_1$patho))),
                                                   low = list(rep(0.5,length(BrS_object_1$patho)))))),
                          SS = list(info  = "informative",
                                    input = "match_range",
                                    val   = list(
                                      MSS1 = list(up = list(rep(0.5,length(SS_object_1$patho))),
                                                  low = list(rep(0.01,length(SS_object_1$patho)))))) 
                          )  # <---- TPR prior.
                        )
                      )     

assign_model(model_options,data_nplcm)

## ----eval=T--------------------------------------------------------------
# parent directory for testing code (LOCAL):
working_dir <- tempdir()

# date stamp for analysis:
Date <- gsub("-", "_", Sys.Date())
# include stratification information in file name:
dated_strat_name <- file.path(working_dir,paste0(Date,"_discrete_predictor"))
# create folder
result_folder <- dated_strat_name
dir.create(result_folder)

# options for MCMC chains:
mcmc_options <- list(debugstatus = TRUE,
                     n.chains   = 1,
                     n.itermcmc = as.integer(50),
                     n.burnin   = as.integer(10),
                     n.thin     = 10,
                     individual.pred = !TRUE,
                     ppd             = !TRUE,
                     get.pEti        = !TRUE,
                     result.folder = result_folder,
                     bugsmodel.dir = result_folder,
                     jags.dir = "",
                     use_jags = TRUE)

# Record the settings of current analysis:
dput(data_nplcm,file.path(mcmc_options$result.folder,"data_nplcm.txt")) # <-- in case covariate data X are needed for plotting.

rjags::load.module("glm")

gs <- nplcm(data_nplcm,model_options,mcmc_options)


## ----eval=T,fig.width=7,fig.height=3-------------------------------------
# JAGS:
DIR_NPLCM <- result_folder

new_env   <- new.env()
source(file.path(DIR_NPLCM,"jagsdata.txt"),local=new_env)
bugs.dat <- as.list(new_env)
rm(new_env)

res_nplcm <- coda::read.coda(file.path(DIR_NPLCM,"CODAchain1.txt"),
                             file.path(DIR_NPLCM,"CODAindex.txt"),
                             quiet=TRUE)

print_res <- function(x) for (i in grep(x,colnames(res_nplcm))) plot(res_nplcm[,i],main=paste0('Trace of ',colnames(res_nplcm)[i]))
get_res   <- function(x) res_nplcm[,grep(x,colnames(res_nplcm))]

print_res("pEti")
print_res("psiBS")

## ----eval=TRUE-----------------------------------------------------------

# structure the posterior samples:
JBrS_1        <- ncol(bugs.dat$MBS_1) # number of pathogens in MBS1
n_samp_kept   <- nrow(res_nplcm)      # number of posterior sample after burn-in
Jcause        <- bugs.dat$Jcause      # number of causes
Nd            <- bugs.dat$Nd          # case size
Nu            <- bugs.dat$Nu          # control size
n_unique_Eti_level <- bugs.dat$n_unique_Eti_level  # number of stratums


# etiology:
plotid_Eti <- which(data_nplcm$Y==1) # <--- specifies who to look at.
Eti_prob_scale <- array(get_res("pEti"),c(n_samp_kept,n_unique_Eti_level,Jcause))

# posterior etiology mean for each cause for each site
Eti_mean <- apply(Eti_prob_scale,c(2,3),mean)   
# posterior etiology quantiles for each cause for each site
Eti_q    <- apply(Eti_prob_scale,c(2,3),quantile,c(0.025,0.975))

# marginalized posteior etiology ignoring site
Eti_overall <- apply(Eti_prob_scale,c(3,1),mean)
# posteior etiology mean for each cause across all sites
Eti_overall_mean <- rowMeans(Eti_overall)
# posteior etiology quantiles for each cause across all sites
Eti_overall_q    <- apply(Eti_overall,1,quantile,c(0.025,0.975))

## ----eval=TRUE,fig.width=25,fig.height=20--------------------------------
# weight to marginalize posterior etiology distributions
user_weight <- rep(1/N.SITE,N.SITE) # c(0.3,0.2,0.1,0.1,0.1,0.1,0.1)
# marginalized posterior etiology over all sites using user-defined weights
Eti_overall_usr_weight <- apply(Eti_prob_scale,1,function(S) t(S)%*%matrix(user_weight,ncol=1))
# marginalized posterior etiology mean using user-defined weights
Eti_overall_mean_usr_weight <- rowMeans(Eti_overall_usr_weight)
# marginalized posterior etiology quantiles using user-defined weights
Eti_overall_q_usr_weight    <- apply(Eti_overall_usr_weight,1,quantile,c(0.025,0.975))

# plot posterior distribution for etiology probability
par(mfcol=c(1+n_unique_Eti_level,Jcause),mar=c(3,15,1,0))
for (j in 1:Jcause){
  for (site in 1:n_unique_Eti_level){
    hist(Eti_prob_scale[,site,j],xlim=c(0,1),breaks="Scott",freq=FALSE,main="",xlab="")
    mtext(text = paste0('SITE',levels(as.factor(data_nplcm$X$SITE))[site],": ",cause_list[j]),3,-1,cex=2,adj = 0.9)
    if (j==1){
      mtext(paste0(round(user_weight[site],4)),2,5,cex=2,col="blue",las=1)
      if (site==5) {mtext("User-specified weight towards overall pie:", 2,12, cex=3)}
    }
  }
  hist(Eti_overall_usr_weight[j,],xlim=c(0,1),breaks="Scott",freq=FALSE,main="",col="blue",
       xlab="Etiology")
  mtext(text = paste0("Overall: ",cause_list[j]),3,adj=0.9,cex=2,col="blue")
}

