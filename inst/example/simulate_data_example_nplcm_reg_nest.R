### The functions below specify the CSCF and subclass weight generating processes for simulating this data. 


# 1. cause-specific case fractions, or "etiology" - to weight TPRs:
compute_pEti <- function(X, L){# output N by L:
  # a few functions of a continuous covariate:
  f1 <- function(x) {1*sin(8*pi*(x-0.5)/7)}
  f2 <- function(x) {-0.5*sin(2*pi*x)/(pi*x)}
  f3 <- function(x) {-0.5*sin(8*pi*(x-0.5)/7)}
  f4 <- function(x) {0}
  f5 <- function(x) {4*(exp(3*x)/(1+exp(3*x))-0.5)}
  # a list of 8 functions, each of which is used as the true linear predictor 
  # in etiology regression for log(pi_k/pi_K):
  f_list <- make_list(f1,f5,f4,f4,f4,f4,f4,f4)
  
  K_here <- L              # number of causes
  x      <- scale(X$date)  # standardize the dates on which subject were observed.
  N      <- nrow(X)        # number of all subjects.
  Mu     <- matrix(NA,nrow=N,ncol=K_here-1) # linear predictors for all subjects.
  for (k in 1:(K_here-1)){  # only need K_here-1 for identifiability.
    Mu[,k] <- f_list[[k]](x)
  }
  
  Alpha      <- cbind(Mu) #+ matrix(rnorm(N*K.true,0,1),nrow=N,ncol=K.true-1)
  t(apply(cbind(expit(Alpha),1),1,tsb))
}

# 2. subclass weights among cases:
compute_eta  <- function(X, K){# output Nd by K:
  f1 <- function(x) {4*(exp(3*x)/(1+exp(3*x))-0.5)}
  f2 <- function(x) {0}
  f_list <- make_list(f1,f2)
  
  x  <- X[,1]
  N  <- nrow(X)
  Mu <- matrix(NA,nrow=N,ncol=K-1) # K is the number of subclasses, which may differ from K_here.
  for (k in 1:(K-1)){Mu[,k] <- f_list[[k]](x)}
  
  Alpha      <- cbind(Mu) 
  t(apply(cbind(expit(Alpha),1),1,tsb)) # `tsb` calculates the probabilities from
  # parameters in a truncated stick-breaking formulation.
}

# 3. subclass weights among controls:
compute_nu   <- function(X, K){ # output Nu by K:
  f1_rev <- function(x) {4*(exp(-3*x)/(1+exp(-3*x))-0.5)} # this is f1's mirror image around x=0.
  f2     <- function(x) {0}
  f_list <- make_list(f1_rev,f2)
  
  x  <- X[,1]
  N  <- nrow(X)
  Mu <- matrix(NA,nrow=N,ncol=K-1)
  for (k in 1:(K-1)){Mu[,k] <- f_list[[k]](x)}
  
  Alpha      <- cbind(Mu) 
  t(apply(cbind(expit(Alpha),1),1,tsb)) # `tsb` calculates the probabilities from
  # parameters in a truncated stick-breaking formulation.
}


### FINISHED SETTING UP REGRESSION FUNCTIONS.


#### The following simulate data, again by using `simulate_nplcm()` in a different way.
#### We iterate over each subject when simulating data. Because `simulate_nplcm()` simulates
#### a pair of data (case-control by design), we can simulate two hypothetical subjects'
#### data. For an actual case subject, we discard the hypothetical control's simulated data; for
#### an actual control subject, we discard the hypothetical case's simulated data. This
#### is implemeneted as follows.


# create covariate dataset (site ID and enrollment date)
X   <- data.frame(siteID = rep(1:2,times=c(Nd+Nu,Nd+Nu)),
                  date   = sample(1:300,sum(c(Nd+Nu,Nd+Nu)),replace=TRUE)) # date ids.

# create vector of cases and controls 
Y <- rep(c(1,0,1,0),times=c(Nd,Nu,Nd,Nu))

# standardize the dates to simplify the regression fitting 
X$std_date <- dm_Rdate_FPR(X$date,Y,effect = "fixed") # standardized dates from the date ids.

# simulate the etiology probabilities
simu_pEti  <- compute_pEti(X, L.CAUSE)

## simulate the weights for cases and controls
simu_nu      <- compute_nu(X$std_date,K)
simu_eta     <- compute_eta(X$std_date,K)  
simu_PR_ctrl <- simu_nu%*%t(PsiBS0) # marginal positive rates, averaged over control subclasses.

out_list     <- vector("list",length=N)
for (i in 1:N){ # a total sample size of 2N.
  # setup parameters for the present individual:
  set_parameter <- list(
    cause_list      = cause_list,    # the vector of characters of causes
    etiology        = simu_pEti[i,], # etiologies; only meaningful for cases 
    pathogen_BrS    = LETTERS[1:J.BrS], # names of the BrS measurements
    pathogen_SS     = LETTERS[1:J.SS], # names of the SS measurements
    meas_nm         = list(MBS = c("MBS1"),MSS=c("MSS1")),
    Lambda          = simu_nu[i,],    # control subclass weight for BrS   
    Eta             = t(replicate(J.BrS,simu_eta[i,])),  # case subclass weight for BrS
    PsiBS           = cbind(rep(0.5,J.BrS),rep(0.05,J.BrS)), # FPR for BrS
    PsiSS           = cbind(rep(0,J.SS)),                     # FPR for SS measures; perfect specificity.
    ThetaBS         = cbind(c(0.95,0.90,0.95,0.95,0.95,0.99),    # TPR for BrS
                            c(0.95,0.90,0.95,0.85,0.95,0.95)
    ),
    ThetaSS         = cbind(rep(0.2,J.SS)),                   # TPR for SS
    Nd      =     1,  # the current individual.
    Nu      =     1 
  )
  simu_out   <- simulate_nplcm(set_parameter)
  out        <- simu_out$data_nplcm
  out$X      <- rbind(X[i,],X[i,])
  out_list[[i]] <- out
}

# extract cases and controls and combine all the data into one:
data_nplcm_list <- lapply(1:N, function(s) subset_data_nplcm_by_index(out_list[[s]],2-Y[s]))
# the above keeps the hypothetical subject's data (the first subject if Y[s]=1; the second subject if Y[s]=0).
data_nplcm_unordered      <- combine_data_nplcm(data_nplcm_list) 
# put cases on top:
data_nplcm_reg_nest <- subset_data_nplcm_by_index(data_nplcm_unordered,
                                                  order(-data_nplcm_unordered$Y)) 

# Add the continuous covariate to the X data frame
# (get date values based on the simulated "date" variable): 
ENRLDATE.seq <- seq(as.Date('2010-01-06'),as.Date('2010-11-01'),by=1) # 300 days
data_nplcm_reg_nest$X$ENRLDATE <- ENRLDATE.seq[data_nplcm_reg_nest$X$date] # actual date names.

## rename the covariates
colnames(data_nplcm_reg_nest$X) <- c("SITE","DATE","std_date","ENRLDATE")



# ##save the simulated data to the R package for illustration:
# save(data_nplcm_reg_nest, file = "data/data_nplcm_reg_nest.rda", compress = "xz")
