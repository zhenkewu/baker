set.seed(2)
# true TPR and FPR for measurement data 
PsiBS0   <- cbind(c(0.25,0.25,0.2,0.15,0.15,0.15), c(0.2,0.2,0.25,0.1,0.1,0.1))

ThetaBS0 <- cbind(c(0.95,0.9,0.9,0.9,0.9,0.9), c(0.95,0.9,0.9,0.9,0.9,0.9))

ThetaSS_withNA <- c(0.15,0.1,NA,NA,NA,NA)
PsiSS_withNA <- c(0,0,NA,NA,NA,NA)

# true subclass weights 
lambda   <- c(0.5,0.5) 
eta      <- c(0,1)

# the following paramter names are set using names in the 'baker' package:
set_parameter_noreg <- list(
  cause_list      = cause_list,
  etiology        = c(0.5,0.2,0.15,0.05,0.05,0.05), 
  pathogen_BrS    = LETTERS[1:L.CAUSE],
  pathogen_SS     = LETTERS[1:L.CAUSE][!is.na(ThetaSS_withNA)],
  meas_nm         = list(MBS = c("MBS1"),MSS=c("MSS1")), # a single slice of Bronze Standard (BrS) data
  Lambda          = lambda,           
  Eta             = t(replicate(J.BrS,eta)), 
  PsiBS           = PsiBS0,
  ThetaBS         = ThetaBS0,
  PsiSS           = PsiSS_withNA[!is.na(PsiSS_withNA)], # FPR for SS measures; 0 for perfect specificity
  ThetaSS         = ThetaSS_withNA[!is.na(ThetaSS_withNA)], # TPRS for MSS1
  Nd              =     Nd, 
  Nu              =     Nu 
)

simu_out_noreg   <- simulate_nplcm(set_parameter_noreg)
data_nplcm_noreg <- simu_out_noreg$data_nplcm


##save the simulated data to the R package for illustration:
save(data_nplcm_noreg, file = "data/data_nplcm_noreg.rda", compress = "xz")