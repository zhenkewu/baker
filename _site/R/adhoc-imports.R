# copied from https://github.com/suyusung/R2jags/blob/master/R/jags.sims.R
# copy from R2WinBUGS
.decode.parameter.name <- function (a) 
{
    left.bracket <- regexpr("[[]", a)
    if (left.bracket == -1) {
        root <- a
        dimension <- 0
        indexes <- NA
    }
    else {
        root <- substring(a, 1, left.bracket - 1)
        right.bracket <- regexpr("[]]", a)
        a <- substring(a, left.bracket + 1, right.bracket - 1)
        indexes <- as.numeric(unlist(strsplit(a, ",")))
        dimension <- length(indexes)
    }
    list(root = root, dimension = dimension, indexes = indexes)
}


#' @importFrom stats var
#' @importFrom utils read.table
jags.sims <- function (parameters.to.save, n.chains, n.iter, n.burnin, n.thin, 
  DIC = TRUE) 
{
  
  #require(R2WinBUGS)
  sims.files <- paste("CODAchain", 1:n.chains, ".txt", sep = "")
  index <- read.table("CODAindex.txt", header = FALSE)#, sep = "\t")
  if (is.R()) {
      parameter.names <- as.vector(index[, 1])
      n.keep <- index[1, 3] - index[1, 2] + 1
  }
  else {
      parameter.names <- row.names(index)
      n.keep <- index[1, 2] - index[1, 1] + 1
  }
  n.parameters <- length(parameter.names)
  n.sims <- n.keep * n.chains
  sims <- matrix(, n.sims, n.parameters)
  sims.array <- array(NA, c(n.keep, n.chains, n.parameters))
  root.long <- character(n.parameters)
  indexes.long <- vector(n.parameters, mode = "list")
  for (i in 1:n.parameters) {
      temp <- .decode.parameter.name(parameter.names[i])
      root.long[i] <- temp$root
      indexes.long[[i]] <- temp$indexes
  }
  n.roots <- length(parameters.to.save)
  left.bracket.short <- as.vector(regexpr("[[]", parameters.to.save))
  right.bracket.short <- as.vector(regexpr("[]]", parameters.to.save))
  root.short <- ifelse(left.bracket.short == -1, parameters.to.save, 
      substring(parameters.to.save, 1, left.bracket.short - 
          1))
  dimension.short <- rep(0, n.roots)
  indexes.short <- vector(n.roots, mode = "list")
  n.indexes.short <- vector(n.roots, mode = "list")
  long.short <- vector(n.roots, mode = "list")
  length.short <- numeric(n.roots)
  for (j in 1:n.roots) {
      long.short[[j]] <- (1:n.parameters)[root.long == root.short[j]]
      length.short[j] <- length(long.short[[j]])
      if (length.short[j] == 0) {
          stop(paste("parameter", root.short[[j]], "is not in the model"))
      }
      else if (length.short[j] > 1) {
          dimension.short[j] <- length(indexes.long[[long.short[[j]][1]]])
          n.indexes.short[[j]] <- numeric(dimension.short[j])
          for (k in 1:dimension.short[j]){
            n.indexes.short[[j]][k] <- length(unique(unlist(lapply(indexes.long[long.short[[j]]], 
              .subset, k))))
          }
          length.short[j] <- prod(n.indexes.short[[j]])
          if (length(long.short[[j]]) != length.short[j]){ 
              stop(paste("error in parameter", root.short[[j]], 
                "in parameters.to.save"))
          }
          indexes.short[[j]] <- as.list(numeric(length.short[j]))
          for (k in 1:length.short[j]){
            indexes.short[[j]][[k]] <- indexes.long[[long.short[[j]][k]]]
          }
      }
  }
  rank.long <- unlist(long.short)
  for (i in 1:n.chains) {
      if (is.R()) {
          sims.i <- scan(sims.files[i], quiet = TRUE)[2 * (1:(n.keep * 
              n.parameters))]
      }
      else {
          sims.i <- scan(sims.files[i])[2 * (1:(n.keep * n.parameters))]
      }
      sims[(n.keep * (i - 1) + 1):(n.keep * i), ] <- sims.i
      sims.array[, i, ] <- sims.i
  }
  dimnames(sims) <- list(NULL, parameter.names)
  dimnames(sims.array) <- list(NULL, NULL, parameter.names)
  summary <- monitor(sims.array, n.chains, keep.all = TRUE)
  last.values <- as.list(numeric(n.chains))
  for (i in 1:n.chains) {
      n.roots.0 <- if (DIC){ n.roots - 1}
      else {n.roots}
      last.values[[i]] <- as.list(numeric(n.roots.0))
      names(last.values[[i]]) <- root.short[1:n.roots.0]
      for (j in 1:n.roots.0) {
          if (dimension.short[j] <= 1) {
              last.values[[i]][[j]] <- sims.array[n.keep, i, 
                long.short[[j]]]
              names(last.values[[i]][[j]]) <- NULL
          }
          else{
            last.values[[i]][[j]] <- aperm(array(sims.array[n.keep, 
                i, long.short[[j]]], rev(n.indexes.short[[j]])), 
                dimension.short[j]:1)
          }
      }
  }
  sims <- sims[sample(n.sims), , drop = FALSE]
  sims.list <- summary.mean <- summary.sd <- summary.median <- vector(n.roots, 
      mode = "list")
  names(sims.list) <- names(summary.mean) <- names(summary.sd) <- names(summary.median) <- root.short
  for (j in 1:n.roots) {
      if (length.short[j] == 1) {
          sims.list[[j]] <- sims[, long.short[[j]]]
          summary.mean[[j]] <- summary[long.short[[j]], "mean"]
          summary.sd[[j]] <- summary[long.short[[j]], "sd"]
          summary.median[[j]] <- summary[long.short[[j]], "50%"]
      }
      else {
        sims.list[[j]] <- array(sims[, long.short[[j]]], c(n.sims, rev(n.indexes.short[[j]])))#, c(1, (dimension.short[j] + 1):2))
        #sims.list[[j]] <- sims[, long.short[[j]]]
        summary.mean[[j]] <- array(summary[long.short[[j]],"mean"],n.indexes.short[[j]])
        summary.sd[[j]] <- array(summary[long.short[[j]],"sd"],n.indexes.short[[j]])
        summary.median[[j]] <- array(summary[long.short[[j]],"50%"],n.indexes.short[[j]])
#          temp2 <- dimension.short[j]:1
#          sims.list[[j]] <- aperm(array(sims[, long.short[[j]]], 
#              c(n.sims, rev(n.indexes.short[[j]]))), c(1, (dimension.short[j] + 
#              1):2))
#          summary.mean[[j]] <- aperm(array(summary[long.short[[j]], 
#              "mean"], rev(n.indexes.short[[j]])), temp2)
#          summary.sd[[j]] <- aperm(array(summary[long.short[[j]], 
#              "sd"], rev(n.indexes.short[[j]])), temp2)
#          summary.median[[j]] <- aperm(array(summary[long.short[[j]], 
#              "50%"], rev(n.indexes.short[[j]])), temp2)
      }
  }
  summary <- summary[rank.long, ]
  all <- list(n.chains = n.chains, n.iter = n.iter, n.burnin = n.burnin, 
      n.thin = n.thin, n.keep = n.keep, n.sims = n.sims, sims.array = sims.array[, 
          , rank.long, drop = FALSE], sims.list = sims.list, 
      sims.matrix = sims[, rank.long], summary = summary, mean = summary.mean, 
      sd = summary.sd, median = summary.median, root.short = root.short, 
      long.short = long.short, dimension.short = dimension.short, 
      indexes.short = indexes.short, last.values = last.values)
  if (DIC) {
    deviance <- all$sims.array[, , "deviance", drop = FALSE]
    dimnames(deviance) <- NULL
    dim(deviance) <- dim(deviance)[1:2]
    pD <- numeric(n.chains)
    DIC <- numeric(n.chains)
    for (i in 1:n.chains) {
      pD[i] <- var(deviance[, i])/2
      DIC[i] <- mean(deviance[, i]) + pD[i]
    }
    all <- c(all, list(isDIC = TRUE, DICbyR = TRUE, pD = mean(pD), 
                DIC = mean(DIC)))
  }
  else {
    all <- c(all, isDIC = FALSE)
  }
  all
}

if(!is.R()) .subset <- function(x, index) x[index]


#' function to write bugs model (copied from R2WinBUGS)
#' 
#' @param model R / S-PLUS function containing the BUGS model in the BUGS model language, for minor differences see Section Details.
#' @param con passed to writeLines which actually writes the model file
#' @param digits number of significant digits used for WinBUGS input, see formatC
#' 
write.model <- function(model, con = "model.bug", digits = 5)
{
  if (is.R()){
    model.text <- c("model", replaceScientificNotationR(body(model), digits = digits))
    # "[\+\-]?\d*\.?[Ee]?[\+\-]?\d*"
  } else {
    ## In S-PLUS the source code of a function can be obtained with
    ## as.character(function_name).  This omits the "function_name <- function()" piece
    model.text <- paste("model", as.character(model))
  }
  model.text <- gsub("%_%", "", model.text)
  writeLines(model.text, con = con)
}


replaceScientificNotationR <- function(bmodel, digits = 5){
  env <- new.env()
  assign("rSNRidCounter", 0, envir=env)
  replaceID <- function(bmodel, env, digits = 5){
    for(i in seq_along(bmodel)){
      if(length(bmodel[[i]]) == 1){
        if(as.character(bmodel[[i]]) %in% c(":", "[", "[[")) return(bmodel)
        if((typeof(bmodel[[i]]) %in% c("double", "integer")) && ((abs(bmodel[[i]]) < 1e-3) || (abs(bmodel[[i]]) > 1e+4))){
          counter <- get("rSNRidCounter", envir=env) + 1
          assign("rSNRidCounter", counter, envir=env)
          id <- paste("rSNRid", counter, sep="")
          assign(id, formatC(bmodel[[i]], digits=digits, format="E"), envir=env)
          bmodel[[i]] <- id
        }
      } else {
        bmodel[[i]] <- replaceID(bmodel[[i]], env, digits = digits)
      }
    }
    bmodel
  }
  bmodel <- deparse(replaceID(bmodel, env, digits = digits), control = NULL)
  for(i in ls(env)){
    bmodel <- gsub(paste('"', i, '"', sep=''), get(i, envir=env), bmodel, fixed=TRUE)
  }
  bmodel
}



#' @importFrom stats sd
"monitor" <-
  function (a, n.chains=dim(a)[2], trans=NULL, keep.all=FALSE, Rupper.keep=FALSE) {
    
    ## If keep.all=T:  a is a n x m x k array:
    ##   m sequences of length n, k variables measured
    ## If keep.all=F:  a is a 2n x m x k array (first half will be discarded)
    ##
    ## trans is a vector of length k:  "" if no transformation, or "log" or "logit"
    ## (If trans is not defined, it will be set to "log" for parameters that
    ## are all-positive and 0 otherwise.)
    ##
    ## If Rupper.keep=TRUE:  keep Rupper.  (Otherwise don't display it.)
    invlogit <- function (x) {1 / (1 + exp(-x))}
    nparams <- if(length(dim(a)) < 3) 1 else dim(a)[length(dim(a))]
    # Calculation and initialization of the required matrix "output"
    output <- matrix( , ncol = if(n.chains > 1){if(Rupper.keep) 10 else 9} else 7, nrow = nparams)
    if (length(dim(a))==2) a <- array (a, c(dim(a),1))
    if (!keep.all){
      n <- floor(dim(a)[1]/2)
      a <- a[(n+1):(2*n), , , drop = FALSE]
    }
    if (is.null(trans))
      trans <- ifelse ((apply (a<=0, 3, sum))==0, "log", "")
    for (i in 1:nparams){
      # Rupper.keep:  discard Rupper (nobody ever uses it)
      ai <- a[ , , i, drop = FALSE]
      if (trans[i]=="log"){
        conv.p <- conv.par(log(ai), n.chains, Rupper.keep=Rupper.keep) # reason????
        conv.p <- list(quantiles = exp(conv.p$quantiles),
                       confshrink = conv.p$confshrink, n.eff = conv.p$n.eff)
      }
      else if (trans[i]=="logit"){
        if (!is.R()){
          logit <- function (x) { log(x /(1- x)) }
        }    
        conv.p <- conv.par(logit(ai), n.chains, Rupper.keep=Rupper.keep)
        conv.p <- list(quantiles = invlogit(conv.p$quantiles),
                       confshrink = conv.p$confshrink, n.eff = conv.p$n.eff)
      }
      else conv.p <- conv.par(ai, n.chains, Rupper.keep=Rupper.keep)
      output[i, ] <- c(mean(ai), sd(as.vector(ai)),
                       conv.p$quantiles, 
                       if(n.chains > 1) conv.p$confshrink, 
                       if(n.chains > 1) round(conv.p$n.eff, min(0, 1 - floor(log10(conv.p$n.eff))))
      )
    }
    if(n.chains > 1)
      dimnames(output) <- list(dimnames(a)[[3]], c("mean","sd",
                                                   "2.5%","25%","50%","75%","97.5%", "Rhat", if(Rupper.keep) "Rupper","n.eff"))
    else
      dimnames(output) <- list(dimnames(a)[[3]], c("mean","sd",
                                                   "2.5%","25%","50%","75%","97.5%"))
    return (output)
  }

#' @importFrom stats qf
"conv.par" <-
  function (x, n.chains, Rupper.keep = TRUE) {
    m <- ncol(x)
    n <- nrow(x)
    
    # We compute the following statistics:
    #
    #  xdot:  vector of sequence means
    #  s2:  vector of sequence sample variances (dividing by n-1)
    #  W = mean(s2):  within MS
    #  B = n*var(xdot):  between MS.
    #  muhat = mean(xdot):  grand mean; unbiased under strong stationarity
    #  varW = var(s2)/m:  estimated sampling var of W
    #  varB = B^2 * 2/(m+1):  estimated sampling var of B
    #  covWB = (n/m)*(cov(s2,xdot^2) - 2*muhat*cov(s^2,xdot)):
    #                                               estimated sampling cov(W,B)
    #  sig2hat = ((n-1)/n))*W + (1/n)*B:  estimate of sig2; unbiased under
    #                                               strong stationarity
    #  quantiles:  emipirical quantiles from last half of simulated sequences
    
    xdot <- apply(x,2,mean)
    muhat <- mean(xdot)
    s2 <- apply(x,2,var)
    W <- mean(s2)
    quantiles <- quantile (as.vector(x), probs=c(.025,.25,.5,.75,.975))
    
    if ((W > 1.e-8) && (n.chains > 1)) {            # non-degenerate case
      
      B <- n*var(xdot)
      varW <- var(s2)/m
      varB <- B^2 * 2/(m-1)
      covWB <- (n/m)*(var(s2, xdot^2) - 2*muhat*var(s2, xdot))
      sig2hat <- ((n-1)*W + B)/n
      
      # Posterior interval post.range combines all uncertainties
      # in a t interval with center muhat, scale sqrt(postvar),
      # and postvar.df degrees of freedom.
      #
      #       postvar = sig2hat + B/(mn):  variance for the posterior interval
      #                               The B/(mn) term is there because of the
      #                               sampling variance of muhat.
      #       varpostvar:  estimated sampling variance of postvar
      
      postvar <- sig2hat + B/(m*n)
      varpostvar <- max(0, 
                        (((n-1)^2) * varW + (1 + 1/m)^2 * varB + 2 * (n-1) * (1 + 1/m) * covWB) / n^2)
      post.df <- min(2*(postvar^2/varpostvar), 1000)
      
      # Estimated potential scale reduction (that would be achieved by
      # continuing simulations forever) has two components:  an estimate and
      # an approx. 97.5% upper bound.
      #
      # confshrink = sqrt(postvar/W),
      #     multiplied by sqrt(df/(df-2)) as an adjustment for the
      ###      CHANGED TO sqrt((df+3)/(df+1))
      #     width of the t-interval with df degrees of freedom.
      #
      # postvar/W = (n-1)/n + (1+1/m)(1/n)(B/W); we approximate the sampling dist.
      # of (B/W) by an F distribution, with degrees of freedom estimated
      # from the approximate chi-squared sampling dists for B and W.  (The
      # F approximation assumes that the sampling dists of B and W are independent;
      # if they are positively correlated, the approximation is conservative.)
      
      confshrink.range <- postvar/W
      if(Rupper.keep){
        varlo.df <- 2*(W^2/varW) 
        confshrink.range <- c(confshrink.range, 
                              (n-1)/n + (1+1/m)*(1/n)*(B/W) * qf(.975, m-1, varlo.df))
      }
      confshrink.range <- sqrt(confshrink.range * (post.df+3) / (post.df+1))
      
      # Calculate effective sample size:  m*n*min(sigma.hat^2/B,1)
      # This is a crude measure of sample size because it relies on the between
      # variance, B, which can only be estimated with m degrees of freedom.
      
      n.eff <- m*n*min(sig2hat/B,1)
      list(quantiles=quantiles, confshrink=confshrink.range,
           n.eff=n.eff)
      
    }
    else {      # degenerate case:  all entries in "data matrix" are identical
      list (quantiles=quantiles, confshrink = rep(1, Rupper.keep + 1),
            n.eff=1)
      
    }
  }


