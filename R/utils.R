


#' Create new folder name
#'
#' @param parent_path The parent directory where to put the new folder
#' @param parameter_names The parameters that distinguish this folder's scenario
#' @param parameter_vals The actual parameter values
#'
#' @return A string for folder name
#' @export
#'
make_foldername <-
  function(parent_path,parameter_names,parameter_vals) {
    subfolder <-
      paste(parameter_names,parameter_vals,collapse = "_",sep = "=")
    res       <- paste(parent_path,subfolder,sep = "\\")
    res
  }

#' Create new file name
#'
#'
#' @param parameter_names The parameters that distinguish this folder's scenario
#' @param parameter_vals The actual parameter values
#' @param format The suffix ".XXX" in the end to specify the file format
#'
#'
#' @return A string for file name
#'
#'
#' @export
#'
make_filename <-
  function(parameter_names,parameter_vals,format) {
    res1 <- paste(parameter_names,parameter_vals,collapse = "_",sep = "=")
    res  <- paste(res1,format,sep = ".")
    res
  }

#' logit function
#'
#' @param p Probability between 0 and 1
#' @return A real number
#'
#' @export
logit <- function(p)
  log(p) - log(1 - p)

#' expit function
#'
#' @param x A real number
#' @return a Probability between 0 and 1
#' @export
expit <- function(x)
  1 / (1 + exp(-x))




#' Sample a vector of Bernoulli variables.
#'
#' Sample a vector of Bernoulli variables with higher speed
#' (same length with \code{"p"}).
#' The Bernoulli random variables can have different means.
#'
#' @param p A vector of probabilities, each being the head probability
#' of an independent coin toss
#'
#' @return A vector of 1s (head) and 0s (tail)
#' @export
rvbern <-
  function(p) {
    U  <- runif(length(p),0,1)
    res <- (U < p) + 0
    res
  }






#' Convert 0/1 binary coded sequence into decimal digits
#'
#' Useful when try to list all the binary patterns. One can group the binary
#' sequences according to their equivalent decimal values.
#'
#' @param binary_vector a binary number
#' @return a decimal number
#' @export
#'
bin2dec <- function(binary_vector) {
  sum(2 ^ (which(rev(binary_vector) == TRUE) - 1))
}




#' Convert names of pathogen/combinations into 0/1 coding
#'
#' @param pathogen_name The allowed pathogen name (can be a combination of pathogens in "pathlist")
#' @param pathogen_list The complete list of pathogen names
#'
#' @return A 1 by length(pathlist) matrix of binary code (usually for pathogen presence/absence)
#'
#' @examples
#' symb2I("A",c("A","B","C"))
#' symb2I("A+B",c("A","B","C"))
#' symb2I("NoA",c("A","B","C"))
#' symb2I(c("A","B+C"),c("A","B","C")) # gives a 2 by 3 matrix.
#'
#'
#' @export
#'
symb2I <-
  function(pathogen_name,pathogen_list) {
    J <- length(pathogen_list)
    splited <- strsplit(pathogen_name,split = "+",fixed = TRUE)
    deploy  <- function(inst,J) {
      if (inst[1] == "NoA") {
        rep(0,J)
      }else{
        sapply(inst,grep,x = pathogen_list)
      }
    }
    res <- lapply(splited,deploy,J = J)
    for (l in seq_along(res)) {
      any_integer_0 <- sum(sapply(res[[l]],function(v)
        length(v) == 0)) > 0
      if (any_integer_0) {
        stop(
          paste0(
            "\n==",l,"-th cause, \n",pathogen_name[l],",
            \n has pathogen(s) not in the overall pathogen list!=="
          )
          )
      }
      }
    
    nc  <- length(res)
    matres <- t(sapply(1:nc,function(i) {
      tempres <- rep(0,J)
      tempres[res[[i]]] <- 1
      tempres
    }))
    matres
    }

#symb2I("A+D",c("A","B","C")) # will return error.

#' Convert 0/1 coding to pathogen/combinations
#'
#' Reverse to \code{\link{symb2I}}
#' @param binary_code Binary indictors for pathogens
#' @param pathogen_list The complete list of pathogen names
#'
#' @return The name of pathogen or pathogen combination indicated by "code"
#'
#' @examples
#' I2symb("001",c("A","B","C"))
#' I2symb("000",c("A","B","C"))
#'
#' @export
I2symb <- function(binary_code,pathogen_list) {
  ind <- grep("1",strsplit(binary_code,split = "")[[1]])
  res <-
    ifelse(length(ind) == 0,"NoA",paste(pathogen_list[ind],collapse = "+"))
  res
}



#' Convert a matrix of binary indicators to categorial variables
#'
#' @param binary_mat The matrix of binary indicators. Rows for subjects, columns for pathogens in the \code{"pathogen.list"}
#' @param cause_list The list of causes
#' @param pathogen_list The complete list of pathogen names
#'
#' @return A vector of categorical variables. Its length equals the length of \code{"allowed.list"}
#'
#' @examples
#'
#' Imat2cat(rbind(diag(3),c(1,1,0),c(0,0,0)),c("A","B","C","A+B","NoA"),c("A","B","C"))
#' @export
Imat2cat <- function(binary_mat,cause_list,pathogen_list) {
  known_code = apply(binary_mat,1,function(v)
    paste(v,collapse = ""))
  known_symb = sapply(known_code,I2symb,pathogen_list)
  if (sum(known_symb %in% cause_list == FALSE) > 0) {
    stop("Some binary pattern in 'binary_mat' is not included by 'cause_list'.")
  }else {
    known_Icat = sapply(known_symb,function(s)
      which(cause_list == s))
    return(known_Icat)
  }
}

#' Plot beta density
#' @param a The first parameter
#' @param b The second parameter
#' @return None
#'
#' @export
beta_plot = function(a,b) {
  x = seq(0,1,by = 0.001)
  y = dbeta(x,a,b)
  plot(x,y,type = "l",main = paste0("a=",a,",b=",b))
}


#' Get package from CRAN website
#'
#' @param pckg package name
#'
#' @return None
#' @export
getPckg <- function(pckg) {
  install.packages(pckg, repos = "http://cran.r-project.org")
}


#' Convert factor to numeric without losing information on the label
#'
#' @param f A factor
#'
#' @return A numeric vector
#'
#' @export
unfactor <- function(f) {
  as.numeric(levels(f))[f]
}


#' Order rows of data frame according to variable combinations.
#'
#'
#' Author: Kevin Wright\cr
#' http://tolstoy.newcastle.edu.au/R/help/04/09/4300.html
#' Some ideas from Andy Liaw\cr
#' http://tolstoy.newcastle.edu.au/R/help/04/07/1076.html\cr
#' Use + for ascending, - for decending.\cr
#' Sorting is left to right in the formula\cr
#' Useage is either of the following:\cr
#' sort.data.frame(~Block-Variety,Oats)\cr
#' sort.data.frame(Oats,~-Variety+Block)\cr
#'
#'
#'@param form A Formula. See example in \strong{Description}.
#'@param dat Data frame to be ordered.
#'
#'@return Ordered data frame.
#'
#'@export

sort_data_frame <- function(form,dat) {
  # Author: Kevin Wright
  # http://tolstoy.newcastle.edu.au/R/help/04/09/4300.html
  # Some ideas from Andy Liaw
  # http://tolstoy.newcastle.edu.au/R/help/04/07/1076.html
  # Use + for ascending, - for decending.
  # Sorting is left to right in the formula
  # Useage is either of the following:
  # sort.data.frame(~Block-Variety,Oats)
  # sort.data.frame(Oats,~-Variety+Block)
  
  # If dat is the formula, then switch form and dat
  if (inherits(dat,"formula")) {
    f = dat
    dat = form
    form = f
  }
  if (form[[1]] != "~") {
    stop("Formula must be one-sided.")
  }
  # Make the formula into character and remove spaces
  formc <- as.character(form[2])
  formc <- gsub(" ","",formc)
  # If the first character is not + or -, add +
  if (!is.element(substring(formc,1,1),c("+","-"))) {
    formc <- paste("+",formc,sep = "")
  }
  # Extract the variables from the formula
  vars <- unlist(strsplit(formc, "[\\+\\-]"))
  vars <- vars[vars != ""] # Remove spurious "" terms
  # Build a list of arguments to pass to "order" function
  calllist <- list()
  pos = 1 # Position of + or -
  for (i in 1:length(vars)) {
    varsign <- substring(formc,pos,pos)
    pos <- pos + 1 + nchar(vars[i])
    if (is.factor(dat[,vars[i]])) {
      if (varsign == "-")
        calllist[[i]] <- -rank(dat[,vars[i]])
      else
        calllist[[i]] <- rank(dat[,vars[i]])
    }
    else {
      if (varsign == "-")
        calllist[[i]] <- -dat[,vars[i]]
      else
        calllist[[i]] <- dat[,vars[i]]
    }
  }
  dat[do.call("order",calllist),]
}

#' Reorder the measurement dimensions to match the order for display
#'
#' @param disp_order The vector of names to be displayed (order matters)
#' @param raw_nm The vector of names from raw measurements (order matters)
#'
#' @return A permuted vector from 1 to \code{length(raw_nm)}. For example, if
#' its first element is 3, it means that the 3rd pathogen in \code{raw_nm}
#' should be arranged to the first in the raw measurements.
#'
#' @examples
#'   disp_order <- c("B","E","D","C","F","A")
#'   raw_nm <- c("C","A","E")
#'   my_reorder(disp_order,raw_nm)
#'
#'
#' @export
my_reorder <- function(disp_order,raw_nm) {
  #disp_order <- display_order
  #raw_nm <- pathogen_BrS
  res <- rep(NA,length(raw_nm))
  for (i in seq_along(raw_nm)) {
    res[i] <- which(disp_order == raw_nm[i])
  }
  order(res)
}



#' calculate pairwise log odds ratios
#'
#' Case at upper triangle; control at lower triangle
#'
#' @param MBS.case Case Bronze-Standard (BrS) data
#' @param MBS.ctrl Control Bronze-Standard (BrS) data
#'
#' @export
#'

logOR <- function(MBS.case,MBS.ctrl) {
  JBrS <- ncol(MBS.case)
  logORmat       <- matrix(NA,nrow = JBrS,ncol = JBrS)
  logORmat.se    <- matrix(NA,nrow = JBrS,ncol = JBrS)
  for (j2 in 1:(JBrS - 1)) {
    #case (j2,j1); ctrl (j1,j2).
    for (j1 in (j2 + 1):JBrS) {
      # cases: (upper triangle)
      x <- MBS.case[,j2]
      y <- MBS.case[,j1]
      
      fit <- glm(y ~ x,family = binomial(link = "logit"))
      
      if ("x" %in% rownames(summary(fit)$coef)) {
        logORmat[j2,j1] = round(summary(fit)$coef["x",1],3)
        logORmat.se[j2,j1] = round(summary(fit)$coef["x",2],3)
      }
      # controls: (lower triangle)
      x <- MBS.ctrl[,j2]
      y <- MBS.ctrl[,j1]
      
      fit <- glm(y ~ x,family = binomial(link = "logit"))
      
      if ("x" %in% rownames(summary(fit)$coef)) {
        logORmat[j1,j2] = round(summary(fit)$coef["x",1],3)
        logORmat.se[j1,j2] = round(summary(fit)$coef["x",2],3)
      }
    }
  }
  
  #cell.num = logORmat/logORmat.se
  tmp       = logORmat
  tmp[abs(logORmat.se) > 10] = NA
  
  tmp.se = logORmat.se
  tmp.se[abs(logORmat.se) > 10] = NA
  
  
  res <- list(tmp,tmp.se)
  names(res) <- c("logOR","logOR.se")
  return(res)
}







#' Visualize matrix for a quantity measured on cases and controls (a single number)
#'
#' Special to case-control visualization: upper right for cases, lower left
#' for controls.
#'
#' @param mat matrix of values: upper for cases, lower for controls;
#' @param dim_names names of the columns, from left to right. It is also the
#' names of the rows, from bottom to top. Default is 1 through \code{ncol(mat)};
#' @param cell_metrics the meaning of number in every cell;
#' @param folding_line Default is \code{TRUE} for adding dashed major diagnoal
#' line.
#' @param axes plot axes; default is \code{FALSE};
#' @param xlab label for x-axis
#' @param ylab label for y-axis
#' @param asp aspect ratio; default is \code{1} to ensure square shape
#' @param title text for the figure
#'
#' @export
visualize_case_control_matrix <- function(mat, dim_names = ncol(mat),
                                          cell_metrics = "",folding_line = TRUE,
                                          axes = FALSE, xlab = "",ylab = "",
                                          asp = 1,title = "") {
  n = nrow(mat)
  J = n
  # size of the numbers in the boxes:
  cex_main = min(2,20 / n)
  cex_se  = min(1.5,15 / n)
  
  par(mar = c(0, 0, 5, 0), bg = "white",xpd = TRUE)
  plot(
    c(0, n + 0.8), c(0, n + 0.8), axes = axes, xlab = "",
    ylab = "", asp = 1, type = "n"
  )
  ##add grid
  segments(rep(0.5, n + 1), 0.5 + 0:n, rep(n + 0.5, n + 1),
           0.5 + 0:n, col = "gray")
  segments(0.5 + 0:n, rep(0.5, n + 1), 0.5 + 0:n, rep(n + 0.5,
                                                      n), col = "gray")
  mat.txt <- round(t(mat)[,n:1],1)
  mat.txt3 <- round(t(mat)[,n:1],3)
  
  # add meaning of the number in a cell:
  text(1,J,cell_metrics,cex = cex_main / 2)
  
  for (i in 1:n) {
    for (j in 1:n) {
      abs.mat <- abs(mat.txt[i,j])
      if ((!is.na(abs.mat)) && abs.mat > 2) {
        text(
          i,j,round(mat.txt3[i,j],1),
          col = ifelse(mat.txt3[i,j] > 0,"red","blue"),cex = cex_main
        )
      }
      
    }
  }
  
  if (folding_line) {
    # diagonal line:
    segments(
      0.5 + 1,.5 + n - 1,.5 + n,0.5,col = "black",lty = 3,lwd = 3
    )
  }
  
  # put pathogen names on rows and columns:
  for (s in 1:J) {
    text(
      -0,J - s + 1,paste0(dim_names[s],":(",s,")"),cex = min(1.5,20 / J),adj =
        1
    )
    text(
      s,J + 0.7,paste0("(",s,"):",dim_names[s]),cex = min(1.5,20 / J),srt = 45,adj =
        0
    )
  }
  # labels for cases and controls:
  text(J + 1,J / 2,"cases",cex = 2,srt = -90)
  text(J / 2,0,"controls",cex = 2)
}

#' convert 'NA' to '.'
#'
#' @param s A string of characters that may contain "NA"
#' @return A string of characters without 'NA'
#' @export
NA2dot <- function(s) {
  gsub("NA",".",s,fixed = TRUE)
}


#' Pick parameters in the Beta distribution to match the specified range
#'
#' \code{beta_parms_from_quantiles} produces prior Beta parameters for
#'  the true positive rates (TPR)
#'
#' @param q A vector of lower and upper bounds, in which Beta distribution
#' will have quantiles specified by \code{p}. For example, \code{q=c(0.5,0.99)}
#' @param p The lower and upper quantiles of the range one wants to specify.
#' @param precision Approximation precisions.
#' @param derivative.epsilon Precision of calculating derivative.
#' @param start.with.normal.approx Default is \code{TRUE}, for normal approximation.
#' @param start Starting values of beta parameters.
#' @param plot Default is \code{FALSE} to suppress plotting of the beta density,
#' otherwise, set to \code{TRUE}.
#'
#' @return A list containing the selected Beta prameters \code{a}, and \code{b}.
#' Other elements of the list include some details about the computations involved
#' in finding \code{a} and \code{b}.
#'
#' @references \url{http://www.medicine.mcgill.ca/epidemiology/
#'                  Joseph/PBelisle/BetaParmsFromQuantiles.html}
#' @export
#'
beta_parms_from_quantiles <- function(q, p = c(0.025,0.975),
                                      precision = 0.001,
                                      derivative.epsilon = 1e-3,
                                      start.with.normal.approx = T,
                                      start = c(1, 1),
                                      plot = F) {
  # Version 1.2.2 (December 2012)
  #
  # Function developed by
  # Lawrence Joseph and Patrick Belisle
  # Division of Clinical Epidemiology
  # Montreal General Hospital
  # Montreal, Qc, Can
  #
  # patrick.belisle@clinepi.mcgill.ca
  # http://www.medicine.mcgill.ca/epidemiology/Joseph/PBelisle/BetaParmsFromQuantiles.html
  #
  # Please refer to our webpage for details on each argument.
  
  f <-
    function(x, theta) {
      dbeta(x, shape1 = theta[1], shape2 = theta[2])
    }
  F.inv <-
    function(x, theta) {
      qbeta(x, shape1 = theta[1], shape2 = theta[2])
    }
  f.cum <-
    function(x, theta) {
      pbeta(x, shape1 = theta[1], shape2 = theta[2])
    }
  f.mode <- function(theta) {
    a <- theta[1]; b <- theta[2];
    mode <- ifelse(a > 1, (a - 1) / (a + b - 2), NA); mode
  }
  theta.from.moments <-
    function(m, v) {
      a <- m * m * (1 - m) / v - m; b <- a * (1 / m - 1); c(a, b)
    }
  plot.xlim <- c(0, 1)
  
  dens.label <- 'dbeta'
  parms.names <- c('a', 'b')
  
  if (length(p) != 2)
    stop("Vector of probabilities p must be of length 2.")
  if (length(q) != 2)
    stop("Vector of quantiles q must be of length 2.")
  p <- sort(p); q <- sort(q)
  
  #_____________________________________________________________________________________________________
  
  print.area.text <- function(p, p.check, q, f, f.cum, F.inv,
                              theta, mode, cex, plot.xlim, M = 30, M0 =
                                50)
  {
    par.usr <- par('usr')
    par.din <- par('din')
    
    p.string <-
      as.character(round(c(0,1) + c(1,-1) * p.check, digits = 4))
    str.width <- strwidth(p.string, cex = cex)
    str.height <- strheight("0", cex = cex)
    
    J <- matrix(1, nrow = M0, ncol = 1)
    
    x.units.1in <- diff(par.usr[c(1,2)]) / par.din[1]
    y.units.1in <- diff(par.usr[c(3,4)]) / par.din[2]
    aspect.ratio <- y.units.1in / x.units.1in
    
    # --- left area  -----------------------------------------------------------
    
    scatter.xlim <- c(max(plot.xlim[1], par.usr[1]), q[1])
    scatter.ylim <- c(0, par.usr[4])
    x <- seq(from = scatter.xlim[1], to = scatter.xlim[2], length = M)
    y <- seq(from = scatter.ylim[1], to = scatter.ylim[2], length = M)
    x.grid.index <- rep(seq(M), M)
    y.grid.index <- rep(seq(M), rep(M, M))
    
    grid.df <- f(x, theta)
    
    # Estimate mass center
    tmp.p <- seq(from = 0, to = p[1], length = M0)
    tmp.x <- F.inv(tmp.p, theta)
    h <- f(tmp.x, theta)
    mass.center <-
      c(mean(tmp.x), sum(h[-1] * diff(tmp.x)) / diff(range(tmp.x)))
    
    # Identify points under the curve
    # (to eliminate them from the list of candidates)
    gridpoint.under.the.curve <-
      y[y.grid.index] <= grid.df[x.grid.index]
    w <- which(gridpoint.under.the.curve)
    x <- x[x.grid.index]; y <- y[y.grid.index]
    if (length(w)) {
      x <- x[-w]; y <- y[-w]
    }
    
    # Eliminate points to the right of the mode, if any
    w <- which(x > mode)
    if (length(w)) {
      x <- x[-w]; y <- y[-w]
    }
    
    # Eliminate points for which the text would fall out of the plot area
    w <-
      which((par.usr[1] + str.width[1]) <= x &
              (y + str.height) <= par.usr[4])
    x <- x[w]; y <- y[w]
    
    # For each height, eliminate the closest point to the curve
    # (we want to stay away from the curve to preserve readability)
    w <- which(!duplicated(y, fromLast = T))
    if (length(w)) {
      x <- x[-w]; y <- y[-w]
    }
    
    # For each point, compute distance from mass center and pick the closest point
    d <-
      ((x - mass.center[1]) ^ 2) + ((y - mass.center[2]) / aspect.ratio) ^ 2
    w <- which.min(d)
    x <- x[w]; y <- y[w]
    
    if (length(x))
    {
      text(
        x, y, labels = p.string[1], adj = c(1,0), col = 'gray', cex = cex
      )
    }
    else
    {
      text(
        plot.xlim[1], mean(par.usr[c(3,4)]), labels = p.string[1],
        col = 'gray', cex = cex, srt = 90, adj = c(1,0)
      )
    }
    
    # --- right area  ----------------------------------------------------------
    
    scatter.xlim <- c(q[2], plot.xlim[2])
    scatter.ylim <- c(0, par.usr[4])
    x <- seq(from = scatter.xlim[1], to = scatter.xlim[2], length = M)
    y <- seq(from = scatter.ylim[1], to = scatter.ylim[2], length = M)
    x.grid.index <- rep(seq(M), M)
    y.grid.index <- rep(seq(M), rep(M, M))
    grid.df <- f(x, theta)
    
    # Estimate mass center
    tmp.p <-
      seq(
        from = p[2], to = f.cum(plot.xlim[2], theta), length = M0
      )
    tmp.x <- F.inv(tmp.p, theta)
    h <- f(tmp.x, theta)
    mass.center <-
      c(mean(tmp.x), sum(h[-length(h)] * diff(tmp.x)) / diff(range(tmp.x)))
    
    # Identify points under the curve
    # (to eliminate them from the list of candidates)
    gridpoint.under.the.curve <-
      y[y.grid.index] <= grid.df[x.grid.index]
    w <- which(gridpoint.under.the.curve)
    x <- x[x.grid.index]; y <- y[y.grid.index]
    if (length(w)) {
      x <- x[-w]; y <- y[-w]
    }
    
    # Eliminate points to the left of the mode, if any
    w <- which(x < mode)
    if (length(w)) {
      x <- x[-w]; y <- y[-w]
    }
    
    # Eliminate points for which the text would fall out of the plot area
    w <-
      which((par.usr[2] - str.width[2]) >= x &
              (y + str.height) <= par.usr[4])
    x <- x[w]; y <- y[w]
    
    # For each height, eliminate the closest point to the curve
    # (we want to stay away from the curve to preserve readability)
    w <- which(!duplicated(y))
    if (length(w)) {
      x <- x[-w]; y <- y[-w]
    }
    
    # For each point, compute distance from mass center and pick the closest point
    d <-
      ((x - mass.center[1]) ^ 2) + ((y - mass.center[2]) / aspect.ratio) ^ 2
    w <- which.min(d)
    x <- x[w]; y <- y[w]
    
    if (length(x))
    {
      text(
        x, y, labels = p.string[2], adj = c(0,0), col = 'gray', cex = cex
      )
    }
    else
    {
      text(
        plot.xlim[2], mean(par.usr[c(3,4)]), labels = p.string[2],
        col = 'gray', cex = cex, srt = -90, adj = c(1,0)
      )
    }
  }
  
  # ......................................................................................................................................
  
  Newton.Raphson <- function(derivative.epsilon, precision,
                             f.cum, p, q, theta.from.moments,
                             start.with.normal.approx, start)
  {
    Hessian <- matrix(NA, 2, 2)
    
    if (start.with.normal.approx)
    {
      # Probably not a very good universal choice,
      # but proved good in most cases in practice
      m <-  diff(q) / diff(p) * (0.5 - p[1]) + q[1]
      v <- (diff(q) / diff(qnorm(p))) ^ 2
      theta <- theta.from.moments(m, v)
    }
    else
      theta <- start
    
    
    change <- precision + 1
    niter <- 0
    # Newton-Raphson multivariate algorithm
    while (max(abs(change)) > precision)
    {
      Hessian[,1] <- (f.cum(q, theta) - f.cum(q, theta -
                                                c(derivative.epsilon, 0))) / derivative.epsilon
      Hessian[,2] <- (f.cum(q, theta) - f.cum(q, theta -
                                                c(0, derivative.epsilon))) / derivative.epsilon
      
      f <- f.cum(q, theta) - p
      change <- solve(Hessian) %*% f
      last.theta <- theta
      theta <- last.theta - change
      
      # If we step out of limits, reduce change
      
      if (any(theta < 0))
      {
        k <- min(last.theta / change)
        theta <- last.theta - k / 2 * change
      }
      
      niter <- niter + 1
    }
    
    list(
      theta = as.vector(theta), niter = niter, last.change = as.vector(change)
    )
  }
  
  # ...............................................................................................................
  
  plot.density <- function(p, q, f, f.cum, F.inv, mode, theta,
                           plot.xlim, dens.label, parms.names, cex)
  {
    if (length(plot.xlim) == 0)
    {
      plot.xlim <- F.inv(c(0, 1), theta)
      
      if (is.infinite(plot.xlim[1]))
      {
        tmp <- min(c(0.001, p[1] / 10))
        plot.xlim[1] <- F.inv(tmp, theta)
      }
      
      if (is.infinite(plot.xlim[2]))
      {
        tmp <- max(c(0.999, 1 - (1 - p[2]) / 10))
        plot.xlim[2] <- F.inv(tmp, theta)
      }
    }
    plot.xlim <- sort(plot.xlim)
    
    
    x <- seq(
      from = min(plot.xlim), to = max(plot.xlim), length = 1000
    )
    h <- f(x, theta)
    x0 <- x; f0 <- h
    ylab <- paste(
      c(
        dens.label, '(x, ', parms.names[1], ' = ',
        round(theta[1], digits = 5), ', ', parms.names[2], ' = ',
        round(theta[2], digits = 5), ')'
      ), collapse = ''
    )
    plot(x, h, type = 'l', ylab = ylab)
    
    # fill in area on the left side of the distribution
    x <- seq(from = plot.xlim[1], to = q[1], length = 1000)
    y <- f(x, theta)
    x <- c(x, q[1], plot.xlim[1]); y <- c(y, 0, 0)
    polygon(x, y, col = 'lightgrey', border = 'lightgray')
    # fill in area on the right side of the distribution
    x <- seq(from = max(plot.xlim), to = q[2], length = 1000)
    y <- f(x, theta)
    x <- c(x, q[2], plot.xlim[2]); y <- c(y, 0, 0)
    polygon(x, y, col = 'lightgrey', border = 'lightgray')
    # draw distrn again
    points(x0, f0, type = 'l')
    h <- f(q, theta)
    points(rep(q[1], 2), c(0, h[1]), type = 'l', col = 'orange')
    points(rep(q[2], 2), c(0, h[2]), type = 'l', col = 'orange')
    # place text on both ends areas
    print.area.text(p, p.check, q, f, f.cum, F.inv, theta, mode, cex, plot.xlim)
    
    xaxp <- par("xaxp")
    x.ticks <- seq(from = xaxp[1], to = xaxp[2], length = xaxp[3] + 1)
    q2print <-
      as.double(setdiff(as.character(q), as.character(x.ticks)))
    
    mtext(
      q2print, side = 1, col = 'orange', at = q2print, cex = 0.6, line = 2.1
    )
    points(q, rep(par('usr')[3] + 0.15 * par('cxy')[2], 2), pch = 17, col =
             'orange')
  }
  
  #________________________________________________________________________________________________________________
  
  
  parms <-
    Newton.Raphson(
      derivative.epsilon, precision, f.cum, p, q,
      theta.from.moments, start.with.normal.approx, start =
        start
    )
  p.check <- f.cum(q, parms$theta)
  
  if (plot)
    plot.density(
      p, q, f, f.cum, F.inv, f.mode(parms$theta),
      parms$theta, plot.xlim, dens.label, parms.names, 0.8
    )
  
  list(
    a = parms$theta[1], b = parms$theta[2], last.change = parms$last.change,
    niter = parms$niter, q = q, p = p, p.check = p.check
  )
}







#' Load or install a package
#'
#'\code{load_or_install} checks if a package is installed,
#' loads it if it is, and installs it if not.
#'
#' @param package_names A vector of package names
#' @param repos URL for downloading the packages. Default is
#' \url{http://lib.stat.cmu.edu/R/CRAN}
#'
#' @references Credit:
#'  \url{http://www.vikparuchuri.com/blog/loading-andor-installing-packages/}
#' @return No message if it successfully loads the specified packages; Error
#'  if such a package does not exist.
#'
#' @export
#'
#'
load_or_install <-
  function(package_names,repos = "http://lib.stat.cmu.edu/R/CRAN") {
    is_installed <-
      function(mypkg)
        is.element(mypkg, installed.packages()[,1])
    
    for (package_name in package_names)
    {
      if (!is_installed(package_name))
      {
        install.packages(package_name,repos)
      }
      library(
        package_name,character.only = TRUE,quietly = TRUE,verbose = FALSE
      )
    }
  }

#' Convert \code{NULL} to zero.
#'
#' \code{null_as_zero} make \code{NULL} to be zero.
#'
#' @param x A number (usually a member of a list) that might be \code{NULL}
#' @return A number
#'
#' @export
#'
#'
null_as_zero <- function(x) {
  if (is.null(x)) {
    return(0)
  } else {
    x
  }
}

#' Stratification setup by covaraites
#'
#' \code{set_strat} makes group indicators based on \code{model_options$X_reg_*}
#'
#' @details the resuls from this function will help stratify etiology or FPR for
#' different strata; the ways of stratification for etiology and FPR can be based
#' on different covariates.
#'
#' @param X   A data frame of covariates
#' @param X_reg The vector of covariates that will stratify the analyses. These
#' variables have to be categorical.
#'
#' @return A list with following elements:
#' \itemize{
#' \item \code{N_group} The number of groups
#' \item \code{group} A vector of group indicator for every observation
#' }
#'
#' @export

set_strat <- function(X,X_reg) {
  if (!is.data.frame(X)) {
    stop("==X is not a data frame. Please transform it into a data frame.==")
  }
  if (!all(X_reg %in% names(X))) {
    stop("==",paste(X_reg,collapse = ", ")," not in X ==")
  }
  X_group        <- X[,X_reg,drop = FALSE]
  #   # dichotomize age variable:
  #   if ("AGECAT" %in% model_options$X_reg_Eti){
  #    X_group$AGECAT <- as.numeric(X_group$AGECAT > 1)+1
  #   }
  X_group$group_names <- apply(X_group,1,paste,collapse = "&")
  X_group$ID          <- 1:nrow(X_group)
  
  form_agg <- as.formula(paste0("cbind(group_names,ID)~",
                                paste(X_reg,collapse = "+")))
  
  grouping <- aggregate(form_agg,X_group,identity)
  
  ## temporary code to get the count of observations in each group:
  #form_agg2 <- as.formula(paste0("cbind(group_names,ID)~",
  #                              paste(c("Y",model_options$X_reg),collapse="+")))
  #aggregate(form_agg2,X_group,length)
  
  group_nm <-
    lapply(grouping$group_names,function(v)
      unique(as.character(v)))
  names(group_nm) <- 1:length(group_nm)
  X_group$grp <- rep(NA,nrow(X_group))
  
  for (l in seq_along(group_nm)) {
    X_group$grp[unfactor(grouping$ID[[l]])] <- l
  }
  
  group <- X_group$grp
  N_vec <- table(group)
  N_grp <- length(N_vec)
  
  list(N_grp = N_grp, group = group)
}


#' Check if covariates are discrete
#'
#' \code{is_discrete} checks if the specified covariates could be regarded as discrete
#' variables.
#'
#' @details Note that this function should be used with caution. It used
#' \deqn{nrow(X)/nrow(unique(X[,X_reg,drop=FALSE]))>10} as an \emph{ad hoc} criterion.
#' It is not the same as \code{\link{is.discrete}}
#'
#' @inheritParams set_strat
#'
#' @return \code{TRUE} for discrete; \code{FALSE} otherwise.
#' @export

is_discrete <- function(X,X_reg) {
  if (!is.data.frame(X)) {
    stop("==X is not a data frame. Please transform it into a data frame.==")
  }
  if (!all(X_reg %in% names(X))) {
    stop("==",paste(X_reg,collapse = ", ")," not in X ==")
  }
  nrow(X) / nrow(unique(X[,X_reg,drop = FALSE])) > 10
}

# is_discrete(X,"ENRLDATE")
# is_discrete(X,c("ENRLDATE","AGECAT"))
# is_discrete(X,c("HIV","AGECAT"))
# is_discrete(X,c("newSITE","AGECAT"))

#' Test for 'try-error' class
#'
#' @param x An object to be test if it is "try-error"
#'
#'
#' @references  \url{http://adv-r.had.co.nz/Exceptions-Debugging.html}
#' @return Logical. \code{TRUE} for "try-error"; \code{FALSE} otherwise
#' @export
is.error <- function(x)
  inherits(x, "try-error")



#' Get unique month from Date
#'
#' \code{unique_month} converts observed dates into unique months
#' to help visualize sampled months
#'
#' @param Rdate standard date format in R
#'
#' @return a vector of characters with \code{month-year}, e.g., \code{4-2012}.
#' @export
#'
#'
#'
unique_month <- function(Rdate) {
  unique(paste(lubridate::month(Rdate),lubridate::year(Rdate),sep = "-"))
}


#' get symmetric difference of months from two vector of R-format dates
#'
#' \code{sym_diff_month} evaluates the symmetric difference between two sets
#' of R-formated date
#'
#' @param Rdate1,Rdate2 R-formated R dates. See \code{\link{as.Date}}
#'
#' @return \code{NULL} if no difference; the set of different months otherwise.
#'
#' @export
#'
sym_diff_month <- function(Rdate1, Rdate2) {
  month1 <- unique_month(Rdate1)
  month2 <- unique_month(Rdate1)
  
  symdiff <- function(x, y) {
    setdiff(union(x, y), intersect(x, y))
  }
  
  res <- symdiff(month1,month2)
  
  if (length(res) == 0) {
    return(NULL)
  }
  
  #cat("==Different months for cases and controls: ", sym_diff_month ,"==")
  return(res)
}



#' Make FPR design matrix for dates with R format.
#'
#' \code{dm_Rdate_FPR} creates desigm matrices for false positive rate regressions.
#'
#' @param Rdate a vector of dates of R format
#' @param Y binary case/control status; 1 for case; 0 for controls
#' @param effect The design matrix for "random" or "fixed" effect; Default
#' is "fixed". When specified as "fixed", it produces standardized R-format dates
#' using control's mean and standard deviation; When specified as "random", it produces
#' \code{num_knots_FPR} columns of design matrix for thin-plate regression splines (TPRS) fitting.
#' One needs both "fixed" and "random" in a FPR regression formula in \code{model_options}
#' to enable TPRS fitting. For example, \code{model_options$X_reg_FPR} can be \cr
#' \cr
#' \code{~ AGECAT+HIV+dm_Rdate_FPR(ENRLDATE,Y,"fixed")+dm_Rdate_FPR(ENRLDATE,Y,"random",10)}\cr
#' \cr
#' means FPR regression with intercept, main effects for 'AGECAT' and 'HIV', and TPRS
#' bases for 'ENRLDATE' using 10 knots placed at 10 equal-probability-spaced sample quantiles.
#' @param num_knots_FPR number of knots for FPR regression; default is \code{NULL}
#' to accomodate fixed effect specification.
#'
#' @seealso \code{\link{nplcm}}
#' @return Design matrix for FPR regression:
#' \itemize{
#' \item \code{Z_FPR_ctrl} transformed design matrix for FPR regression for controls
#' \item \code{Z_FPR_case} transformed design matrix for borrowing FPR
#' regression from controls to cases. It is obtained using control-standardation,
#' and square-root the following matrix (\eqn{\Omega}]) with (\eqn{j_1},\eqn{j_2}) element being
#' \deqn{\Omega_{j_1j_2}=\|knots_{j_1}-knots_{j_2}\|^3}.
#' }
#' @export
dm_Rdate_FPR <- function(Rdate,Y,effect = "fixed",num_knots_FPR = NULL) {
  if (is.null(num_knots_FPR) & effect == "random") {
    stop(
      "==Please specify number of knots for FPR in thin-plate regression spline using 'num_knots_FPR'.=="
    )
  }
  
  if (!is.null(num_knots_FPR) & effect == "fixed") {
    stop("==Don't need 'num_knots_FPR' for fixed effects.==")
  }# standardization:
  df    <- data.frame(Y = Y,num_date = as.numeric(Rdate))
  grp_mean <- c(mean(df$num_date[Y == 0]),mean(df$num_date[Y == 1]))
  grp_sd <- c(sd(df$num_date[Y == 0]),sd(df$num_date[Y == 1]))
  df$ingrp_std_num_date <-
    (df$num_date - grp_mean[df$Y + 1]) / grp_sd[df$Y + 1]
  #outgrp_std_num_date standardizes the cases' dates using controls' mean and sd:
  df$outgrp_std_num_date <- (df$num_date - grp_mean[1]) / grp_sd[1]
  df$outgrp_std_num_date[df$Y == 0] <- NA
  
  if (effect == "random") {
    # for FPR regression in controls:
    ctrl_ingrp_std_num_date <- df$ingrp_std_num_date[df$Y == 0]
    knots_FPR <- quantile(unique(ctrl_ingrp_std_num_date),
                          seq(0,1,length = (num_knots_FPR + 2))[-c(1,(num_knots_FPR +
                                                                        2))])
    
    
    Z_K_FPR_ctrl <-
      (abs(outer(
        ctrl_ingrp_std_num_date,knots_FPR,"-"
      ))) ^ 3
    OMEGA_all_ctrl <- (abs(outer(knots_FPR,knots_FPR,"-"))) ^ 3
    svd.OMEGA_all_ctrl <-  svd(OMEGA_all_ctrl)
    sqrt.OMEGA_all_ctrl <-
      t(svd.OMEGA_all_ctrl$v %*% (t(svd.OMEGA_all_ctrl$u) * sqrt(svd.OMEGA_all_ctrl$d)))
    
    Z_FPR_ctrl <-  t(solve(sqrt.OMEGA_all_ctrl,t(Z_K_FPR_ctrl)))
    
    # for borrowing FPR regression from controls to cases:
    case_outgrp_std_num_date <- df$outgrp_std_num_date[df$Y == 1]
    Z_K_FPR_case <-
      (abs(outer(
        case_outgrp_std_num_date,knots_FPR,"-"
      ))) ^ 3
    Z_FPR_case <- t(solve(sqrt.OMEGA_all_ctrl,t(Z_K_FPR_case)))
    
    ind <- which(Y == 1)
    res <- matrix(NA,nrow = length(Y),ncol = ncol(Z_FPR_case))
    res[ind,] <- Z_FPR_case
    res[-ind,] <- Z_FPR_ctrl
    return(res)
  }
  
  if (effect == "fixed") {
    ind <- which(Y == 1)
    res <- matrix(NA,nrow = length(Y),ncol = 1)
    res[ind,1] <- df$outgrp_std_num_date[ind]
    res[-ind,1] <- df$ingrp_std_num_date[-ind]
    return(res)
  }
}



#' Make etiology design matrix for dates with R format.
#'
#' \code{dm_Rdate_Eti} creates desigm matrices for etiology regressions.
#'
#' @param Rdate a vector of dates of R format
#' @param Y binary case/control status; 1 for case; 0 for controls
#' @param num_knots_Eti number of knots for etiology regression
#' @param basis_Eti the type of basis functions to use for etiology regression. It can be "ncs" (natural
#' cubic splines) or "tprs" (thin-plate regression splines). Default is "ncs". "tprs"
#' will be implemented later.
#'
#' @details It is used in \code{model_options$X_reg_Eti}. For example, one can specify
#' it as: \cr
#' \cr
#' \code{~ AGECAT+HIV+dm_Rdate_Eti(ENRLDATE,Y,5)} \cr
#' \cr
#' to call an etiology regression with intercept, main effects for 'AGECAT' and 'HIV', and
#' natual cubic spline bases for 'ENRLDATE' using 5 knots defined as 5 equal-probability-spaced
#' sample quantiles.
#'
#' @seealso \code{\link{nplcm}}
#' @return Design matrix for etiology regression:
#' \itemize{
#' \item \code{Z_Eti} transformed design matrix for etiology regression
#' }
#' @export
dm_Rdate_Eti <- function(Rdate,Y,num_knots_Eti,basis_Eti = "ncs") {
  #       #
  #       # test:
  #       #
  #       Rdate <- data_nplcm$X$ENRLDATE
  #       Y <- data_nplcm$Y
  #       num_knots_Eti <- 5
  #       basis_Eti = "ncs"
  #       #
  #       #
  #       #
  #
  # standardization:
  df    <- data.frame(Y = Y,num_date = as.numeric(Rdate))
  grp_mean <- c(mean(df$num_date[Y == 0]),mean(df$num_date[Y == 1]))
  grp_sd <- c(sd(df$num_date[Y == 0]),sd(df$num_date[Y == 1]))
  df$ingrp_std_num_date <-
    (df$num_date - grp_mean[df$Y + 1]) / grp_sd[df$Y + 1]
  #outgrp_std_num_date standardizes the cases' dates using controls' mean and sd:
  df$outgrp_std_num_date <- (df$num_date - grp_mean[1]) / grp_sd[1]
  df$outgrp_std_num_date[df$Y == 0] <- NA
  
  if (basis_Eti == "ncs") {
    # for etiology regression:
    case_ingrp_std_num_date <- df$ingrp_std_num_date[df$Y == 1]
    Z_Eti <- splines::ns(case_ingrp_std_num_date,df = num_knots_Eti)
  }
  
  if (basis_Eti == "tprs") {
    stop("==Under development. Please contact maintainer. Thanks.==")
    # for etiology regression:
    case_ingrp_std_num_date <- df$ingrp_std_num_date[df$Y == 1]
    knots_Eti <- quantile(unique(case_ingrp_std_num_date),
                          seq(0,1,length = (num_knots_Eti + 2))[-c(1,(num_knots_Eti +
                                                                        2))])
    
    Z_K_Eti <- (abs(outer(
      case_ingrp_std_num_date,knots_Eti,"-"
    ))) ^ 3
    OMEGA_all_Eti <- (abs(outer(knots_Eti,knots_Eti,"-"))) ^ 3
    svd.OMEGA_all_Eti <- svd(OMEGA_all_Eti)
    sqrt.OMEGA_all_Eti <-
      t(svd.OMEGA_all_Eti$v %*% (t(svd.OMEGA_all_Eti$u) * sqrt(svd.OMEGA_all_Eti$d)))
    Z_Eti  <- t(solve(sqrt.OMEGA_all_Eti,t(Z_K_Eti)))
  }
  
  ind <- which(Y == 1)
  res <- matrix(NA,nrow = length(Y),ncol = ncol(Z_Eti))
  res[ind,]  <- Z_Eti
  res
}

#' create regressor summation equation used in regression for FPR
#'
#' \code{create_bugs_regressor_FPR} creates linear product of coefficients
#' and a row of design matrix used in regression
#'
#' @param n the length of coefficients
#' @param dm_nm name of design matrix; default \code{"dm_FPR"}
#' @param b_nm name of the coefficients; defaul \code{"b"}
#' @param ind_nm name of the coefficient iterator; default \code{"j"}
#' @param sub_ind_nm name of the subject iterator; default \code{"k"}
#'
#' @return a character string with linear product form
#'
#' @export
#'

create_bugs_regressor_FPR <- function(n,dm_nm = "dm_FPR",
                                      b_nm = "b",ind_nm = "j",
                                      sub_ind_nm = "k") {
  summand  <- rep(NA,n)
  for (i in 1:n) {
    summand[i] <- paste0(b_nm,"[",ind_nm,",",i,"]*",
                         dm_nm,"[",sub_ind_nm,",",i + 2,"]")
  }
  paste(summand,collapse = "+")
}

#' create regressor summation equation used in regression for etiology
#'
#' \code{create_bugs_regressor_Eti} creates linear product of coefficients
#' and a row of design matrix used in regression
#'
#' @param n the length of coefficients
#' @param dm_nm name of design matrix; default \code{"dm_Eti"}
#' @param b_nm name of the coefficients; defaul \code{"betaEti"}
#' @param ind_nm name of the coefficient iterator; default \code{"j"}
#' @param sub_ind_nm name of the subject iterator; default \code{"k"}
#'
#' @return a character string with linear product form
#'
#' @export
#'
create_bugs_regressor_Eti <- function(n,dm_nm = "dm_Eti",
                                      b_nm = "betaEti",ind_nm = "j",
                                      sub_ind_nm = "k") {
  summand  <- rep(NA,n)
  for (i in 1:n) {
    summand[i] <- paste0(b_nm,"[",ind_nm,",",i,"]*",
                         dm_nm,"[",sub_ind_nm,",",i,"]")
  }
  paste(summand,collapse = "+")
}


#' Deletes a pattern from the start of a string, or each of a vector of strings.
#'
#' \code{delete_start_with} is used for clean the column names in raw data.
#' For example, R adds "X" at the start of variable names. This function deletes
#' "X_"s from the column names. This can happen if the raw data have column
#' names such as "\code{_CASE_ABX}".
#'
#' @param s the pattern (a single string) to be deleted from the start.
#' @param vec a vector of strings with unwanted starting strings (specified by \code{s}).
#'
#' @return string(s) with deleted patterns from the start.
#'
#' @examples
#' delete_start_with("X_",c("X_hello"))
#' delete_start_with("X_",c("X_hello","hello2"))
#' delete_start_with("X_",c("X_hello","hello2","X_hello3"))
#' @export


delete_start_with = function(s,vec) {
  ind = grep(s,substring(vec,1,nchar(s)))
  old = vec[ind]
  vec[ind] = substring(old,nchar(s) + 1)
  return(vec)
}



#' Takes any number of R objects as arguments and returns a list whose names are
#' derived from the names of the R objects.
#'
#' Roger Peng's listlabeling challenge from
#' \url{http://simplystatistics.tumblr.com/post/11988685443/computing-on-the-language}.
#' Code copied from \url{https://gist.github.com/ajdamico/1329117/0134148987859856fcecbe4446cfd37e500e4272}
#'
#' @param ... any R objects
#'
#' @return a list as described above
#'
#' @examples
#' #create three example variables for a list
#' x <- 1
#' y <- 2
#' z <- "hello"
#' #display the results
#' make_list( x , y , z )
#' @export

#create the function
make_list <- function(...) {
  #put all values into a list
  argument_values <- list(...)
  
  #save all argument names into another list
  argument_names <- as.list(sys.call())
  
  #cycle through the first list and label with the second, ignoring the function itself
  for (i in 2:length(argument_names)) {
    names(argument_values)[i - 1] <- argument_names[i]
  }
  
  #return the newly-labeled function
  argument_values
}

#' Make a list with numbered names
#'
#' To collect multiple measurements within the same category, e.g., bronze-standard.
#'
#' @param ... any R object
#'
#' @return a list with names numbered
#'
#' @export

make_numbered_list <- function(...) {
  #put all values into a list
  argument_values <- list(...)
  
  #save all argument names into another list
  argument_names <- as.list(sys.call())
  
  #cycle through the first list and label with the second, ignoring the function itself
  for (i in 2:length(argument_names)) {
    names(argument_values)[i - 1] <- paste0(argument_names[i],"_",i - 1)
  }
  
  #return the newly-labeled function
  argument_values
}


#' make a mapping template for model fitting
#'
#' \code{make_template} creates a mapping matrix so that a measurement is mapped
#' to inform a particular latent status. Crucial for model fitting.
#'
#' @details The first argument has to be character substrings from the second argument. The
#' second argument can have character strings not matched in the first argument.
#' 
#' @param patho a vector of pathogen names. \code{patho}
#'  must be substring of some \code{cause_list} elements, e.g.,
#'  "PNEU" is a substring of "PNEU_VT13". Also see examples.
#' @param cause_list the list of potential latent status
#'
#' @examples
#'
#' cause_list <- c("HINF","PNEU_VT13","PNEU_NOVT13","SAUR","HMPV_A_B","FLU_A",
#' "PARA_1","PARA_3","PARA_4","PV_EV","RHINO","RSV", "ENTRB","TB")
#'
#' patho_BrS_NPPCR <- c("HINF","PNEU","SAUR","HMPV_A_B","FLU_A","PARA_1",
#' "PARA_3","PARA_4","PV_EV","RHINO","RSV")
#' make_template(patho_BrS_NPPCR,cause_list)
#'
#'
#'  cause = c("A","B1","B2","C","A+C","B+C")
#'  patho = c("A","B","C")
#'  make_template(patho,cause)
#'  
#'  cause = c("A","B1","B2","C","A+C","B+C","other")
#'  patho = c("A","B","C")
#'  make_template(patho,cause)
#' @return a mapping from \code{patho} to \code{cause_list}.
#'\code{NROW = length(cause_list)+1};
#'\code{NCOL = length(patho)}. This value is crucial in model fitting to determine
#'which measurements are informative of a particular category of latent status.
#'
#' @export

make_template <- function(patho, cause_list) {
  # patho must be substring of some cause_list elements, e.g., "PNEU" is a substring of "PNEU_VT13".
  res <- list()
  for (i in seq_along(patho)) {
    res[[i]] <- as.numeric(grepl(patho[i],cause_list))
  }
  template <- t(do.call(rbind,res))
  
  template <- rbind(template,rep(0,ncol(template)))
  
  # give names to columns and rows:
  nm_row <- c(cause_list,"control")
  nm_col <- patho
  
  rownames(template) <- nm_row
  colnames(template) <- nm_col
  
  template
}

#' get position to store in data_nplcm$Mobs:
#' 
#' @details  also works for a vector
#' 
#' @param quality_nm names of quality: can be "BrS", "SS" or "GS"
#' 
#' @return position of the quality name: "BrS"-1; "SS"-2; "GS"-3.
#' 
#' @examples
#' \dontrun{
#' lookup_quality("BrS")
#' lookup_quality("HH")
#' }
#' 
#' @export

lookup_quality <- function(quality_nm) {
  res_pos <- list()
  for (i in seq_along(quality_nm)){
    if (quality_nm[i] == "BrS") {
      res_pos[i] <- 1
    }
    if (quality_nm[i] == "SS") {
      res_pos[i] <- 2
    }
    if (quality_nm[i] == "GS") {
      res_pos[i] <- 3
    }
    if (! (quality_nm[i] %in% c("BrS","SS","GS"))){
      stop("== Please provide measurement quality names within (\"BrS\",\"SS\",\"GS\")! ==")
    }
  }
  unlist(res_pos)
}

#' parse regression components (either false positive rate or etiology regression)
#' for fitting nplcm
#' 
#' @param form regression formula
#' @param data_nplcm data object for nplcm; must contain covariates X and outcome Y.
#' 
#' @return TURE for doing regression; FALSE otherwise.
#' 
#' @export
parse_nplcm_reg <- function(form,data_nplcm){
  if (is.null(data_nplcm$X)) {
    print("\n == There are no covariate data in `data_nplcm` =="); 
    return(FALSE)
  }
  
  res <-
    try(model.matrix(form,data.frame(data_nplcm$X,Y = data_nplcm$Y)))
  if (is.error(res)) {
    stop()
  }
  if (is.empty.model(form)) {
    return(FALSE)
  } else{
    return(TRUE)
  }
}


#' convert one column data frame to a vector
#' 
#' @param x an one-column data.frame
#' 
#' @return a vector
#' 
#' @export
as.matrix_or_vec <- function(x){
  if (ncol(x)==1){
    return(c(as.matrix(x)))
  }  
  as.matrix(x)
}
