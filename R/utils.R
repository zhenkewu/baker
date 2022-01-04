# Taken from https://github.com/suyusung/R2jags/blob/master/R/util.R#L12
repath <- function(x) {
  n.part <- length(x[[1]])
  temp <- rep(NA, n.part)
  for (i in 1:n.part){
    if (nchar(x[[1]][i])> 8){
      temp[i] <- paste(substr(x[[1]][i], 1, 6), "~1/", sep="")
    }
    if (nchar(x[[1]][i])<=8){
      temp[i] <- paste(x[[1]][i],"/", sep="")
    }
  }
  return(temp)
}


# Taken from https://github.com/suyusung/R2jags/blob/master/R/util.R#L27
win2unixdir <- function(windir){
  Dir <- substr(windir, 1, 3)
  tempdir <- substr(windir, 4, 1000)
  tempdir <- gsub(" ", "", tempdir)
  tempdir <- strsplit(tempdir, "/")
  n.part <- length(tempdir[[1]])
  if (n.part>0){
    tempdir <- repath(tempdir)
    path <- ""
    for (i in 1:n.part){
      path <- paste(path, tempdir[i], sep="")
    }
    newpath <- paste(Dir, path, sep="")
  }
  else {
    newpath <- Dir
  }
  return(c(newpath))
}

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
#' (same length with `"p"`).
#' The Bernoulli random variables can have different means.
#'
#' @param p A vector of probabilities, each being the head probability
#' of an independent coin toss
#'
#' @return A vector of 1s (head) and 0s (tail)
#' @export
rvbern <-
  function(p) {
    U  <- stats::runif(length(p),0,1)
    res <- (U < p) + 0
    res
  }



#' Convert a 0/1 binary-coded sequence into decimal digits
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
#' Reverse to [symb2I()]
#' @param binary_code Binary indicators for pathogens
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



#' Convert a matrix of binary indicators to categorical variables
#'
#' @param binary_mat The matrix of binary indicators. Rows for subjects, columns for pathogens in the `"pathogen.list"`
#' @param cause_list The list of causes
#' @param pathogen_list The complete list of pathogen names
#'
#' @return A vector of categorical variables. Its length equals the length of `"allowed.list"`
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
    stop("==Some binary pattern in 'binary_mat' is not included by 'cause_list'!==")
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
  y = stats::dbeta(x,a,b)
  graphics::plot(x,y,type = "l",main = paste0("a=",a,",b=",b))
}


#' Get package from CRAN website
#'
#' @param pckg package name
#'
#' @return None
#' @export
getPckg <- function(pckg) {
  utils::install.packages(pckg, repos = "http://cran.r-project.org")
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
#' Use + for ascending, - for descending.\cr
#' Sorting is left to right in the formula\cr
#' Usage is either of the following:\cr
#' sort.data.frame(~Block-Variety,Oats)\cr
#' sort.data.frame(Oats,~-Variety+Block)\cr
#'
#'
#'@param form A Formula. See example in **Description**.
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
  # Usage is either of the following:
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
#' @return A permuted vector from 1 to `length(raw_nm)`. For example, if
#' its first element is 3, it means that the 3rd pathogen in `raw_nm`
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
      
      fit <- stats::glm(y ~ x,family = stats::binomial(link = "logit"))
      
      if ("x" %in% rownames(summary(fit)$coef)) {
        logORmat[j2,j1] = round(summary(fit)$coef["x",1],3)
        logORmat.se[j2,j1] = round(summary(fit)$coef["x",2],3)
      }
      # controls: (lower triangle)
      x <- MBS.ctrl[,j2]
      y <- MBS.ctrl[,j1]
      
      fit <- stats::glm(y ~ x,family = stats::binomial(link = "logit"))
      
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
#' names of the rows, from bottom to top. Default is 1 through `ncol(mat)`;
#' @param cell_metrics the meaning of number in every cell;
#' @param folding_line Default is `TRUE` for adding dashed major diagonal
#' line.
#' @param axes plot axes; default is `FALSE`;
#' @param xlab label for x-axis
#' @param ylab label for y-axis
#' @param asp aspect ratio; default is `1` to ensure square shape
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
  
  graphics::par(mar = c(0, 0, 5, 0), bg = "white",xpd = TRUE)
  graphics::plot(
    c(0, n + 0.8), c(0, n + 0.8), axes = axes, xlab = "",
    ylab = "", asp = 1, type = "n"
  )
  ##add grid
  graphics::segments(rep(0.5, n + 1), 0.5 + 0:n, rep(n + 0.5, n + 1),
                     0.5 + 0:n, col = "gray")
  graphics::segments(0.5 + 0:n, rep(0.5, n + 1), 0.5 + 0:n, rep(n + 0.5,
                                                                n), col = "gray")
  mat.txt <- round(t(mat)[,n:1],1)
  mat.txt3 <- round(t(mat)[,n:1],3)
  
  # add meaning of the number in a cell:
  graphics::text(1,J,cell_metrics,cex = cex_main / 2)
  
  for (i in 1:n) {
    for (j in 1:n) {
      abs.mat <- abs(mat.txt[i,j])
      if ((!is.na(abs.mat)) && abs.mat > 2) {
        graphics::text(
          i,j,round(mat.txt3[i,j],1),
          col = ifelse(mat.txt3[i,j] > 0,"red","blue"),cex = cex_main
        )
      }
      
    }
  }
  
  if (folding_line) {
    # diagonal line:
    graphics::segments(
      0.5 + 1,.5 + n - 1,.5 + n,0.5,col = "black",lty = 3,lwd = 3
    )
  }
  
  # put pathogen names on rows and columns:
  for (s in 1:J) {
    
    graphics::text(0.25,J-s+1,paste0(dim_names[s],":(",s,")"),cex=min(1.5,20/J),srt=45,adj=1)
    graphics::text(
      s,J + 0.7,paste0("(",s,"):",dim_names[s]),cex = min(1.5,20 / J),srt = 45,adj =
        0
    )
  }
  # labels for cases and controls:
  graphics::text(J + 1,J / 2,"cases",cex = 2,srt = -90)
  graphics::text(J / 2,0,"controls",cex = 2)
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
#' `beta_parms_from_quantiles` produces prior Beta parameters for
#'  the true positive rates (TPR)
#'
#' @param q A vector of lower and upper bounds, in which Beta distribution
#' will have quantiles specified by `p`. For example, `q=c(0.5,0.99)`
#' @param p The lower and upper quantiles of the range one wants to specify.
#' @param precision Approximation precisions.
#' @param derivative.epsilon Precision of calculating derivative.
#' @param start.with.normal.approx Default is `TRUE`, for normal approximation.
#' @param start Starting values of beta parameters.
#' @param plot Default is `FALSE` to suppress plotting of the beta density,
#' otherwise, set to `TRUE`.
#'
#' @return A list containing the selected Beta parameters `a`, and `b`.
#' Other elements of the list include some details about the computations involved
#' in finding `a` and `b`.
#'
#' @references <http://www.medicine.mcgill.ca/epidemiology/Joseph/PBelisle/BetaParmsFromQuantiles.html>
#' @export
#'
beta_parms_from_quantiles <- function(q, p = c(0.025,0.975),
                                      precision = 0.001,
                                      derivative.epsilon = 1e-3,
                                      start.with.normal.approx = TRUE,
                                      start = c(1, 1),
                                      plot = FALSE) {
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
      stats::dbeta(x, shape1 = theta[1], shape2 = theta[2])
    }
  F.inv <-
    function(x, theta) {
      stats::qbeta(x, shape1 = theta[1], shape2 = theta[2])
    }
  f.cum <-
    function(x, theta) {
      stats::pbeta(x, shape1 = theta[1], shape2 = theta[2])
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
  
  dens.label <- 'stats::dbeta'
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
    par.usr <- graphics::par('usr')
    par.din <- graphics::par('din')
    
    p.string <-
      as.character(round(c(0,1) + c(1,-1) * p.check, digits = 4))
    str.width <- graphics::strwidth(p.string, cex = cex)
    str.height <- graphics::strheight("0", cex = cex)
    
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
      graphics::text(
        x, y, labels = p.string[1], adj = c(1,0), col = 'gray', cex = cex
      )
    }
    else
    {
      graphics::text(
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
      graphics::text(
        x, y, labels = p.string[2], adj = c(0,0), col = 'gray', cex = cex
      )
    }
    else
    {
      graphics::text(
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
      v <- (diff(q) / diff(stats::qnorm(p))) ^ 2
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
    
    if (abs(theta[1]-1)<0.0001 && abs(theta[2]-1)<0.0001){
      graphics::plot(x, h, type = 'l', ylab = ylab,ylim=c(0,2))
    }else{
      graphics::plot(x, h, type = 'l', ylab = ylab)
    }
    
    # fill in area on the left side of the distribution
    x <- seq(from = plot.xlim[1], to = q[1], length = 1000)
    y <- f(x, theta)
    x <- c(x, q[1], plot.xlim[1]); y <- c(y, 0, 0)
    graphics::polygon(x, y, col = 'lightgrey', border = 'lightgray')
    # fill in area on the right side of the distribution
    x <- seq(from = max(plot.xlim), to = q[2], length = 1000)
    y <- f(x, theta)
    x <- c(x, q[2], plot.xlim[2]); y <- c(y, 0, 0)
    graphics::polygon(x, y, col = 'lightgrey', border = 'lightgray')
    # draw distrn again
    graphics::points(x0, f0, type = 'l')
    h <- f(q, theta)
    graphics::points(rep(q[1], 2), c(0, h[1]), type = 'l', col = 'orange')
    graphics::points(rep(q[2], 2), c(0, h[2]), type = 'l', col = 'orange')
    # place text on both ends areas
    print.area.text(p, p.check, q, f, f.cum, F.inv, theta, mode, cex, plot.xlim)
    
    xaxp <- graphics::par("xaxp")
    x.ticks <- seq(from = xaxp[1], to = xaxp[2], length = xaxp[3] + 1)
    q2print <-
      as.double(setdiff(as.character(q), as.character(x.ticks)))
    
    graphics::mtext(
      q2print, side = 1, col = 'orange', at = q2print, cex = 0.6, line = 2.1
    )
    graphics::points(q, rep(graphics::par('usr')[3] + 0.15 * graphics::par('cxy')[2], 2), pch = 17, col =
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
      parms$theta,plot.xlim, dens.label, parms.names, 0.8
    )
  
  list(
    a = parms$theta[1], b = parms$theta[2], last.change = parms$last.change,
    niter = parms$niter, q = q, p = p, p.check = p.check
  )
}







#' Load or install a package
#'
#'`load_or_install` checks if a package is installed,
#' loads it if it is, and installs it if not.
#'
#' @param package_names A vector of package names
#' @param repos URL for downloading the packages. 
#'
#' @references Credit:
#'  <http://www.vikparuchuri.com/blog/loading-andor-installing-packages/>
#' @return No message if it successfully loads the specified packages; Error
#'  if such a package does not exist.
#'
#' @export
#'
#'
load_or_install <-
  function(package_names,repos) {
    is_installed <-
      function(mypkg)
        is.element(mypkg, utils::installed.packages()[,1])
    
    for (package_name in package_names)
    {
      if (!is_installed(package_name))
      {
        utils::install.packages(package_name,repos)
      }
      library(
        package_name,character.only = TRUE,quietly = TRUE,verbose = FALSE
      )
    }
  }

#' Convert `NULL` to zero.
#'
#' `null_as_zero` make `NULL` to be zero.
#'
#' @param x A number (usually a member of a list) that might be `NULL`
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

#' Stratification setup by covariates
#'
#' `set_strat` makes group indicators based on `model_options$X_reg_*`
#'
#' @details the results from this function will help stratify etiology or FPR for
#' different strata; the ways of stratification for etiology and FPR can be based
#' on different covariates.
#'
#' @param X   A data frame of covariates
#' @param X_reg The vector of covariates that will stratify the analyses. These
#' variables have to be categorical.
#'
#' @return A list with following elements:
#' \itemize{
#' \item `N_group` The number of groups
#' \item `group` A vector of group indicator for every observation
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
  
  form_agg <- stats::as.formula(paste0("cbind(group_names,ID)~",
                                       paste(X_reg,collapse = "+")))
  
  grouping <- stats::aggregate(form_agg,X_group,identity)
  
  ## temporary code to get the count of observations in each group:
  #form_agg2 <- stats::as.formula(paste0("cbind(group_names,ID)~",
  #                              paste(c("Y",model_options$X_reg),collapse="+")))
  #stats::aggregate(form_agg2,X_group,length)
  
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
#' `is_discrete` checks if the specified covariates could be regarded as discrete
#' variables.
#'
#' @details Note that this function should be used with caution. It used
#' \deqn{nrow(X)/nrow(unique(X[,X_reg,drop=FALSE]))>10} as an *ad hoc* criterion.
#' It is not the same as [plyr::is.discrete()]
#'
#' @param X   A data frame of covariates
#' @param X_reg The vector of covariates that will stratify the analyses. These
#' variables have to be categorical. Or a formula (can be tested by `plyr::is.formula`), 
#' e.g., `~as.factor(SITE8) + as.factor(AGECAT > 1)`.
#'
#' @return `TRUE` for all being discrete; `FALSE` otherwise.
#' @export

is_discrete <- function(X,X_reg) {
  if (!is.data.frame(X)) {
    stop("==[baker]X is not a data frame. Please transform it into a data frame.==")
  }
  
  if (is.character(X_reg)){
    if (!all(X_reg %in% names(X))) {
      stop("==[baker]",paste(X_reg,collapse = ", ")," not in the X data provided ==")
    }
    res <- nrow(X) / nrow(unique(X[,X_reg,drop = FALSE])) > 10
    
  } else{
    X_dm <- stats::model.matrix(X_reg,data.frame(X)) # <--- X for case only.
    res  <- nrow(X) / nrow(unique(X_dm)) > 10
  }
  res
}

# is_discrete(X,"ENRLDATE")
# is_discrete(X,c("ENRLDATE","AGECAT"))
# is_discrete(X,c("HIV","AGECAT"))
# is_discrete(X,c("newSITE","AGECAT"))

# or is_discrete(X,~as.factor(SITE8) + as.factor(AGECAT > 1))

#' Test for 'try-error' class
#'
#' @param x An object to be test if it is "try-error"
#'
#'
#' @references  <http://adv-r.had.co.nz/Exceptions-Debugging.html>
#' @return Logical. `TRUE` for "try-error"; `FALSE` otherwise
#' @export
is.error <- function(x)
  inherits(x, "try-error")



#' Get unique month from Date
#'
#' `unique_month` converts observed dates into unique months
#' to help visualize sampled months
#'
#' @param Rdate standard date format in R
#'
#' @return a vector of characters with `month-year`, e.g., `4-2012`.
#' @export
#'
#'
#'
unique_month <- function(Rdate) {
  unique(paste(lubridate::month(Rdate),lubridate::year(Rdate),sep = "-"))
}


#' get symmetric difference of months from two vector of R-format dates
#'
#' `sym_diff_month` evaluates the symmetric difference between two sets
#' of R-formatted date
#'
#' @param Rdate1,Rdate2 R-formatted R dates. See [as.Date()]
#'
#' @return `NULL` if no difference; the set of different months otherwise.
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
#' `dm_Rdate_FPR` creates design matrices for false positive rate regressions; 
#' can also be used to standardize dates.
#'
#' @param Rdate a vector of dates of R format
#' @param Y binary case/control status; 1 for case; 0 for controls
#' @param effect The design matrix for "random" or "fixed" effect; Default
#' is "fixed". When specified as "fixed", it produces standardized R-format dates
#' using control's mean and standard deviation; When specified as "random", it produces
#' `num_knots_FPR` columns of design matrix for thin-plate regression splines (TPRS) fitting.
#' One needs both "fixed" and "random" in a FPR regression formula in `model_options`
#' to enable TPRS fitting. For example, `model_options$likelihood$FPR_formula` can be \cr
#' \cr
#' `~ AGECAT+HIV+dm_Rdate_FPR(ENRLDATE,Y,"fixed")+dm_Rdate_FPR(ENRLDATE,Y,"random",10)`\cr
#' \cr
#' means FPR regression with intercept, main effects for 'AGECAT' and 'HIV', and TPRS
#' bases for 'ENRLDATE' using 10 knots placed at 10 equal-probability-spaced sample quantiles.
#' @param num_knots_FPR number of knots for FPR regression; default is `NULL`
#' to accommodate fixed effect specification.
#'
#' @seealso [nplcm()]
#' @return Design matrix for FPR regression:
#' \itemize{
#' \item `Z_FPR_ctrl` transformed design matrix for FPR regression for controls
#' \item `Z_FPR_case` transformed design matrix for borrowing FPR
#' regression from controls to cases. It is obtained using control-standardization,
#' and square-root the following matrix (\eqn{\Omega}]) with (\eqn{j_1},\eqn{j_2}) element being
#' \deqn{\Omega_{j_1j_2}=\|knots_{j_1}-knots_{j_2}\|^3}.
#' }
#' @export
dm_Rdate_FPR <- function(Rdate,Y,effect = "fixed",num_knots_FPR = NULL) {
  if (is.null(num_knots_FPR) & effect == "random") {
    stop(
      "==[baker] Please specify number of knots for FPR in thin-plate regression spline using 'num_knots_FPR'.=="
    )
  }
  
  if (!is.null(num_knots_FPR) & effect == "fixed") {
    stop("==Don't need 'num_knots_FPR' for fixed effects.==")
  }
  
  # standardization:
  df       <- data.frame(Y = Y,num_date = as.numeric(Rdate))
  grp_mean <- c(mean(df$num_date[Y == 0]),mean(df$num_date[Y == 1]))
  grp_sd   <- c(stats::sd(df$num_date[Y == 0]),stats::sd(df$num_date[Y == 1]))
  df$ingrp_std_num_date <-
    (df$num_date - grp_mean[df$Y + 1]) / grp_sd[df$Y + 1]
  #outgrp_std_num_date standardizes the cases' dates using controls' mean and stats::sd:
  df$outgrp_std_num_date <- (df$num_date - grp_mean[1]) / grp_sd[1]
  df$outgrp_std_num_date[df$Y == 0] <- NA
  
  if (effect == "random") {
    # for FPR regression in controls:
    ctrl_ingrp_std_num_date <- df$ingrp_std_num_date[df$Y == 0]
    knots_FPR <- stats::quantile(unique(ctrl_ingrp_std_num_date),
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
  
  if (effect == "fixed" && is.null(num_knots_FPR)) {
    ind <- which(Y == 1)
    res <- matrix(NA,nrow = length(Y),ncol = 1)
    res[ind,1] <- df$outgrp_std_num_date[ind]
    res[-ind,1] <- df$ingrp_std_num_date[-ind]
    return(res)
  }
}



#' Make etiology design matrix for dates with R format.
#'
#' `dm_Rdate_Eti` creates design matrices for etiology regressions.
#'
#' @param Rdate a vector of dates of R format
#' @param Y binary case/control status; 1 for case; 0 for controls
#' @param num_knots_Eti number of knots for etiology regression
#' @param basis_Eti the type of basis functions to use for etiology regression. It can be "ncs" (natural
#' cubic splines) or "tprs" (thin-plate regression splines). Default is "ncs". "tprs"
#' will be implemented later.
#'
#' @details It is used in `model_options$likeihood$Eti_formula`. For example, one can specify
#' it as: \cr
#' \cr
#' `~ AGECAT+HIV+dm_Rdate_Eti(ENRLDATE,Y,5)` \cr
#' \cr
#' to call an etiology regression with intercept, main effects for 'AGECAT' and 'HIV', and
#' natural cubic spline bases for 'ENRLDATE' using 5 knots defined as 5 equal-probability-spaced
#' sample quantiles.
#'
#' @seealso [nplcm()]
#' @return Design matrix for etiology regression:
#' \itemize{
#' \item `Z_Eti` transformed design matrix for etiology regression
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
  grp_sd <- c(stats::sd(df$num_date[Y == 0]),stats::sd(df$num_date[Y == 1]))
  df$ingrp_std_num_date <-
    (df$num_date - grp_mean[df$Y + 1]) / grp_sd[df$Y + 1]
  #outgrp_std_num_date standardizes the cases' dates using controls' mean and stats::sd:
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
    knots_Eti <- stats::quantile(unique(case_ingrp_std_num_date),
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
#' `create_bugs_regressor_FPR` creates linear product of coefficients
#' and a row of design matrix used in regression
#'
#' @param n the length of coefficients
#' @param dm_nm name of design matrix; default `"dm_FPR"`
#' @param b_nm name of the coefficients; default `"b"`
#' @param ind_nm name of the coefficient iterator; default `"j"`
#' @param sub_ind_nm name of the subject iterator; default `"k"`
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
#' `create_bugs_regressor_Eti` creates linear product of coefficients
#' and a row of design matrix used in regression
#'
#' @param n the length of coefficients
#' @param dm_nm name of design matrix; default `"dm_Eti"`
#' @param b_nm name of the coefficients; default `"betaEti"`
#' @param ind_nm name of the coefficient iterator; default `"j"`
#' @param sub_ind_nm name of the subject iterator; default `"k"`
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
#' `delete_start_with` is used for clean the column names in raw data.
#' For example, R adds "X" at the start of variable names. This function deletes
#' "X_"s from the column names. This can happen if the raw data have column
#' names such as "`_CASE_ABX`". Check [clean_perch_data()] for 
#' its actual usage.
#'
#' @param s the pattern (a single string) to be deleted from the start.
#' @param vec a vector of strings with unwanted starting strings (specified by `s`).
#'
#' @return string(s) with deleted patterns from the start.
#'
#' @examples
#' delete_start_with("X_",c("X_hello"))
#' delete_start_with("X_",c("X_hello","hello2"))
#' delete_start_with("X_",c("X_hello","hello2","X_hello3"))
#' 
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
#' <http://simplystatistics.tumblr.com/post/11988685443/computing-on-the-language>.
#' Code copied from <https://gist.github.com/ajdamico/1329117/0134148987859856fcecbe4446cfd37e500e4272>
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
#' `make_template` creates a mapping matrix (binary values). Each pathogen 
#' in a measurement slice (e.g., nasal-pharyngeal PCR test) is mapped to inform
#' one category of latent status. All the possible categories (e.g., causes of pneumonia) 
#' remain the same regardless of the measurement slice used (e.g., NPPCR or BCX).
#'
#' @details The first argument has to be character substrings from the second argument. 
#' For example, the two arguments can respectively be `"A"` and `"A_1"`, 
#' or `"A"` and `"A+B"`.The second argument can have character strings not 
#' matched in the first argument. If so, it means some causes of diseases are not 
#' directly measured in the current measurement slice. 
#' For each element of `patho`, the function matches from the start of the strings
#' of `cause_list`. Therefore, make sure that latent statuses from the same family 
#' (e.g., "PNEU_VT13" and "PNEU_NOVT13") need to start with the same family name 
#' (e.g., "PNEU") followed by subcategories (e.g., "_VT13" and "_NOVT13").
#' 
#' @param patho A vector of pathogen names for a particular measurement slice. 
#' `patho` must be a substring of some elements in `cause_list`, e.g.,
#'  "PNEU" is a substring of "PNEU_VT13". Also see Examples for this function.
#'  
#' @param cause_list A vector of characters; Potential categories of latent statuses.
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
#'  
#'  
#'  cause = c("A","B1","B2","X_B","Y_B","C","A+C","B+C","other")
#'  patho = c("A","B","C","X_B","Y_B")
#'  make_template(patho,cause)
#'  
#'  
#' @return a mapping from `patho` to `cause_list`.
#'`NROW = length(cause_list)+1`;
#'`NCOL = length(patho)`. This value is crucial in model fitting to determine
#'which measurements are informative of a particular category of latent status.
#'
#' @export
make_template <- function(patho, cause_list) {
  # patho must be substring of some cause_list elements, e.g., "PNEU" is a substring of "PNEU_VT13".
  res <- list()
  for (i in seq_along(patho)) {
    pat <- eval(paste0("(^|\\+)",patho[i]))
    res[[i]] <- as.numeric(grepl(pat,cause_list))
  }
  template <- t(do.call(rbind,res))
  
  template <- rbind(template,rep(0,ncol(template)))
  
  # give names to columns and rows:
  nm_row <- c(cause_list,"control")
  nm_col <- patho
  
  rownames(template) <- nm_row
  colnames(template) <- nm_col
  
  template <- matrix(as.integer(template),nrow=nrow(template),ncol=ncol(template))
  
  template
}

#' Get position to store in data_nplcm$Mobs:
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
#' @seealso [extract_data_raw()]
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
#' for fitting npLCM; Only use this when formula is not `NULL`.
#' 
#' @param form regression formula
#' @param data_nplcm data object for [nplcm()]; may contain covariates X; 
#' must have case-control status Y.
#' @param silent Default is `TRUE` for no message about covariates; 
#' `FALSE` otherwise.
#' @return `TRUE` for doing regression; `FALSE` otherwise.
#' 
#' @export
parse_nplcm_reg <- function(form,data_nplcm,silent=TRUE){
  if (is.null(data_nplcm$X)) {
    if(!silent){print(" ==[baker] There are no covariate data in `data_nplcm`. ==\n")}; 
    return(FALSE)
  }
  res <-
    try(stats::model.matrix(form,data.frame(data_nplcm$X,Y = data_nplcm$Y)))
  if (is.error(res)) {
    stop("==[baker] Cannot parse regression formula.==\n")
  }
  
  if (stats::is.empty.model(form)) {
    stop("==[baker] An empty regression model is specified. Please replace it by 
         `~ 1` for no regression or `~1+X1+X2` for regression with intercept, X1 and X2 (additive). ==\n")
  } else if (ncol(res)==1 & all(res==1)){
    return(FALSE)
  } else{
    return(TRUE)
  }
}


#' check if the formula is intercept only
#' 
#' outputs logical values for a formula; to identify intercept-only formula.
#' 
#' @param form Regression formula
#' 
#' @return `TRUE` for intercept-only; `FALSE` otherwise
#' 
#' @export
is_intercept_only <- function(form){
  form_remove_space <- gsub(" ","",form,fixed=TRUE)
  formula_parts <- strsplit(form_remove_space[2],"+",fixed=TRUE)[[1]]
  res <- length(formula_parts)==1 && formula_parts=="1"
  return(res)
}



#' convert one column data frame to a vector
#' 
#' @details WinBUGS/JAGS cannot accept a dataframe with one column; This function
#' converts it to a vector, which WinBUGS/JAGS will allow.
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

#' get index of latent status
#' 
#' @param cause_list see mode_options in [nplcm]
#' @param ord order of cause_list according to posterior mean
#' @param select_latent Default is NULL
#' @param exact Default is TRUE
#' 
#' @return a vector of indices
#' @export
get_latent_seq <- function(cause_list, ord,select_latent=NULL,exact=TRUE){
  cause_list_ord <- cause_list[ord]
  latent_seq <- 1:length(cause_list)
  original_num <- ord
  if (!is.null(select_latent)){
    original_index <- sapply(paste("^",select_latent,"$",sep=""),grep,cause_list)
    original_num   <- original_index[my_reorder(cause_list_ord,select_latent)]
    latent_seq <- sapply(paste("^",select_latent,"$",sep=""),grep,cause_list_ord)[my_reorder(cause_list[ord],select_latent)]
    if (!exact){
      original_index <- which(rowSums(make_template(select_latent,cause_list))>0)
      original_num <- original_index[my_reorder(cause_list_ord,cause_list[original_index])]
      latent_seq <- which(rowSums(make_template(select_latent,cause_list_ord))>0)
    }
    
    if (length(latent_seq)!=length(select_latent)){warning("==Some of `select_latent` are not matched by `cause_list` in `model_options$likelihood` ! Please check you supplied correct names in `select_latent`.==")}
  }
  return(make_list(latent_seq,original_num))
}



#' Shannon entropy for multivariate discrete data
#' 
#' @param px a vector of positive numbers sum to 1
#' 
#' @return a non-negative number
#' @export
#' 
H <- function(px) {
  -sum(px * log(px))
}

#' Shannon entropy for binary data
#' 
#' @param m_px a number between 0 and 1
#' 
#' @return a non-negative number
#' 
#' @export
#' 
marg_H <- function(m_px){-m_px*log(m_px)-(1-m_px)*log(1-m_px)}


#' load an object from .RDATA file
#' 
#' @param objName the name of the object
#' @param file the file path
#' @param envir environment; default is calling environment: [parent.frame]
#' @param assign.on.exit default is TRUE
#' 
#' @return a new environment
#' @export
loadOneName <- function(objName, file, envir = parent.frame(),
                        assign.on.exit = TRUE) {
  tempEnv <- new.env()
  load(file, envir = tempEnv)
  stopifnot(objName %in% ls(tempEnv))
  if(assign.on.exit) {
    assign(objName, tempEnv[[objName]], envir = envir)
    return(invisible(tempEnv[[objName]]))
  }
  tempEnv[[objName]]
}


#' generate stick-breaking prior (truncated) from a vector of random probabilities
#' 
#' @param u a vector of probabilities, with the last element 1.
#' 
#' @return a vector of the same length as u; sum to 1.
#' 
#' @examples 
#' 
#' graphics::par(mfrow=c(3,3),oma=c(0,1,5,0),
#'    mar=c(1,2,1,1))
#' for (iter in 1:9){
#'  u   <- c(rbeta(9,1,0.8),1)
#'  res <- tsb(u)
#'  barplot(res,ylim=c(0,1),main=paste0("Random Sample #", iter),ylab="Probability")
#' }
#' graphics::mtext("Truncated Stick-Breaking Dist. (10 segments)",3,
#'      outer=TRUE,cex=1.5,line=1.5)
#' @export
#' 
tsb <- function(u){
  K <- length(u)
  if (u[K]!=1) {stop("==The last element of u must be 1 for truncated stick-breaking!==\n")}
  w <- rep(NA,K)
  w[1] <- u[1]
  for (k in 2:(K)){
    w[k] <- w[k-1]/u[k-1]*(1-u[k-1])*u[k]
  }
  w
}

#' Show function dependencies
#' 
#' @param fname Character string for one function
#' @param pckg Package name; default is `"package:baker"`
#' @param ... Other parameters accepted by [mvbutils::foodweb()]
#' @return A figure showing function dependencies
#' @importFrom mvbutils foodweb
#' 
#' @examples
#' \dontrun{
#' show_dep("nplcm",ancestor=FALSE)
#' show_dep("nplcm",ancestor=FALSE)
#' show_dep("nplcm_fit_NoReg",ancestor=FALSE)
#' show_dep("nplcm_fit_NoReg")
#' }
#' 
#' @export

show_dep <- function(fname,pckg="package:baker",...){
  suppressWarnings(mvbutils::foodweb(where = pckg, prune = fname,
                                     #border = TRUE,
                                     #expand.xbox = 2, 
                                     boxcolor = "#FC6512",
                                     textcolor = "black", cex = 1.0, lwd=2,...))
  graphics::mtext(paste0("The ",fname," function foodweb"))
  
  # do not uncomment - just to say there is another cool way to visualize:
  # library(DependenciesGraphs) # if not installed, try this-- devtools::install_github("datastorm-open/DependenciesGraphs")
  # library(QualtricsTools) # devtools::install_github("emmamorgan-tufts/QualtricsTools")
  # dep <- funDependencies('package:baker','nplcm')
  # plot(dep)
}


#' check existence and create folder if non-existent
#' 
#' @param path Folder path to check and create if not there.
#' 
#' @export
check_dir_create <- function(path){
  if (file.exists(path)){
    return(NULL)
  } 
  dir.create(path)
}


nevercalled <- function(){
  ignored <- shinyFiles::getVolumes()
  ignored2 <- shinydashboard::box()
}


total_loc <- function(DIR="R"){
  files <- list.files(DIR,pattern="*.R", 
                      recursive=TRUE, full.names=TRUE)
  N <- 0
  for (f in files){N <- N+ length(readLines(f))}
  N
}

#' test if a formula has other terms not created by s_date_Eti or s_date_FPR
#' 
#' @param form a formula
#' @examples 
#' form1 <- as.formula(~ -1+s_date_FPR(DATE,Y,basis = "ps",10) + as.factor(SITE))
#' form2 <- as.formula(~ -1+s_date_FPR(DATE,Y,basis = "ps",10))
#' form3 <- as.formula(~ s_date_FPR(DATE,Y,basis = "ps",10))
#' 
#' has_non_basis(form1)
#' has_non_basis(form2)
#' has_non_basis(form3)
#' @export
has_non_basis <- function(form){
  out <- stats::terms(form)
  outlab <- attr(out,"term.labels")
  (attr(out,"intercept")>0) ||
    (length(grep("^s_",outlab))>=1 && length(outlab[-grep("^s_",outlab)])>=1 && attr(out,"intercept")==0) ||
     (length(outlab)>=1 && length(grep("^s_",outlab))==0 && attr(out,"intercept")==0) 
}

#' Make Etiology design matrix for dates with R format.
#'
#' `s_date_Eti` creates design matrices for etiology regressions; 
#' 
#' @param Rdate a vector of dates of R format
#' @param Y Binary case/control status; 1 for case; 0 for controls
#' @param basis "ncs" for natural cubic splines; "ps" for penalized-splines based
#' on B-spline basis functions (NB: baker does not recommend setting ncs using 
#' this function; use splines::ns)
#' @param dof Degree-of-freedom for the bases. For "ncs" basis, `dof` is
#' the number of columns; For "ps" basis,  the number of columns is `dof`
#' if `intercept=TRUE`; `dof-1` if `FALSE`.
#' @param ... Other arguments as in [splines::bs()]
#' @seealso [nplcm()]
#' @return 
#' \itemize{
#' \item `Z_Eti` design matrix for etiology regression on dates.
#' }
#' @importFrom stats quantile
#' @export
s_date_Eti <- function(Rdate,Y,basis = "ps",dof=ifelse(basis=="ncs",5,10),...) {
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
  grp_sd <- c(stats::sd(df$num_date[Y == 0]),stats::sd(df$num_date[Y == 1]))
  df$ingrp_std_num_date <-
    (df$num_date - grp_mean[df$Y + 1]) / grp_sd[df$Y + 1]
  #outgrp_std_num_date standardizes the cases' dates using controls' mean and stats::sd:
  df$outgrp_std_num_date <- (df$num_date - grp_mean[1]) / grp_sd[1]
  df$outgrp_std_num_date[df$Y == 0] <- NA
  
  case_ingrp_std_num_date <- df$ingrp_std_num_date[df$Y == 1]
  if (basis == "ncs") {
    # for etiology regression:
    Z_Eti <- splines::ns(case_ingrp_std_num_date,df = dof,Boundary.knots=
                           quantile(case_ingrp_std_num_date,c(0.05,0.95),...))
  }
  
  if (basis == "ps"){ # for penalized-splines based on B-spline basis:
    # for etiology regression:
    x       <- case_ingrp_std_num_date
    myknots <- stats::quantile(x,seq(0,1,length = (dof - 2))[-c(1,(dof - 2))])
    Z_Eti      <- matrix(splines::bs(x,knots= myknots,...),nrow=length(x))
    # if intercept=FALSE, then ncol(Z_Eti) = dof-1.
  }
  ind <- which(Y == 1)
  res <- matrix(NA,nrow = length(Y),ncol = ncol(Z_Eti))
  res[ind,]  <- Z_Eti
  # colnames(res) <- paste("date.basis.eti",1:ncol(Z_Eti),sep=".")
  res
}

#' Make false positive rate (FPR) design matrix for dates with R format.
#'
#' `s_date_FPR` creates design matrices for FPR regressions; 
#' 
#' @param Rdate a vector of dates of R format
#' @param Y Binary case/control status; 1 for case; 0 for controls
#' @param basis "ps" for penalized-splines based
#' on B-spline basis functions
#' @param dof Degree-of-freedom for the bases.For "ps" basis,  
#' the number of columns is `dof`
#' if `intercept=TRUE`; `dof-1` if `FALSE`.
#' @param ... Other arguments as in [splines::bs()]
#'
#' @importFrom stats quantile
#' @seealso [nplcm()]
#' @return Design matrix for FPR regression, with cases' rows on top of
#' controls'.
#' @export
s_date_FPR <- function(Rdate,Y,basis="ps",dof=10,...) {
  # standardization:
  df       <- data.frame(Y = Y,num_date = as.numeric(Rdate))
  grp_mean <- c(mean(df$num_date[df$Y == 0]),mean(df$num_date[df$Y == 1]))
  grp_sd   <- c(stats::sd(df$num_date[df$Y == 0]),stats::sd(df$num_date[df$Y == 1]))
  df$ingrp_std_num_date <-
    (df$num_date - grp_mean[df$Y + 1]) / grp_sd[df$Y + 1]
  #outgrp_std_num_date standardizes the cases' dates using controls' mean and stats::sd:
  df$outgrp_std_num_date <- (df$num_date - grp_mean[1]) / grp_sd[1]
  df$outgrp_std_num_date[df$Y == 0] <- NA
  
  ctrl_ingrp_std_num_date <- df$ingrp_std_num_date[df$Y == 0]
  # borrowing FPR regression from controls to cases:
  case_outgrp_std_num_date <- df$outgrp_std_num_date[df$Y == 1]
  if (basis=="ps"){
    x       <- ctrl_ingrp_std_num_date
    myknots <- stats::quantile(unique(x),seq(0,1,length = (dof - 2))[-c(1,(dof - 2))])
    Z_FPR_ctrl      <- matrix(splines::bs(x,knots= myknots),nrow=length(x))
    
    # borrowing FPR regression from controls to cases:
    x       <- case_outgrp_std_num_date
    Z_FPR_case      <- matrix(splines::bs(x,knots= myknots),nrow=length(x))
    
    ind <- which(df$Y == 1)
    res <- matrix(NA,nrow = length(df$Y),ncol = ncol(Z_FPR_case))
    res[ind,]  <- Z_FPR_case
    res[-ind,] <- Z_FPR_ctrl
  }
  # } else if (basis=="ncs"){
  #   x       <- ctrl_ingrp_std_num_date
  #   ncs_for_control <- splines::ns(x,df = dof,
  #               Boundary.knots = quantile(x,c(0.05,0.95)),...)
  #   Z_FPR_ctrl      <- matrix(ncs_for_control,nrow=length(x))
  #   
  #   # borrowing FPR regression from controls to cases:
  #   x       <- case_outgrp_std_num_date
  #   Z_FPR_case      <- matrix(splines::ns(x,knots=attr(ncs_for_control,"knots"),
  #           Boundary.knots = quantile(x,c(0.05,0.95)),...),nrow=length(x))
  #   
  #   ind <- which(df$Y == 1)
  #   res <- matrix(NA,nrow = length(df$Y),ncol = ncol(Z_FPR_case))
  #   res[ind,]  <- Z_FPR_case
  #   res[-ind,] <- Z_FPR_ctrl
  # }
  return(res)
}

# dudue_s_date_FPR <- function(Rdate,Y,basis="ncs",dof=10,...) {
#   # standardization:
#   df       <- data.frame(Y = Y,num_date = as.numeric(Rdate))
#   grp_mean <- c(mean(df$num_date[Y == 0]),mean(df$num_date[Y == 1]))
#   grp_sd   <- c(stats::sd(df$num_date[Y == 0]),stats::sd(df$num_date[Y == 1]))
#   df$ingrp_std_num_date <-
#     (df$num_date - grp_mean[df$Y + 1]) / grp_sd[df$Y + 1]
#   #outgrp_std_num_date standardizes the cases' dates using controls' mean and stats::sd:
#   df$outgrp_std_num_date <- (df$num_date - grp_mean[1]) / grp_sd[1]
#   df$outgrp_std_num_date[df$Y == 0] <- NA
#   
#   ctrl_ingrp_std_num_date <- df$ingrp_std_num_date[df$Y == 0]
#   # borrowing FPR regression from controls to cases:
#   case_outgrp_std_num_date <- df$outgrp_std_num_date[df$Y == 1]
#   if (basis=="ncs"){
#     x       <- ctrl_ingrp_std_num_date
#     Z_FPR_ctrl <- splines::ns(x,df = dof,Boundary.knots=
#                            quantile(x,c(0.05,0.95),...))
#     
#     # borrowing FPR regression from controls to cases:
#     x       <- case_outgrp_std_num_date
#     Z_FPR_case <- splines::ns(x,df = dof,Boundary.knots=
#                                 quantile(x,c(0.05,0.95),...))
#     
#     ind <- which(Y == 1)
#     res <- matrix(NA,nrow = length(Y),ncol = ncol(Z_FPR_case))
#     res[ind,]  <- Z_FPR_case
#     res[-ind,] <- Z_FPR_ctrl
#     return(res)
#   }
# }

#' Run ‘JAGS’ from R
#' 
#' The jags function takes data and starting values as input. 
#' It automatically writes a jags script, calls the model, and saves the 
#' simulations for easy access in R. Check the R2jags::jags2 for details about
#' the argument. 
#' 
#' This modifies the jags2 function in R2jags package. 
#' 
#' @inheritParams R2jags::jags
#' @import R2jags
#' @seealso [R2jags::jags()]
#' @export
jags2_baker <- function (data, inits, parameters.to.save, model.file = "model.bug", 
                         n.chains = 3, n.iter = 2000, n.burnin = floor(n.iter/2), 
                         n.thin = max(1, floor((n.iter - n.burnin)/1000)), DIC = TRUE, 
                         jags.path = "", working.directory = NULL, clearWD = TRUE, 
                         refresh = n.iter/50) 
{
  if (!is.null(working.directory)) {
    working.directory <- path.expand(working.directory)
    savedWD <- getwd()
    setwd(working.directory)
    on.exit(setwd(savedWD))
  }
  else {
    savedWD <- getwd()
    working.directory <- savedWD
  }
  redo <- ceiling(n.iter - n.burnin)
  if (is.list(data)) {
    data.list <- data
    lapply(data.list, dump, append = TRUE, file = "jagsdata.txt", 
           envir = parent.frame(1))
  }
  else {
    if (!(length(data) == 1 && is.vector(data) && is.character(data) && 
          (regexpr("\\.txt$", data) > 0))) {
      data.list <- lapply(as.list(data), get, pos = parent.frame(1))
      names(data.list) <- as.list(data)
    }
    else {
      if (all(basename(data) == data)) {
        try(file.copy(file.path(working.directory, data), 
                      data, overwrite = TRUE))
      }
      if (!file.exists(data)) {
        stop("File", data, "does not exist.")
      }
      data.list <- data
    }
  }
  lapply(names(data.list), dump, append = TRUE, file = "jagsdata.txt")
  data <- read.jagsdata("jagsdata.txt")
  if (is.function(model.file)) {
    temp <- tempfile("model")
    temp <- if (is.R() || .Platform$OS.type != "windows") {
      paste(temp, "txt", sep = ".")
    }
    else {
      gsub("\\.tmp$", ".txt", temp)
    }
    write.model(model.file, con = temp)
    model.file <- gsub("\\\\", "/", temp)
    if (!is.R()) 
      on.exit(file.remove(model.file), add = TRUE)
  }
  jags.call <- if (jags.path == "") {
    "jags"
  }
  else {
    jags.path <- win2unixdir(jags.path)
    paste(jags.path, "jags", sep = "")
  }
  no.inits <- FALSE
  inits.files <- NULL
  chain.names <- NULL
  if (is.null(inits)) {
    no.inits <- TRUE
  }
  else if (is.function(inits)) {
    for (i in 1:n.chains) {
      initial.values <- inits()
      inits.files <- c(inits.files, paste("jagsinits", 
                                          i, ".txt", sep = ""))
      chain.names <- c(chain.names, paste("chain(", i, 
                                          ")", sep = ""))
      curr_init_txt_file <- paste("jagsinits", i, ".txt", sep = "")
      with(initial.values, dump(names(initial.values), 
                                file = curr_init_txt_file))
      
      # fix dimension problem.... convert say 7:6 to c(7,6) (an issue for a dumped matrix):
      inits_fnames <- list.files(pattern = "^jagsinits[0-9]+.txt",
                                 full.names = TRUE)
      for (fiter in seq_along(inits_fnames)){
        curr_inits_txt_file <- inits_fnames[fiter]
        bad_jagsinits_txt <- readLines(curr_inits_txt_file)
        good_jagsinits_txt <- gsub( "([0-9]+):([0-9]+)", "c(\\1L,\\2L)", bad_jagsinits_txt,fixed = FALSE)
        writeLines(good_jagsinits_txt, curr_inits_txt_file)
      }
      
      
    }
  }
  else if (is.list(inits)) {
    if (length(inits) == n.chains) {
      for (i in 1:n.chains) {
        initial.values <- inits[[i]]
        inits.files <- c(inits.files, paste("jagsinits", 
                                            i, ".txt", sep = ""))
        chain.name <- c(chain.names, paste("chain(", 
                                           i, ")", sep = ""))
        with(initial.values, dump(names(initial.values), 
                                  file = paste("jagsinits", i, ".txt", sep = "")))
      }
    }
    else {
      stop(message = "initial value must be specified for all of chains")
    }
  }
  if (DIC) {
    parameters.to.save <- c(parameters.to.save, "deviance")
  }
  cat("model clear\ndata clear\n", if (DIC) {
    "load dic\n"
  }, "model in ", "\"", model.file, "\"", "\n", "cd ", "\"", 
  working.directory, "\"", "\n", "data in ", "\"jagsdata.txt\"", 
  "\n", "compile, nchains(", n.chains, ")", "\n", if (!no.inits) {
    paste("inits in \"", inits.files, "\", ", chain.names, 
          "\n", sep = "")
  }, "initialize", "\n", "update ", n.burnin, ", by(", 
  refresh, ")\n", paste("monitor ", parameters.to.save, 
                        ", thin(", n.thin, ")\n", sep = ""), "update ", redo, 
  ", by(", refresh, ")\n", "coda *\n", sep = "", file = "jagsscript.txt")
  system(paste(jags.call, "jagsscript.txt"))
  # fit <- R2jags:::jags.sims(parameters.to.save = parameters.to.save, 
  #                  n.chains = n.chains, n.iter = n.iter, n.burnin = n.burnin, 
  #                  n.thin = n.thin, DIC = DIC)
  fit <- jags.sims(parameters.to.save = parameters.to.save, 
                   n.chains = n.chains, n.iter = n.iter, n.burnin = n.burnin, 
                   n.thin = n.thin, DIC = DIC)  
  if (clearWD) {
    file.remove(c("jagsdata.txt", "CODAindex.txt", inits.files, 
                  "jagsscript.txt", paste("CODAchain", 1:n.chains, 
                                          ".txt", sep = "")))
  }
  class(fit) <- "bugs"
  return(fit)
}

#' log sum exp trick
#' 
#' @param x a vector of numbers
#' @export
logsumexp <- function (x) {
  y = max(x)
  y + log(sum(exp(x - y)))
}


#' softmax
#' 
#' uses logsumexp trick to prevent numerical overflow
#' 
#' @param x a vector of numbers
#' @examples
#' 
#' softmax2 <- function(x) exp(x) / sum(exp(x))
#' softmax(c(1, 2, 3) * 1000)  # NaN NaN NaN
#' softmax2(c(1, 2, 3) * 1000)  # 0 0 1
#' 
#' @export
softmax <- function (x) {
  exp(x - logsumexp(x))
}


#' subset data from the output of [clean_perch_data()]
#'
#' It is particularly useful in simulating data from a regression model where one
#' generates a case and control at a particular covariate value, and just choose
#' a case or control to retain in the simulated data.
#'
#' @param data_nplcm data for fitting nplcm; See [nplcm]
#' @param index a vector of indices indicating the observations you hope to subset;
#' it will subset in all the sublists of data_nplcm
#'
#' @return a list with the requested data, in the order determined by 'index'
#'
#' @examples 
#' 
#'J = 3                          # number of causes
#'cause_list = c(LETTERS[1:J])   # cause list
#'K = 2                          # number of subclasses
#'lambda = c(1,0)                # subclass weights for control group
#'eta = c(1,0)                   # subclass weights for case group
#'
#'# setup parameters for the present individual:
#'set_parameter <- list(
#'  cause_list      = cause_list,
#'  etiology        = c(0.5,0.2,0.3), # only meaningful for cases 
#'  pathogen_BrS    = LETTERS[1:J],
#'  pathogen_SS     = LETTERS[1:2],
#'  meas_nm         = list(MBS = c("MBS1"),MSS=c("MSS1")),
#'  Lambda          = lambda,         # for BrS   
#'  Eta             = t(replicate(J,eta)),  # case subclass weight for BrS
#'  PsiBS           = cbind(c(0.15,0.3,0.35),   
#'                          c(0.25,0.2,0.15)), # FPR
#'  PsiSS           = cbind(rep(0,J),rep(0,J)),
#'  ThetaBS         = cbind(c(0.95,0.9,0.85),    # TPR
#'                          c(0.95,0.9,0.85)),
#'  ThetaSS         = cbind(c(0.25,0.10),
#'                          c(0.25,0.10)),
#'  Nd      =     5,
#'  Nu      =     3 
#')
#'simu_out   <- simulate_nplcm(set_parameter)
#'out <- simu_out$data_nplcm
#'out
#'subset_data_nplcm_by_index(out,c(1,4,5))
#'subset_data_nplcm_by_index(out,2)
#'
#' @family data operation functions
#' @export
subset_data_nplcm_by_index <- function(data_nplcm,index) {
  n_data_quality <- length(data_nplcm$Mobs)
  Mobs  <- list()
  for (i in seq_along(data_nplcm$Mobs)) {
    if (is.null(data_nplcm$Mobs[[i]])){
      Mobs[[i]] <- NULL
    } else{
      Mobs[[i]] <- lapply(data_nplcm$Mobs[[i]],function(df)
        df[index,])
    }
  }
  names(Mobs) <- c("MBS","MSS","MGS")[1:length(Mobs)]
  X_df <- data_nplcm$X[index,,drop=FALSE]
  names(X_df) <- names(data_nplcm$X)
  list(Mobs = Mobs,
       Y    = data_nplcm$Y[index],
       X    = X_df)
}


#' For a list of many sublists each of which has matrices as its member, 
#' we combine across the many sublists to produce a final list
#' 
#' @param list_of_lists a list of sublists
#' 
#' @examples 
#' DT1 = list(A=1:3,B=letters[1:3])
#' DT2 = list(A=4:5,B=letters[4:5])
#' DT3 = list(A=1:4,B=letters[1:4])
#' DT4 = list(A=4:7,B=letters[4:7])
#' l = list(DT1,DT2);names(l) <- c("haha","hihi")
#' l2 = list(DT3,DT4);names(l2) <- c("haha","hihi")
#' listoflists <- list(l,l2);names(listoflists) <- c("dude1","dude2")
#' listoflists
#' merge_lists(listoflists)
#' @family data operation functions
#' @return a list after merge
#' @export
merge_lists <- function(list_of_lists){
  if (length(unique(lapply(list_of_lists,length)))>1){
    stop("==[baker] the list of lists have distinct lengths.==")
  }
  len_one_list <- length(list_of_lists[[1]])
  res <- lapply(1:len_one_list,function(l) {
    do.call(rbind,lapply(1:length(list_of_lists),
                         function(s) list_of_lists[[s]][[l]]))})
  names(res) <- names(list_of_lists[[1]])
  res
}

#' combine multiple data_nplcm (useful when simulating data from regression models)
#' 
#' @param data_nplcm_list a list of data_nplcm in [nplcm]
#' 
#' @examples 
#' 
#' \dontrun{
#' rm(list=ls())
#' N=10000
#' Y = rep(c(1,0),times=5000) # simulate two cases and two controls.
#' out_list <- vector("list",length=N)
#' J = 3                          # number of causes
#' cause_list = c(LETTERS[1:J])   # cause list
#' K = 2                          # number of subclasses
#' lambda = c(.8,.2)                # subclass weights for control group
#' eta = c(.9,.1)                   # subclass weights for case group
#' 
#' for (i in 1:N){
#'   #setup parameters for the present individual:
#'   set_parameter <- list(
#'     cause_list      = cause_list,
#'     etiology        = c(0.5,0.2,0.3), # only meaningful for cases 
#'     pathogen_BrS    = LETTERS[1:J],
#'     pathogen_SS     = LETTERS[1:2],
#'     meas_nm         = list(MBS = c("MBS1"),MSS=c("MSS1")),
#'     Lambda          = lambda,         # for BrS   
#'     Eta             = t(replicate(J,eta)),  # case subclass weight for BrS
#'     PsiBS           = cbind(c(0.15,0.3,0.35),   
#'                             c(0.25,0.2,0.15)), # FPR
#'     PsiSS           = cbind(rep(0,J),rep(0,J)),
#'     ThetaBS         = cbind(c(0.95,0.9,0.85),    # TPR
#'                             c(0.95,0.9,0.85)),
#'     ThetaSS         = cbind(c(0.25,0.10),
#'                             c(0.25,0.10)),
#'     Nd      =     1,
#'     Nu      =     1 
#'   )
#'   simu_out   <- simulate_nplcm(set_parameter)
#'   out <- simu_out$data_nplcm
#'   out_list[[i]] <- out
#' }
#' 
#' # extract cases and controls and combine all the data into one:
#' data_nplcm_list <- lapply(1:N, function(s) subset_data_nplcm_by_index(out_list[[s]],2-Y[s]))
#' data_nplcm_unordered      <- combine_data_nplcm(data_nplcm_list) 
#' data_nplcm_unordered
#' 
#' colMeans(data_nplcm_unordered$Mobs$MBS$MBS1[Y==1,])
#' colMeans(data_nplcm_unordered$Mobs$MBS$MBS1[Y==0,])
#' 
#' set_parameter$PsiBS%*%matrix(lambda,ncol=1)
#' set_parameter$PsiBS%*%matrix(eta,ncol=1)*(1-set_parameter$etiology)+
#'   set_parameter$ThetaBS%*%matrix(eta,ncol=1)*set_parameter$etiology
#' 
#' # data_nplcm <- subset_data_nplcm_by_index(data_nplcm_unordered,
#' #                    order(-data_nplcm_unordered$Y)) #put cases on top.
#' 
#' }
#' 
#' @family data operation functions
#' @export
combine_data_nplcm <- function(data_nplcm_list){
  if (length(data_nplcm_list)==1) {stop("==[baker]==list is of length 1, no need to combine.")}
  n_data_quality <- length(data_nplcm_list[[1]]$Mobs)
  Mobs <- list()
  for (i in 1:n_data_quality) {
    Mobs[[i]] <- merge_lists(lapply(lapply(data_nplcm_list,
                                           function(l) l$Mobs),function(s) s[[i]]))
  }
  names(Mobs) <- c("MBS","MSS","MGS")[1:n_data_quality]
  X_df <- do.call(rbind,lapply(data_nplcm_list,function(l) l$X))
  names(X_df) <- names(data_nplcm_list[[1]]$X)
  list(Mobs = Mobs,
       Y    = c(unlist(lapply(data_nplcm_list,function(l) l$Y))),
       X    = X_df)
}



#' convert line to user coordinates
#' 
#' Here's a version that works with log-scale and linear scale axes. 
#' The trick is to express line locations in npc coordinates rather than user coordinates, 
#' since the latter are of course not linear when axes are on log scales.
#' 
#' `par('cin')[2] * par('cex') * par('lheight')` returns the current line height
#'  in inches, which we convert to user coordinates by multiplying by 
#'  `diff(grconvertX(0:1, 'inches', 'user'))`, the length of an inch in user
#'   coordinates (horizontally, in this case - if interested in the vertical
#'    height of a line in user coords we would use 
#'    `diff(grconvertY(0:1, 'inches', 'user')))`.
#'    
#'    
#' @param line integer
#' @param side integer; 1-4
#' @references <https://stackoverflow.com/questions/29125019/get-margin-line-locations-mgp-in-user-coordinates>
#' @export
#' @examples 
#' 
#' 
#' setup_plot <- function(log = "") {
#'   par(mar = c(2, 10, 2, 2), oma = rep(2, 4))
#'   plot.new()
#'   plot.window(xlim = c(1, 10), ylim = c(1, 10), log = log)
#'   box(which = "plot", lwd = 2, col = "gray40")
#'   box(which = "figure", lwd = 2, col = "darkred")
#'   box(which = "outer", lwd = 2, col = "darkgreen")
#'   text(x = 0.5, y = 0.5, 
#'        labels = "Plot Region", 
#'        col = "gray40", font = 2)
#'   mtext(side = 3, text = "Figure region", line = 0.5, col = "darkred", font = 2)
#'   mtext(side = 3, text = "Device region", line = 2.5, col = "darkgreen", font = 2)
#'   for (i in 0:9) {
#'     mtext(side = 2, col = "darkred", text = paste0("Line", i), line = i)
#'   }
#' }
#' # And here are a couple of examples, applied to your setup_plot with mar=c(5, 5, 5, 5):
#' setup_plot()
#' axis(1, line=5)
#' axis(2, line=5)
#' abline(h=line2user(0:4, 1), lty=3, xpd=TRUE)
#' abline(v=line2user(0:4, 2), lty=3, xpd=TRUE)
#' abline(h=line2user(0:4, 3), lty=3, xpd=TRUE)
#' abline(v=line2user(0:4, 4), lty=3, xpd=TRUE)
#' 
#' setup_plot(log='x')
#' axis(1, line=5)
#' axis(2, line=5)
#' abline(h=line2user(0:4, 1), lty=3, xpd=TRUE)
#' abline(v=line2user(0:4, 2), lty=3, xpd=TRUE)
#' abline(h=line2user(0:4, 3), lty=3, xpd=TRUE)
#' abline(v=line2user(0:4, 4), lty=3, xpd=TRUE)
#' 
#' 
#' setup_plot(log='y')
#' axis(1, line=5)
#' axis(2, line=5)
#' abline(h=line2user(0:4, 1), lty=3, xpd=TRUE)
#' abline(v=line2user(0:4, 2), lty=3, xpd=TRUE)
#' abline(h=line2user(0:4, 3), lty=3, xpd=TRUE)
#' abline(v=line2user(0:4, 4), lty=3, xpd=TRUE)
#' 
#' setup_plot(log='xy')
#' axis(1, line=5)
#' axis(2, line=5)
#' abline(h=line2user(0:4, 1), lty=3, xpd=TRUE)
#' abline(v=line2user(0:4, 2), lty=3, xpd=TRUE)
#' abline(h=line2user(0:4, 3), lty=3, xpd=TRUE)
#' abline(v=line2user(0:4, 4), lty=3, xpd=TRUE)
line2user <- function(line, side) {
  lh <- par('cin')[2] * par('cex') * par('lheight')
  x_off <- diff(grconvertX(c(0, lh), 'inches', 'npc'))
  y_off <- diff(grconvertY(c(0, lh), 'inches', 'npc'))
  switch(side,
         `1` = grconvertY(-line * y_off, 'npc', 'user'),
         `2` = grconvertX(-line * x_off, 'npc', 'user'),
         `3` = grconvertY(1 + line * y_off, 'npc', 'user'),
         `4` = grconvertX(1 + line * x_off, 'npc', 'user'),
         stop("Side must be 1, 2, 3, or 4", call.=FALSE))
}


##################################################################
# get-metric.R
##################################################################

#' Obtain Integrated Squared Aitchison Distance, Squared Bias and Variance (both on 
#' Central Log-Ratio transformed scale) that measure the discrepancy of a posterior
#' distribution of pie and a true pie. 
#'
#' The result is equivalent to Euclidean-type calculation after the
#' compositional vector (e.g., etiologic fraction) is centered-log-ratio (CLRB) transformed.
#' For simulation only.
#'
#' @param DIR_NPLCM File path where Bayesian results are stored
#' @param truth True etiologic fraction vector (must sum to 1) 
#' used to generate data
#' @importFrom robCompositions cenLR
#' @return a vector of (Integrated Squared Aitchison Distance (ISAD), bias-squared, variance, truth)
#' @export
#'
get_metric <- function(DIR_NPLCM,truth){
  use_jags <- is_jags_folder(DIR_NPLCM)
  my.mse <- function(A,B){
    if (nrow(A)!=nrow(B)){
      stop("not equal number of rows!")
    }else{
      aa=cenLR(A)$x.clr
      bb=cenLR(B)$x.clr
      res <- rowSums((aa-bb)^2)
      mean(res)
    }
  }
  
  my.bias.sq <- function(A,B){
    if (nrow(A)!=nrow(B)){
      stop("not equal number of rows!")
    }else{
      aa=cenLR(A)$x.clr
      bb=cenLR(B)$x.clr
      sum((colMeans(aa-bb))^2)
    }
  }
  
  #A = as.data.frame(true_pEti)
  my.var <- function(A){
    res <- cenLR(A)$x.clr
    sum(diag(stats::cov(res)))
  }
  
  #read in data from the folder directory:
  
  out <- nplcm_read_folder(DIR_NPLCM)
  model_options <- out$model_options
  bugs.dat      <- out$bugs.dat
  res_nplcm     <- out$res_nplcm
  #some data preparation:
  Nd <- bugs.dat$Nd
  Nu <- bugs.dat$Nu
  
  Y <-  out$Y
  
  cause_list     <- model_options$likelihood$cause_list
  Jcause         <- length(grep("^pEti\\[",colnames(res_nplcm)))
  if (Jcause != length(cause_list)){
    stop("=='model_options' and actual posterior results 'pEti' have different number of causes!==")
  }
  
  # extract and process some data and posterior samples:
  SubVarName <- rep(NA,Jcause)
  for (j in 1:Jcause){
    SubVarName[j] = paste("pEti","[",j,"]",sep="")
  }
  #get etiology fraction MCMC samples:
  pEti_mat   <- as.matrix(res_nplcm[,SubVarName])
  true_pEti  <- as.data.frame(t(replicate(nrow(pEti_mat),truth)))
  
  # Aitchison distances:
  # aDist(true_pEti,pEti_mat)
  mse     <- my.mse(true_pEti,pEti_mat)
  bias.sq <- my.bias.sq(true_pEti,pEti_mat)
  vv      <- my.var(as.data.frame(pEti_mat))
  
  res <- make_list(mse, bias.sq,vv,truth)
  res
}

#' Obtain direct bias that measure the discrepancy of a posterior
#' distribution of pie and a true pie. 
#'
#' @param DIR_list The list of  where Bayesian results are stored
#' @param truth True etiologic fraction vector (must sum to 1)  used to generate data;
#' Default is `NULL`. If a vector is supplied, then only the first path in `DIR_LIST`
#' is used.
#' @param silent Default is FALSE. To suppress printing messages, set to TRUE.
#' 
#' @return a list of length two. `diff` is the direct differences; 
#' `prb` is the percent relative bias.
#' @export
#'
get_direct_bias <- function(DIR_list,truth=NULL,silent=FALSE){
  
  symdiff <- function( x, y) { setdiff( union(x, y), intersect(x, y))}
  
  if (is.null(truth)){
    # read from folders:
    out_list <- vector("list",length(DIR_list))
    base_nm  <- lapply(DIR_list,basename)
    names(out_list) <- base_nm
    
    pEti_samp_list  <- list()
    for (i in seq_along(DIR_list)){
      out_list[[i]]        <- nplcm_read_folder(DIR_list[[i]])
      pEti_samp_list[[i]]  <- get_pEti_samp(out_list[[i]]$res_nplcm,out_list[[i]]$model_options)
    }
    
    different_names <- symdiff(out_list[[1]]$model_options$likelihood$cause_list,out_list[[2]]$model_options$likelihood$cause_list)
    if (length(different_names) > 0 ){stop("==Two folders have different latent category names!==")}
    
    if(!silent){
      print("==The first folder in 'DIR_LIST' is used as reference for calculating relative bias. ==")
    }
    # plain difference:
    diff_mat <- pEti_samp_list[[2]]$pEti_mat - pEti_samp_list[[1]]$pEti_mat
    diff     <- colMeans(diff_mat)
    names(diff) <- pEti_samp_list[[1]]$latent_nm
    # percent relative bias:
    prb     <- colMeans(diff_mat)/colMeans(pEti_samp_list[[1]]$pEti_mat)
    names(prb) <- pEti_samp_list[[1]]$latent_nm
    
    res <- make_list(diff,prb)
    return(res)
  }
  
  if (!is.null(truth)){
    if(!silent){print("==Only the first folder in 'DIR_LIST' is used to compare with 'truth'!==")}
    # read from folders:
    out_list <- vector("list",length(DIR_list))
    base_nm  <- lapply(DIR_list,basename)
    names(out_list) <- base_nm
    
    pEti_samp_list  <- list()
    for (i in 1){
      out_list[[i]]        <- nplcm_read_folder(DIR_list[[i]])
      pEti_samp_list[[i]]  <- get_pEti_samp(out_list[[i]]$res_nplcm,out_list[[i]]$model_options)
    }
    
    if (length(out_list[[1]]$model_options$likelihood$cause_list) != length(truth)){
      stop("==Results and truth have different number of latent categories! ==")
    }
    
    # plain difference:
    diff_mat <- pEti_samp_list[[1]]$pEti_mat - t(replicate(nrow(pEti_samp_list[[1]]$pEti_mat),truth))
    diff     <- colMeans(diff_mat)
    names(diff) <- pEti_samp_list[[1]]$latent_nm
    # percent relative bias:
    prb     <- colMeans(diff_mat)/truth
    names(prb) <- pEti_samp_list[[1]]$latent_nm
    
    res <- make_list(diff,prb)
    return(res)
  }
}

#' Obtain coverage status from a result folder
#'
#' @param DIR_NPLCM Path to where Bayesian results are stored
#' @param truth True etiologic fraction vector (must sum to 1)  used to generate data.
#' 
#' @examples 
#' \dontrun{
#' DIR_NPLCM <- "~/downloads/rep_1_kfit_2/"  
#' truth     <- c(0.5,0.2,0.15,0.1,0.05)
#' get_coverage(DIR_NPLCM,truth)
#' }
#' @return A logic vector of length as `truth`. 1 for covered; 0 for not.
#' @export
#'
get_coverage <- function(DIR_NPLCM,truth){
  # read from folders:
  out       <- nplcm_read_folder(DIR_NPLCM)
  pEti_samp <- get_pEti_samp(out$res_nplcm,out$model_options)
  
  if (length(out$model_options$likelihood$cause_list) != length(truth)){
    stop("==Results and truth have different number of latent categories! ==")
  }
  
  # plain difference:
  diff_mat <- pEti_samp$pEti_mat - t(replicate(nrow(pEti_samp$pEti_mat),truth))
  UL <- apply(diff_mat,2,stats::quantile,0.975)
  LL <- apply(diff_mat,2,stats::quantile,0.025)
  res <- (UL>=0 & LL <=0)
  return(res)
}


#' Obtain posterior standard deviation from a result folder
#'
#' @param DIR_NPLCM Path to where Bayesian results are stored
#' 
#' @examples 
#' \dontrun{
#' DIR_NPLCM <- "~/downloads/rep_1_kfit_2/"  
#' get_postsd(DIR_NPLCM)
#' }
#' @return a vector of positive numbers
#' @export
#'
get_postsd <- function(DIR_NPLCM){
  # read from folders:
  out       <- nplcm_read_folder(DIR_NPLCM)
  pEti_samp <- get_pEti_samp(out$res_nplcm,out$model_options)
  apply(pEti_samp$pEti_mat, 2, sd)
}



## moved from plot_etiology_side_by_side.R




#' get etiology samples by names (no regression)
#' 
#' @inheritParams order_post_eti
#' 
#' @return A list:
#' \itemize{
#'    `pEti_mat`: a matrix of posterior samples (iteration by cause); overall etiology
#'    `latent_nm`: a vector of character strings representing the names of the causes
#' }
#' 
#' @export

get_pEti_samp <- function(res_nplcm,model_options){
  cause_list <- model_options$likelihood$cause_list
  # total no. of causes:
  Jcause     <- length(cause_list)
  # extract and process some data and posterior samples:
  SubVarName <- rep(NA,Jcause)
  for (j in 1:Jcause){
    SubVarName[j] = paste("pEti","[",j,"]",sep="")
  }
  # get etiology fraction MCMC samples:
  pEti_mat   <- as.data.frame(res_nplcm[,SubVarName,drop=FALSE])
  latent_nm  <- model_options$likelihood$cause_list
  make_list(pEti_mat,latent_nm)
}



#' Match latent causes that might have the same combo but
#' different specifications
#' 
#'  @details In our cause_list, "A+B" represents the same cause
#'   as "B+A". It is used for plotting side-by-side posterior sample comparisons
#' 
#' @param pattern a vector of latent cause names, e.g., from a particular fit
#' @param vec a vector of latent cause names, e.g., usually a union of cause names
#' from several model fits. Usually, it is also the display order that one wants to 
#' show.
#' 
#' @return A vector of length `length(vec)`; `NA` means no pattern matches
#' vec; 1 at position 10 means the first element of `pattern` matches the 
#' 10th element of `vec`.
#' 
#' 
#' @examples 
#' 
#' pattern <- c("X+Y","A+Z","C")
#' vec     <- c(LETTERS[1:26],"Y+Z","Y+X","Z+A")
#' match_cause(pattern,vec)
#' 
#' @export
match_cause <- function(pattern, vec){
  has_plus_sign <- grepl("\\+",pattern)
  vec_split <- strsplit(vec,"\\+")
  res <- rep(NA,length(vec))
  for (p in seq_along(pattern)){
    tmp <- pattern[p]
    if (has_plus_sign[p]) {
      tmp <- strsplit(pattern[p],"\\+")[[1]]
    }
    find_match <- lapply(vec_split,setequal,y=tmp)
    res[which(find_match==TRUE)] <- p
  }
  res
}



#' get unique causes, regardless of the actual order in combo
#' 
#' @param cause_vec a vector of characters with potential combo repetitions
#' written in scrambled orders separated by "+"
#' 
#' @examples 
#' x <- c("A","B","A","CC+DD","DD+CC","E+F+G","B")
#' unique_cause(x)
#' 
#' @return a vector of characters with unique meanings for latent causes
#' 
#' @export

unique_cause <- function(cause_vec){
  cause_split <- strsplit(cause_vec,"\\+")
  num <- lapply(cause_split,length)
  res <- list()
  res[[1]] <- cause_split[[1]] 
  
  if (length(cause_vec)==1){return(cause_vec)}
  
  j <- 1
  for (i in 2:length(cause_vec)){
    got_new <- TRUE
    for (iter in 1:j){
      if (setequal(cause_split[[i]],res[[iter]])){got_new <- FALSE}
    }
    if (got_new){j <- j+1;res[[j]]<-cause_split[[i]]}
  }
  unlist(lapply(res,paste,collapse="+"))
}






