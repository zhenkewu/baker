#' Create new folder name
#'
#' @param parent_path The parent directory where to put the new folder
#' @param parameter_names The parameters that distinguish this folder's scenario
#' @param parameter_vals The actual parameter values
#'
#' @return A string for folder name
#'
#' @export
#'
make_foldername <-
function(parent_path,parameter_names,parameter_vals){
  subfolder <- paste(parameter_names,parameter_vals,collapse="_",sep="=")
  res       <- paste(parent_path,subfolder,sep="\\")
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
make_filename<-
function(parameter_names,parameter_vals,format){
  res1 <- paste(parameter_names,parameter_vals,collapse="_",sep="=")
  res  <- paste(res1,format,sep=".")
  res
}




#' logit function
#'
#' @param p Probability between 0 and 1
#' @return A real number
#'
#' @export
logit <- function(p) log(p)-log(1-p)

#' expit function
#'
#' @param x A real number
#' @return a Probability between 0 and 1
#' @export
expit <- function(x) 1/(1+exp(-x))





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
rvbern<-
function(p){
  U  <-runif(length(p),0,1)
  res<-(U<p)+0
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
bin2dec <-
function(binary_vector) {
  sum(2^(which(rev(binary_vector)==TRUE)-1))
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
function(pathogen_name,pathogen_list){
        J <- length(pathogen_list)
  splited <- strsplit(pathogen_name,split="+",fixed=TRUE)
  deploy  <- function(inst,J) {
    if (inst[1] == "NoA"){
      rep(0,J)
    }else{
      sapply(inst,grep,x=pathogen_list)
    }
  }
  res <- lapply(splited,deploy,J=J)
  for (l in seq_along(res)){
    any_integer_0 <- sum(sapply(res[[l]],function(v) length(v)==0))>0
    if (any_integer_0){
      stop(paste0("\n==",l,"-th cause, \n",pathogen_name[l],",
                  \n has pathogen(s) not in the overall pathogen list!=="))
    }
  }

  nc  <- length(res)
  matres <- t(sapply(1:nc,function(i) {
    tempres <- rep(0,J)
    tempres[res[[i]]]<-1
    tempres}))
  matres
}

#symb2I("A+D",c("A","B","C")) # will return error.

#' Convert 0/1 coding to pathogen/combinations
#'
#' Reverse to \code{\link[nplcm]{symb2I}}
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
I2symb <- function(binary_code,pathogen_list){
  ind <- grep("1",strsplit(binary_code,split="")[[1]])
  res <- ifelse(length(ind)==0,"NoA",paste(pathogen_list[ind],collapse="+"))
  res
}



#' Convert a matrix of binary indicators to categorial variables
#'
#' @param binary_mat The matrix of binary indicators. Rows for subjects, columns for pathogens in the \code{"pathogen.list"}
#' @param allowed_list The list of allowed single pathogen or pathogen combination
#' @param pathogen_list The complete list of pathogen names
#'
#' @return A vector of categorical variables. Its length equals the length of \code{"allowed.list"}
#'
#' @examples
#'
#' Imat2cat(rbind(diag(3),c(1,1,0),c(0,0,0)),c("A","B","C","A+B","NoA"),c("A","B","C"))
#' @export
Imat2cat <- function(binary_mat,allowed_list,pathogen_list){
  known_code = apply(binary_mat,1,function(v) paste(v,collapse=""))
  known_symb = sapply(known_code,I2symb,pathogen_list)
  if (sum(known_symb %in% allowed_list==FALSE)>0){
  stop("Some binary pattern in 'binary_mat' is not included by 'allowed_list'.")
  }else {
  known_Icat = sapply(known_symb,function(s) which(allowed_list==s))
  return(known_Icat)
  }
}

#' Plot beta density
#' @param a The first parameter
#' @param b The second parameter
#' @return None
#'
#' @export
beta_plot=function(a,b){
  x= seq(0,1,by=0.001)
  y = dbeta(x,a,b)
  plot(x,y,type="l",main=paste0("a=",a,",b=",b))
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
unfactor <- function(f){ as.numeric(levels(f))[f]}


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

sort_data_frame <- function(form,dat){
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
  if(inherits(dat,"formula")){
    f=dat
    dat=form
    form=f
  }
  if(form[[1]] != "~") {
    stop("Formula must be one-sided.")
  }
  # Make the formula into character and remove spaces
  formc <- as.character(form[2])
  formc <- gsub(" ","",formc)
  # If the first character is not + or -, add +
  if(!is.element(substring(formc,1,1),c("+","-"))) {
    formc <- paste("+",formc,sep="")
  }
  # Extract the variables from the formula
  vars <- unlist(strsplit(formc, "[\\+\\-]"))
  vars <- vars[vars!=""] # Remove spurious "" terms
  # Build a list of arguments to pass to "order" function
  calllist <- list()
  pos=1 # Position of + or -
  for(i in 1:length(vars)){
    varsign <- substring(formc,pos,pos)
    pos <- pos+1+nchar(vars[i])
    if(is.factor(dat[,vars[i]])){
      if(varsign=="-")
        calllist[[i]] <- -rank(dat[,vars[i]])
      else
        calllist[[i]] <- rank(dat[,vars[i]])
    }
    else {
      if(varsign=="-")
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
my_reorder <- function(disp_order,raw_nm){
        #disp_order <- display_order
        #raw_nm <- pathogen_BrS
        res <- rep(NA,length(raw_nm))
        for (i in seq_along(raw_nm)){
            res[i]<-which(disp_order==raw_nm[i])
        }
        order(res)
}

# disp_order <- c("B","E","D","C","F","A")
#   #union_nm   <- c("A","B","C","D","E")
#   raw_nm <- c("C","A","E")
#   my_reorder(disp_order,raw_nm)
# # list(c("C","E"),
# #      c("D","B"),



logOR <- function(MBS.case,MBS.ctrl,JBrS){  
  JBrS <- ncol(MBS.case)
  logORmat       <- matrix(NA,nrow=JBrS,ncol=JBrS)
  logORmat.se    <- matrix(NA,nrow=JBrS,ncol=JBrS)
  for (j2 in 1:(JBrS-1)){ #case (j2,j1); ctrl (j1,j2).
    for (j1 in (j2+1):JBrS){
      
      # cases: (upper triangle)
      x = MBS.case[,j2]
      y = MBS.case[,j1]
      
      fit = glm(y~x,family = binomial(link="logit"))
      
      if ("x" %in% rownames(summary(fit)$coef)){
        logORmat[j2,j1] = round(summary(fit)$coef["x",1],3)
        logORmat.se[j2,j1] = round(summary(fit)$coef["x",2],3)
      }
      # controls: (lower triangle)
      x = MBS.ctrl[,j2]
      y = MBS.ctrl[,j1]
      
      fit = glm(y~x,family = binomial(link="logit"))
      
      if ("x" %in% rownames(summary(fit)$coef)){
        logORmat[j1,j2] = round(summary(fit)$coef["x",1],3)
        logORmat.se[j1,j2] = round(summary(fit)$coef["x",2],3)
      }
    }
  }
  
  #cell.num = logORmat/logORmat.se
  tmp       = logORmat
  tmp[abs(logORmat.se)>10]=NA
  
  tmp.se = logORmat.se
  tmp.se[abs(logORmat.se)>10] = NA
  
  
  res <- list(tmp,tmp.se)
  names(res) <- c("logOR","logOR.se")
  return(res)
}

#' Visualize matrix for a quantity measured on cases and controls 
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
#' 
#' @export
visualize_case_control_matrix <- function(mat, dim_names=ncol(mat), 
                                   cell_metrics="",folding_line=TRUE,
                                   axes = FALSE, xlab = "",ylab = "",
                                   asp = 1,title="",...){
  n = nrow(mat)
  J = n
  # size of the numbers in the boxes:
  cex_main= min(2,20/n)
  cex_se  = min(1.5,15/n)
  
  par(mar = c(0, 0, 5, 0), bg = "white",xpd=TRUE)
  plot(c(0, n + 0.8), c(0, n + 0.8), axes = axes, xlab = "",
       ylab = "", asp = 1, type = "n")
  ##add grid
  segments(rep(0.5, n + 1), 0.5 + 0:n, rep(n + 0.5, n + 1),
           0.5 + 0:n, col = "gray")
  segments(0.5 + 0:n, rep(0.5, n + 1), 0.5 + 0:n, rep(n + 0.5,
                                                      n), col = "gray")  
  mat.txt<- round(t(mat)[,n:1],1)
  mat.txt3<- round(t(mat)[,n:1],3)
  
  # add meaning of the number in a cell:
  text(1,J,cell_metrics,cex=cex_main/2)
  
  for (i in 1:n){
    for (j in 1:n){
      abs.mat <- abs(mat.txt[i,j])
      if ((!is.na(abs.mat)) && abs.mat>2){
        text(i,j,round(mat.txt3[i,j],1),
             col=ifelse(mat.txt3[i,j]>0,"red","blue"),cex=cex_main)
      }
      
    }
  }
  
  if (folding_line){
  # diagonal line:
    segments(0.5+1,.5+n-1,.5+n,0.5,col="black",lty=3,lwd=3)
  }

  #mtext(title,3,cex=1,line=1)
  
  # put pathogen names on rows and columns:
  for (s in 1:J){
    text(-0,J-s+1,paste0(dim_names[s],":(",s,")"),cex=min(1.5,20/J),adj=1)
    text(s,J+0.7,paste0("(",s,"):",dim_names[s]),cex=min(1.5,20/J),srt=45,adj=0)
  }
  # labels for cases and controls:
  text(J+1,J/2,"cases",cex=2,srt=-90)
  text(J/2,0,"controls",cex=2)
}


NA2dot <- function(s){
  gsub("NA",".",s,fixed=TRUE)
}
