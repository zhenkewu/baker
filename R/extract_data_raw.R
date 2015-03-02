#' Extract data on pathogen measurements and covariates
#'
#'
#' @param Pathogen The vector of pathogen names. It will be reordered according
#'        according to pathogen_type (alphabetical order, "B"->"F"->"V")
#' @param Specimen The vector of specimen names
#' @param Test     The vector of test names
#' @param X        The vector of covariate names to separately extract data.
#'    For example, in PERCH data cleaning, \code{X = c("newSITE","CASECONT")}
#' @param Xval     The list of covariate values to stratify data.
#'    Each element corresponds to that in \code{X}. For example, in PERCH
#'    data cleaning, \code{Xval = list("02GAM","1")}
#' @param MeasDir  The directory to the data set (.csv)
#' @param PathCatDir The directory to the pathogen category specificaton file
#' (.csv)
#' @param extra_covariates The vector of covariate name for regression purposes.
#'   The default is NULL, which means no such covariates are necessary.
#' @param silent TURE/FALSE. If TRUE, the function will not print anything on
#' the screen. The default is TRUE.
#'
#' @return A list of data. Each element is either a measurement matrix, or
#' a vector of covariate values
#'
#' @export
#'
#'
#'
#'


extract_data_raw <-function(Pathogen,Specimen,Test,
                             X,Xval,
                             MeasDir,PathCatDir,
                             extra_covariates=NULL,
                             silent=TRUE){

  pathogen_type = read.csv(PathCatDir)
  rownames(pathogen_type) = pathogen_type[,1]
  typeOrder = order(pathogen_type[Pathogen,2])
  Pathogen = Pathogen[typeOrder]

  if (!silent){
    #show the ordered pathogen
    cat("Pathogens included:")
    print(t(rbind(Pathogen,as.character(pathogen_type[Pathogen,2]))))
  }

  datraw = read.csv(MeasDir)
  #clean column names if the column names start with "_":
  delete_start_with = function(s,vec){
    ind = grep(s,substring(vec,1,nchar(s)))
    old = vec[ind]
    vec[ind] = substring(old,nchar(s)+1)
    return(vec)
  }
  cleanName=delete_start_with("X_",names(datraw))
  dat0 = datraw
  colnames(dat0) = cleanName

  indX = 1:nrow(dat0)
  for (j in 1:length(X)){
    indX = indX[which(dat0[indX,X[j]]==Xval[[j]] & !is.na(dat0[indX,X[j]]))]
  }
  dat = dat0[indX,]

  pstGrid = apply(expand.grid(paste(Pathogen,"_",sep=""),Specimen,Test),
                  1,paste,collapse="")

  # the output will be for each (specimen, test) pair, with pathogens aligned
  stGrid = apply(expand.grid(Specimen,Test),1,paste,collapse="")

  pstTable = matrix(NA,length(Pathogen),length(stGrid))
  colnames(pstTable) = stGrid
  rownames(pstTable) = Pathogen

  for (i in 1:length(stGrid)){
    for (j in 1:length(Pathogen)){
      tempName = paste(Pathogen[j],stGrid[i],sep="_")
      if (tempName%in%cleanName){
        if (sum(is.na(dat[,tempName]))!=nrow(dat)){
          pstTable[j,i] = TRUE
        }
      }
    }
  }

  notindata = which(rowSums(!is.na(pstTable))==0)
  if (length(notindata)>0) {
    stop(Pathogen[notindata],"can't be found in the dataset!","\n")
  }
  naColumns = which(colSums(is.na(pstTable))==length(Pathogen))
  if (length(naColumns)==length(stGrid)) {
    stop("No test has available results on selected pathogens! Try other pathogens.","\n")
  } else{
    if (length(naColumns)==0){
      actualpstTable = pstTable
    } else{
      actualpstTable = pstTable[,-naColumns,drop=FALSE]
    }

    resdat = list()

    for (j in 1:ncol(actualpstTable)){
      if (sum(is.na(actualpstTable[,j]))==0){
        tempnm = paste(Pathogen,colnames(actualpstTable)[j],sep="_")
        resdat[[j]]=dat[,tempnm]
      } else{
        dftemp = as.data.frame(matrix(NA,nrow=nrow(dat),ncol=nrow(actualpstTable)))
        colnames(dftemp) = paste(Pathogen,colnames(actualpstTable)[j],sep="_")
        dftemp[,!is.na(actualpstTable[,j])] =
          dat[,paste(Pathogen[!is.na(actualpstTable[,j])],
                     colnames(actualpstTable)[j],sep="_")]
        resdat[[j]]=dftemp
      }
    }


    # get extra covariates used for regression purposes:
    if (!is.null(extra_covariates)){
      for (i in seq_along(extra_covariates)){
        if (!extra_covariates[i]%in% colnames(dat)){
          stop(extra_covariates[i]," is not in the data set!","\n")
        } else {
          resdat[[ncol(actualpstTable)+i]]=dat[,extra_covariates[i]]
        }
      }
    }
#    curr_len = length(resdat)
#     for ( j in seq_along(X)){
#       resdat[[curr_len+j]] = rep(Xval[[j]],nrow(dat))
#     }
#    names(resdat) = c(colnames(actualpstTable),extra_covariates,X)
    names(resdat) = c(colnames(actualpstTable),extra_covariates)
    return(resdat)
  }
}
