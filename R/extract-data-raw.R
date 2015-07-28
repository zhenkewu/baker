#' Convert raw data sheet to analyzable format
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
#' @param meas_dir  The directory to the data set (.csv)
#' @param extra_covariates The vector of covariate name for regression purposes.
#'   The default is NULL, which means no such covariate is necessary.
#' @param silent Default is \code{TRUE}: the function will not print anything on
#' the screen; otherwise specify as \code{FALSE}.
#' @param individual index for an individual to print his/her measurements
#'
#' @return A list of data. Each element is either a measurement matrix, or
#' a vector of covariate values
#'
#' @export

extract_data_raw <-function(Pathogen,Specimen,Test,
                             X,Xval,
                             meas_dir,
                             extra_covariates=NULL,
                             silent=TRUE,
                             individual=NULL){

#   #
#   # test:
#   #
#   
#   Pathogen = pathogen_BrS
#   Specimen = Specimen
#   Test = Test
#   X = c(X_strat,case_def)
#   Xval = append(X_strat_val,case_def_val)
#   meas_dir = meas_dir
#   PathCatDir = PathCatDir
#   extra_covariates = X_extra
#   silent=TRUE
#   individual=c(1,10)
#   
#   #
#   #
#   #
#   

  datraw = read.csv(meas_dir,header=TRUE,stringsAsFactors=FALSE)
  dat0   = datraw
  cleanName = colnames(dat0)

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
 
  # pathogens not in the data set:
  notindata = which(rowSums(!is.na(pstTable))==0) 
  if (length(notindata)>0) {
    stop("==",Pathogen[notindata],"can't be found in the dataset.",
         " delete this pathogen or put in measurements for this pathogen.","\n")
  }
  
  # (specimen,test) pairs that have no measurements on any pathogens:
  naColumns = which(colSums(is.na(pstTable))==length(Pathogen))
  if (length(naColumns)==length(stGrid)) {# if no (specimen,test) pair has measurements
    stop("==No (specimen, test) pair has available results on specified pathogens! Try other pathogens.==","\n")
  } else{
    if (length(naColumns)==0){
      actualpstTable = pstTable
    } else{
      actualpstTable = pstTable[,-naColumns,drop=FALSE]
    }

    resdat = list()

    for (j in 1:ncol(actualpstTable)){
      if (sum(is.na(actualpstTable[,j]))==0){# if all rows for this (pathogen, test) pair are in the data set:
        tempnm = paste(Pathogen,colnames(actualpstTable)[j],sep="_")
        resdat[[j]]=dat[,tempnm]
      } else{# if some pathogens are not in the data set:
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
          stop("==",extra_covariates[i]," is not in the data set. Delete this covariate.==","\n")
        } else {
          resdat[[ncol(actualpstTable)+i]]=dat[,extra_covariates[i]]
        }
      }
    }

    
    names(resdat) = c(colnames(actualpstTable),extra_covariates)
    
    #
    # display an individual's measurements: rows are pathogens, columns are
    # (specimen,test) pairs:
    #
    if (!is.null(individual)){
          if (nrow(actualpstTable)>0 & ncol(actualpstTable)>0){
            for (i in seq_along(individual)){
              dat_sheet           <- matrix(NA,nrow=nrow(actualpstTable),ncol=ncol(actualpstTable))
              rownames(dat_sheet) <- rownames(actualpstTable)
              colnames(dat_sheet) <- colnames(actualpstTable)
              
              for (p in 1:nrow(actualpstTable)){
                for (j in 1:ncol(actualpstTable)){
                  if (!is.na(actualpstTable[p,j])){
                    tmpname <- paste0(rownames(actualpstTable)[p],"_",colnames(actualpstTable)[j])
                    dat_sheet[p,j] <- resdat[[colnames(actualpstTable)[j]]][individual[i],tmpname]
                 }
               }
              }
              print(individual[i])
              print(dat_sheet)
           }
          }
    }
    #
    #
    #
    
    return(resdat)
  }
}

# # add the following if we want to organize bacteria first:
# @param PathCatDir The directory to the pathogen category specificaton file
# (.csv)
#  #add right after the starting{
#   pathogen_type = read.csv(PathCatDir)
#   rownames(pathogen_type) = pathogen_type[,1]
#   typeOrder = order(pathogen_type[Pathogen,2])
#   Pathogen = Pathogen[typeOrder]
#   if (!silent){
#     #show the ordered pathogen
#     cat("Pathogens included:")
#     print(t(rbind(Pathogen,as.character(pathogen_type[Pathogen,2]))))
#   }
