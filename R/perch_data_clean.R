#' Clean PERCH data
#'
#' It deletes cases with two positives for BCX measures. To order cases
#' and controls separately according to site and enrollment date, we suggest put
#' \code{c("newSITE","ENRLDATE")} in \code{X_extra} and \code{X_order_obs} as listed
#' below.
#'
#' @param clean_options The list of options for cleaning PERCH data.
#' The specific elements are as follows:
#' \itemize{
#' \item{\code{case_def}}{: variable name for case definition;}
#' \item{\code{case_def_val}}{: The value of the case-definition variable;}
#' \item{\code{ctrl_def}}{: variable name for control definition;}
#' \item{\code{ctrl_def_val}}{: value of the control-definition variable;}
#' \item{\code{X_strat}}{: A vector of strings, each defining the variables used to
#' stratify the data;}
#' \item{\code{X_strat_val}}{: A list of actual values that X_strat should take to
#' get stratified data sets;}
#' \item{\code{pathogen_BrS}}{: The vector of pathogen names (arbitrary order).
#' It has to be a subset of pathogen category information in the \code{PathCatDir}.}
#' \item{\code{X_extra}}{: A vector of strings, each being the covariate name that one
#' wants to have in the data set to be analyzed;}
#' \item{\code{X_order_obs}}{: A vector of strings, each being the covariate name that one
#' wants to order the rows of the observations. This helps observations with same
#' site to be put together, or helps order dates from earlier to later;}
#' \item{\code{RawMeasDir}}{: The file path to the raw data set;}
#' \item{\code{write_newSite}}{: Must set to \code{TRUE} if the raw data set is changed;}
#' \item{\code{newSite_write_Dir}}{: The file path where a new/cleaned data set for actual
#'  analyses will be written;}
#' \item{\code{MeasDir}}{: The file path to the cleaned data set;}
#' \item{\code{PathCatDir}}{: The file path to the pathogen category list (.csv). This list
#' should be as complete as possible to display all pathogens used in an actual
#' analysis;}
#' \item{\code{allow_missing}}{: TRUE for the model to using an observation that has either
#' BrS missing, or SS missing, or both missing. It has to be set to \code{TRUE}
#' in order to get all the BCX information available, albeit some cases have missing
#' NPPCR data;}
#'}
#'
#' @return A List: \code{list(Mobs,Y,X,JSS,pathogen_MSS_ordered,pathogen_cat)}, or
#' with additional \code{JSSonly, pathogen_SSonly_cat} if silver-
#' standard only pathogens are supplied.
#' \itemize{
#' \item \code{Mobs} A list of bronze- (\code{MBS}), silver- (\code{MSS}),
#' and gold-standard (\code{MGS}, if available) measurements. Here if all
#' pathogens have BrS measures, MSS has the same number of columns as in MBS;
#' if some pathogens only have SS measures, then MSS will have extra columns
#' for these measurements;
#' \item \code{Y} Binary indicator for cases (1) and controls (0);
#' \item \code{X} Data frame of covariates for cases and controls. The names
#' are specified in \code{X_extra};
#' \item \code{JSS} Number of pathogens with both silver- and bronze-standard
#' data;
#' \item \code{pathogen_MSS_ordered} Ordered vector of pathogen names. Pathogens
#' with both SS and BrS pathogens are ordered first and then those with only BrS
#' measurements. Pathogens with only silver-standard measures are not included.
#' \item \code{pathogen_cat} Pathogen categories ordered according to
#' \code{pathogen_MSS_ordered}.
#' \item \code{JSSonly} Number of pathogens with only silver-standard data;
#' \item \code{pathogen_SSonly_cat} Category of pathogens with only silver-standard
#' data.
#' }
#' The function does not order silver-standard data only pathogens.
#' @import lubridate
#' @export


perch_data_clean <- function(clean_options){

  case_def     <- clean_options$case_def
  case_def_val <- clean_options$case_def_val
  ctrl_def     <- clean_options$ctrl_def
  ctrl_def_val <- clean_options$ctrl_def_val
  X_strat      <- clean_options$X_strat
  X_strat_val  <- clean_options$X_strat_val
  pathogen_BrS     <- clean_options$pathogen_BrS
  X_extra      <- clean_options$X_extra
  RawMeasDir   <- clean_options$RawMeasDir
  write_newSite<- clean_options$write_newSite
  newSite_write_Dir <- clean_options$newSite_write_Dir
  MeasDir           <- clean_options$MeasDir
  PathCatDir        <- clean_options$PathCatDir
  X_order_obs       <- clean_options$X_order_obs


  # if no silver-only pathogen measurements
  if (!is.null(clean_options$pathogen_SSonly)){
          pathogen_SSonly   <- clean_options$pathogen_SSonly
  }

  # combine two sub-sites:
  # 06NTH and 07STH  --> THA,
  # 08MBA and 09DBA  --> BAN:
  PERCH_data_with_newSITE <- combine_subsites(RawMeasDir,
                               subsites_list = list(c("06NTH","07STH"),
                                                    c("08MBA","09DBA")),
                               newsites_vec  = c("06THA","07BAN"))

  if (write_newSite){
      write.csv(PERCH_data_with_newSITE, newSite_write_Dir)
  }

  # list the pathogen categories:
  pathogen_cat_lookup <- read.csv(PathCatDir)

  Specimen  <- c("NP","B")
  Test      <- c("PCR","CX")

  # extract_data_raw() will order pathogen_BrS according to B->F->V:
  if (is.null(X_strat) && is.null(X_strat_val)){
    # no stratification:
    datacase <- extract_data_raw(pathogen_BrS,Specimen,Test,
                                 c(case_def),list(case_def_val),
                                 extra_covariates = X_extra,
                                 MeasDir,PathCatDir,silent=TRUE)

    datactrl <- extract_data_raw(pathogen_BrS,Specimen,Test,
                                 c(ctrl_def),list(ctrl_def_val),
                                 extra_covariates = X_extra,
                                 MeasDir,PathCatDir,silent=TRUE)
  }else{
    # stratify by levels of X_strat:
    datacase <- extract_data_raw(pathogen_BrS,Specimen,Test,
                                 c(X_strat,case_def),
                                 append(X_strat_val,case_def_val),
                                 extra_covariates = X_extra,
                                 MeasDir,PathCatDir,silent=TRUE)

    datactrl <- extract_data_raw(pathogen_BrS,Specimen,Test,
                                 c(X_strat,ctrl_def),
                                 append(X_strat_val,ctrl_def_val),
                                 extra_covariates = X_extra,
                                 MeasDir,PathCatDir,silent=TRUE)
  }

  #write.csv(datacase,"C:/package_test/datacase.csv")
  #write.csv(datactrl,"C:/package_test/datactrl.csv")

  # get pathogens that have many BcX measurements:
  SS_index  <- which(colMeans(is.na(datacase$BCX))<.9)
  cat("==Pathogens with both blood culture and NPPCR measures:==","\n",
      pathogen_BrS[SS_index],"\n")
  JSS       <- length(SS_index)
  JBrS      <- length(pathogen_BrS)
  # order pathogens based on BcX+NPPCR --> NPPCR:
  MSS_avail_order_index <- c(SS_index,(1:JBrS)[-SS_index])
  # (1:JBrS) necessary to avoid arranging SS_only pathogens.

  if (!is.null(pathogen_SSonly)){
    # for silver-standard only Bacteria:
    datacase_SSonly <- extract_data_raw(pathogen_SSonly,Specimen,Test,
                                        c(X_strat,case_def),append(X_strat_val,case_def_val),
                                        extra_covariates = X_extra,
                                        MeasDir,PathCatDir,silent=TRUE)
    SSonly_index <- which(colMeans(is.na(datacase_SSonly$BCX))<.9)
    cat("==Pathogens with ONLY blood culture measure:==","\n",
        pathogen_SSonly[SSonly_index],"\n")
    JSSonly       <- length(SSonly_index)
  }

  # get case/control indices who have complete observations on
  # NPPCR & BCX(cases), or NPPCR(ctrls):
  complete_case_index <- which(rowMeans(is.na(datacase$NPPCR))==0 &
                               rowMeans(is.na(datacase$BCX[,SS_index]))==0)
  case_index_complete_NP_incomp_BCX <- which(rowMeans(is.na(datacase$NPPCR))==0 &
                                             rowMeans(is.na(datacase$BCX[,SS_index]))>0)
  case_index_complete_BCX_incomp_NP <- which(rowMeans(is.na(datacase$NPPCR))>0 &
                                               rowMeans(is.na(datacase$BCX[,SS_index]))==0)
  case_index_incomplete             <- which(rowMeans(is.na(datacase$NPPCR))> 0 &
                                               rowMeans(is.na(datacase$BCX[,SS_index]))>0)
  complete_ctrl_index <- which(rowMeans(is.na(datactrl$NPPCR))==0)

  # print incomplete indices:
  if (length(case_index_complete_NP_incomp_BCX)>0){
    cat("==Case indices with complete NP, incomplete BCX:==","\n",
        case_index_complete_NP_incomp_BCX,"\n")
    print(data.frame(datacase$patid,datacase$NPPCR)[case_index_complete_NP_incomp_BCX,])
    cat("==","\n")
    print(data.frame(datacase$patid,datacase$BCX)[case_index_complete_NP_incomp_BCX,])
  }

  if (length(case_index_complete_BCX_incomp_NP)>0){
    cat("==Case indices with complete BCX, incomplete NP:==","\n",
        case_index_complete_BCX_incomp_NP,"\n")
    print(data.frame(datacase$patid,datacase$NPPCR)[case_index_complete_BCX_incomp_NP,])
    cat("==","\n")
    print(data.frame(datacase$patid,datacase$BCX)[case_index_complete_BCX_incomp_NP,])
  }

  if (length(case_index_incomplete)>0){
    cat("==Case indices with incomplete NP, incomplete BCX:==","\n",
        case_index_incomplete,"\n")
    print(data.frame(datacase$patid,datacase$BCX)[case_index_incomplete,])
    cat("==","\n")
    print(data.frame(datacase$patid,datacase$NPPCR)[case_index_incomplete,])
  }

  if (!is.null(clean_options$allow_missing) &&
                clean_options$allow_missing==TRUE){
    complete_case_index <- 1:nrow(datacase$NPPCR)
    complete_ctrl_index <- 1:nrow(datactrl$NPPCR)
  }

  # actual numbers of complete cases/controls:
  Nd  <- length(complete_case_index)
  Nu  <- length(complete_ctrl_index)

  M_NPPCR    <- rbind(datacase$NPPCR[complete_case_index,],
                      datactrl$NPPCR[complete_ctrl_index,])
  Y          <- c(rep(1,Nd),rep(0,Nu))
  M_BCX_ctrl <- as.data.frame(matrix(NA,nrow=Nu,ncol=length(pathogen_BrS)))
  colnames(M_BCX_ctrl) <- paste(pathogen_BrS,"BCX",sep="_")
  M_BCX     <- rbind(datacase$BCX[complete_case_index,],M_BCX_ctrl)
  if (!is.null(pathogen_SSonly)){
    M_BCX_SSonly_ctrl <- as.data.frame(matrix(NA,nrow=Nu,ncol=length(pathogen_SSonly)))
    colnames(M_BCX_SSonly_ctrl) <- paste(pathogen_SSonly,"BCX",sep="_")
    M_BCX_SSonly      <- rbind(datacase_SSonly$BCX[complete_case_index,],M_BCX_SSonly_ctrl)
  }

  # get data sets as indicated by complete_case_index:
  datcase_X_extra_subset <- lapply(datacase[X_extra],"[",c(complete_case_index))
  datctrl_X_extra_subset <- lapply(datactrl[X_extra],"[",c(complete_ctrl_index))

  # combine specimen measurements, and covariates into data frames
  if (is.null(pathogen_SSonly)){
    # if no SS only pathogen is supplied:
    datobs <- data.frame(M_NPPCR,M_BCX,
                         Y = c(rep(1,Nd),rep(0,Nu)),
                         rbind(do.call(data.frame,datcase_X_extra_subset),
                               do.call(data.frame,datctrl_X_extra_subset)))
  } else {
    # if there exists SS only pathogen:
    datobs <- data.frame(M_NPPCR,cbind(M_BCX,M_BCX_SSonly),
                         Y = c(rep(1,Nd),rep(0,Nu)),
                         rbind(do.call(data.frame,datcase_X_extra_subset),
                               do.call(data.frame,datctrl_X_extra_subset)))

  }
  # 1. if ENRLDATE is in the dataset, transform the date format
  if ("ENRLDATE" %in% X_extra){
    Rdate.case <- as.Date(datcase_X_extra_subset$ENRLDATE, "%d%B%Y")
    Rdate.ctrl <- as.Date(datctrl_X_extra_subset$ENRLDATE, "%d%B%Y")

    uniq.month.case <- unique(paste(month(Rdate.case),year(Rdate.case),sep="-"))
    uniq.month.ctrl <- unique(paste(month(Rdate.ctrl),year(Rdate.ctrl),sep="-"))

    #symm.diff.dates <- as.set(uniq.month.case)%D%as.set(uniq.month.ctrl)
    #if (length(symm.diff.dates)!=0){
    #  cat("Cases and controls have different enrollment months:","\n")
    #  print(symm.diff.dates)
    #}
    datobs$ENRLDATE <- c(Rdate.case, Rdate.ctrl)
  }
  if ("patid" %in% X_extra){
    datobs$patid <- as.character(datobs$patid)
  }
  if ("newSITE" %in% X_extra){
    datobs$newSITE <- as.character(datobs$newSITE)
  }
  #     sort_data_frame <- function(x, decreasing=FALSE, by=1, ... ){
  #       f <- function(...) order(...,decreasing=decreasing)
  #       i <- do.call(f,x[by])
  #       x[i,,drop=FALSE]
  #     }
  #   # order observations in cases/controls defined by the order_obs variable:
  if (!is.null(X_order_obs)){
    # control first; then by X_order_obs:
    form_ord <- as.formula(paste0("~",paste(c("-Y",X_order_obs),collapse="+")))
    datobs <- sort_data_frame(form_ord,datobs)
  } else{
    form_ord <- as.formula("~-Y")
    datobs <- sort_data_frame(form_ord,datobs)
  }

  #datobs[datobs$patid=="G00493",]

  # create the list of measurements;
  # look for pathogens that have more than one positives in BCX:
  if (!is.null(pathogen_SSonly)){
    # if some SSonly pathogen is supplied:
    cat("==Total positives in silver-standard:==","\n")
    SS_index_all <- sapply(paste(c(pathogen_BrS[SS_index],pathogen_SSonly),"BCX",sep="_"),
                           grep,names(datobs))

    print(table(rowSums(datobs[1:Nd,SS_index_all])))
    BCX_more_than_one_index <- which(rowSums(datobs[1:Nd,SS_index_all])>1)
    if (length(BCX_more_than_one_index)>0){
      cat("==Removed case(s) with >1 positive in BCX:==","\n")
      print(datobs$patid[BCX_more_than_one_index])
      print(datobs[BCX_more_than_one_index,SS_index_all])


      Mobs   <- list(MBS = datobs[-BCX_more_than_one_index, grep("_NPPCR",names(datobs))[MSS_avail_order_index]],
                     MSS = datobs[-BCX_more_than_one_index, c(paste(pathogen_BrS,"BCX",sep="_")[MSS_avail_order_index],
                                                              paste(pathogen_SSonly,"BCX",sep="_"))],
                     MGS = NA)
      Y      <- datobs$Y[-BCX_more_than_one_index]
      X      <- datobs[-BCX_more_than_one_index,X_extra]
      Nd     <- Nd - length(BCX_more_than_one_index)
    } else {
      Mobs   <- list(MBS = datobs[,grep("_NPPCR",names(datobs))[MSS_avail_order_index]],
                     MSS = datobs[,c(paste(pathogen_BrS,"BCX",sep="_")[MSS_avail_order_index],
                                     paste(pathogen_SSonly,"BCX",sep="_"))],
                     MGS = NA)
      Y      <- datobs$Y
      X      <- datobs[,X_extra]
    }
  } else{

    # if no SSonly pathogen is supplied:
    cat("==Total positives in silver-stand:==","\n")
    SS_index_all <- sapply(paste(c(pathogen_BrS[SS_index]),"BCX",sep="_"),grep,names(datobs))

    print(table(rowSums(datobs[1:Nd,SS_index_all])))
    BCX_more_than_one_index <- which(rowSums(datobs[1:Nd,SS_index_all])>1)
    if (length(BCX_more_than_one_index)>0){
      cat("==Removed case(s) with >1 positive in BCX:==","\n")
      print(datobs$patid[BCX_more_than_one_index])
      print(datobs[BCX_more_than_one_index,SS_index_all])


      Mobs   <- list(MBS = datobs[-BCX_more_than_one_index, grep("_NPPCR",names(datobs))[MSS_avail_order_index]],
                     MSS = datobs[-BCX_more_than_one_index, c(paste(pathogen_BrS,"BCX",sep="_")[MSS_avail_order_index])],
                     MGS = NA)
      Y      <- datobs$Y[-BCX_more_than_one_index]
      X      <- datobs[-BCX_more_than_one_index,X_extra]
      Nd     <- Nd - length(BCX_more_than_one_index)
    } else {
      Mobs   <- list(MBS = datobs[,grep("_NPPCR",names(datobs))[MSS_avail_order_index]],
                     MSS = datobs[,c(paste(pathogen_BrS,"BCX",sep="_")[MSS_avail_order_index])],
                     MGS = NA)
      Y      <- datobs$Y
      X      <- datobs[,X_extra]
    }
  }

  # order the whole pathogen list (excluing SS-only pathogens),
  # so that those have SS data comes first
  pathogen_MSS_ordered <- pathogen_BrS[MSS_avail_order_index]

  # get the list of pathogen categories for each of the above ordered pathogens:
  path_cat_ind <- sapply(1:JBrS,function(i)
         which(pathogen_cat_lookup$X==pathogen_MSS_ordered[i]))

  if (!is.null(pathogen_SSonly)){
    pathogen_SSonly_cat <- pathogen_cat_lookup[sapply(1:JSSonly,function(i)
                         which(pathogen_cat_lookup$X==pathogen_SSonly[i])),]
    list(Mobs = Mobs,
         Y = Y,
         X = X,
         JSS  = JSS,
         pathogen_MSS_ordered = pathogen_MSS_ordered,
         pathogen_cat         = pathogen_cat_lookup[path_cat_ind,],
         JSSonly = JSSonly,
         pathogen_SSonly_cat = pathogen_SSonly_cat)
  } else{
    list(Mobs = Mobs,
         Y = Y,
         X = X,
         JSS  = JSS,
         pathogen_MSS_ordered = pathogen_MSS_ordered,
         pathogen_cat         = pathogen_cat_lookup[path_cat_ind,])
  }

}
