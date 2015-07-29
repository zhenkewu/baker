#' Clean PERCH data
#'
#' \code{clean_perch_data} transforms a raw data table (row for subject, column
#' for variable - usually \code{\{pathogen name\}_\{specimen\}} and other covariate
#' names) into a list. It is specific for PERCH data format.
#'
#' @details It deletes cases (\code{Y==1}) having two positives 
#' for BCX measures. We suggest put \code{c("newSITE","ENRLDATE")} in 
#' \code{X_extra} and \code{X_order_obs} to order cases and controls separately 
#' according to site and enrollment date. In current implementation, 
#' the raw data must have both BrS and SS measurements. 
#' @param clean_options The list of options for cleaning PERCH data.
#' Its elements are defined as follows:
#' \itemize{
#' \item{\code{raw_meas_dir}}{: The file path to the raw data;}
#' \item{\code{case_def}}{: variable name in raw data for case definition;}
#' \item{\code{case_def_val}}{: The value for case definition;}
#' \item{\code{ctrl_def}}{: variable name in raw data for control definition;}
#' \item{\code{ctrl_def_val}}{: The value for control definition;}
#' \item{\code{X_strat}}{: A vector of variable names for stratifying the data
#' to perform SEPARATE analyses;}
#' \item{\code{X_strat_val}}{: A list of values for \code{X_strat}. The output
#' data will only correspond to those with \code{identical(X_strat,X_strat_val)==TRUE}.
#' To perform analysis on a single site, say \code{"02GAM"}, use \code{X_strat="newSITE"} and 
#' \code{X_strat_val=list("02GAM")};}
#' \item{\code{pathogen_BrS_anyorder}}{: The vector of pathogen names (arbitrary 
#' order) that have bronze-standard (BrS) measurments (cf. Wu et al. (2015) for 
#' definitions and examplels of BrS, SS, and GS). It has to be a subset of 
#' pathogens listed in taxonomy information at the file path \code{patho_taxo_dir}.}
#' \item{\code{pathogen_SSonly}}{: A vector of pathogens that only have SS data;}
#' \item{\code{X_extra}}{: A vector of covariate names for regression
#' or visualization;}
#' \item{\code{X_order_obs}}{: A vector of variable names for ordering observations. 
#' For example, it can include site names or enrollment dates. It must be a 
#' subset of X_extra;}
#' \item{\code{patho_taxo_dir}}{: The file path to the pathogen category or taxonomy 
#' information (.csv). The information should be as complete as possible to 
#' display all pathogens considered in an actual study;}
#' \item{\code{date_formats}}{possible formats of date; default is 
#' \code{c("\%d\%B\%Y","\%d\%B\%y")}.
#' See \link[lubridate]{parse_date_time} for a complete list of date formats.}
#' \item{\code{allow_missing}}{: \code{TRUE} for using an observation that has 
#' either BrS missing, or SS missing. Set it to \code{TRUE} if we want to use 
#' the SS information from some cases who missed BrS measurements. 
#' In other words, all the subjects' data will be used if \code{allow_missing} is 
#' set to \code{TRUE}.}
#' \item{\code{extra_meas_nm}}{: a list of (pathogen,specimen,test) names, each of which
#' is considered extra measurements informative for etiology.}
#'}
#'
#' @return A List: \code{list(Mobs,Y,X,JSS,pathogen_BrS_ordered_by_MSS,pathogen_BrS_cat)}, or
#' with additional \code{JSSonly, pathogen_SSonly_cat} if silver-
#' standard only pathogens are supplied.
#' \itemize{
#' \item \code{Mobs} A list of bronze- (\code{MBS}), silver- (\code{MSS}),
#' and gold-standard (\code{MGS}, if available) measurements. Here if all
#' pathogens have BrS measures, MSS has the same number of columns as in MBS;
#' if some pathogens only have SS measures, then MSS will have extra columns;
#' \item \code{Y} 1 for case; 0 for control;
#' \item \code{X} Data frame of covariates for cases and controls. The covariate
#' names are specified in \code{X_extra};
#' \item \code{JSS} Number of pathogens having both silver- and bronze-standard
#' data;
#' \item \code{pathogen_BrS_ordered_by_MSS} Ordered vector of pathogen names. Pathogens
#' with both SS and BrS pathogens are ordered first and then those with only BrS
#' measurements. Note that for a pathogen name vector of arbitrary order 
#' (\code{pathogen_BrS_anyorder}), this function just picks out those pathogens
#' with BrS+SS measures and puts them at the front. Other pathogens with only 
#' BrS measures are not reordered. Pathogens with only silver-standard measures 
#' are not included.
#' \item \code{pathogen_BrS_cat} Pathogen categories ordered according to
#' \code{pathogen_BrS_ordered_by_MSS}.
#' \item \code{JSSonly} Number of pathogens with only silver-standard measures;
#' \item \code{pathogen_SSonly_cat} Category of pathogens with only silver-standard
#' data.
#' \item \code{extra_Mobs} extra measurements as requested by \code{extra_meas_nm}.
#' }
#' This function does not re-order pathogens that only have silver-standard data.
#' 
#' @seealso \link[lubridate]{parse_date_time}
#' @export

clean_perch_data <- function(clean_options){

  raw_meas_dir <- clean_options$raw_meas_dir
  case_def     <- clean_options$case_def
  case_def_val <- clean_options$case_def_val
  ctrl_def     <- clean_options$ctrl_def
  ctrl_def_val <- clean_options$ctrl_def_val
  X_strat      <- clean_options$X_strat
  X_strat_val  <- clean_options$X_strat_val
  pathogen_BrS_anyorder     <- clean_options$pathogen_BrS_anyorder
  X_extra      <- clean_options$X_extra
  patho_taxo_dir        <- clean_options$patho_taxo_dir
  X_order_obs       <- clean_options$X_order_obs
  extra_meas_nm  <- clean_options$extra_meas_nm
  
  # clean dates: use specified date format; if not specified, try the formats
  # specified in the "else" sub-clause:
  if (!is.null(clean_options$date_formats)){
    date_formats      <- clean_options$date_formats
  }else{
    date_formats      <- c("%d%B%Y","%d%B%y")
  }

  # if there are silver-only pathogen measurements:
  if (!is.null(clean_options$pathogen_SSonly)){
          pathogen_SSonly   <- clean_options$pathogen_SSonly
  }

#   # check if the data specified in 'raw_meas_dir' has been cleaned by the pacakge:
#   # Here "prepared" means revoming X_ from the colmn names and combined Bangladesh and
#   # Thailand subsites. The final output of this whole function will be called "cleaned".
#   raw_data <- read.csv(raw_meas_dir)
#   if (is.null(attr(raw_data,"prepared_by_package"))){
#     do_prepare <- TRUE
#   } else{
#     if (attr(raw_data,"prepared_by_package")== TRUE){
#       do_prepare <- FALSE
#       meas_dir <- raw_meas_dir
#     } else{
#       do_prepare <- TRUE
#     }
#   }
#   rm(raw_data)
#   
#   if (do_prepare){
    # combine two sub-sites:
    # 06NTH and 07STH  --> THA,
    # 08MBA and 09DBA  --> BAN:
    PERCH_data_with_newSITE <- clean_combine_subsites(raw_meas_dir,
                                 subsites_list = list(c("06NTH","07STH"),
                                                      c("08MBA","09DBA")),
                                 newsites_vec  = c("06THA","07BAN"))
    
    cleanName  <- delete_start_with("X_",names(PERCH_data_with_newSITE))
    colnames(PERCH_data_with_newSITE) <- cleanName
    
    #attr(PERCH_data_with_newSITE,"prepared_by_package") <- TRUE
    
    # write cleaned data into working directory:
    meas_dir <- file.path(dirname(raw_meas_dir),paste0("prepared_",basename(raw_meas_dir)))
    write.csv(PERCH_data_with_newSITE, meas_dir, row.names=FALSE)
    rm(PERCH_data_with_newSITE)
  #}
  
  # create the pathogen category lookup table:
  pathogen_cat_lookup <- read.csv(patho_taxo_dir,stringsAsFactors=FALSE)

  Specimen  <- c("NP","B")
  Test      <- c("PCR","CX")
  #Specimen  <- c("NP","B","IS","LA","PF")
  #Test      <- c("PCR","CX","CX2")

  # extract_data_raw() will NOT order pathogen_BrS_anyorder according to B->F->V:
  if (is.null(X_strat) && is.null(X_strat_val)){
    # no stratification:
    X_strat_nm_tmp  <- c(case_def)
    X_strat_val_tmp <- list(case_def_val)
  }else{
    # stratify by levels of X_strat:
    X_strat_nm_tmp  <- c(X_strat, case_def)
    X_strat_val_tmp <- list(X_strat_val, case_def_val)
  }
  datacase <- extract_data_raw(pathogen_BrS_anyorder,Specimen,Test,
                               X_strat_nm_tmp,
                               X_strat_val_tmp,
                               extra_covariates = X_extra,
                               meas_dir,silent=TRUE)
  
  datactrl <- extract_data_raw(pathogen_BrS_anyorder,Specimen,Test,
                               X_strat_nm_tmp,
                               X_strat_val_tmp,
                               extra_covariates = X_extra,
                               meas_dir,silent=TRUE)
  
  
  if (!is.null(extra_meas_nm)){
      # read in extra measurements, e.g., 2nd or 3rd bronze-standard:
      prepared_data <- read.csv(meas_dir,header=TRUE,stringsAsFactors=FALSE)
      
      indX = 1:nrow(prepared_data)
      for (j in 1:length(X_strat_nm_tmp)){
        indX = indX[which(prepared_data[indX,X_strat_nm_tmp[j]]==X_strat_val_tmp[[j]] & 
                            !is.na(prepared_data[indX,X_strat_nm_tmp[j]]))]
      }
      extracted_prepared_data = prepared_data[indX,]
      
      extra_Mobs <- list()
      for (i in seq_along(extra_meas_nm)){
        extra_Mobs <- append(extra_Mobs, list(extracted_prepared_data[extra_meas_nm[[i]]]))
      }
      names(extra_Mobs) <- names(extra_meas_nm)
  }
  
  #write.csv(datacase,"C:/package_test/datacase.csv")
  #write.csv(datactrl,"C:/package_test/datactrl.csv")

  flag_BCX_in_datacase <- "BCX"%in%names(datacase)
  # if there is BCX measurements:
  if (flag_BCX_in_datacase){
    # get pathogens that have many BcX measurements:
    SS_index  <- which(colMeans(is.na(datacase$BCX))<.9)
    cat("==Pathogens with both 'blood culture' and 'NPPCR' measures:==","\n",
        pathogen_BrS_anyorder[SS_index],"\n")
    JSS       <- length(SS_index)  
    JBrS      <- length(pathogen_BrS_anyorder)
    # order pathogens based on BcX+NPPCR --> NPPCR:
    MSS_avail_order_index <- c(SS_index,(1:JBrS)[-SS_index])
    # (1:JBrS) necessary to avoid arranging SS_only pathogens.
  }else{
#     stop("==This dataset does not have pathogens that have both BrS and SS measurements.
#          Please contact maintainer of this package for an update. Thanks. ==")
  
    JBrS      <- length(pathogen_BrS_anyorder)
    JSS       <- 0
    # order pathogens based on BcX+NPPCR --> NPPCR:
    MSS_avail_order_index <- 1:JBrS
    datacase$BCX <- datacase$NPPCR+NA
    colnames(datacase$BCX) <- paste(pathogen_BrS_anyorder,"BCX",sep="_")
    datactrl$BCX <- datactrl$NPPCR+NA
    colnames(datactrl$BCX) <- paste(pathogen_BrS_anyorder,"BCX",sep="_")
  }
  
  if (!is.null(pathogen_SSonly)){
    # for silver-standard only Bacteria:
    datacase_SSonly <- extract_data_raw(pathogen_SSonly,Specimen,Test,
                                        c(X_strat,case_def),append(X_strat_val,case_def_val),
                                        extra_covariates = X_extra,
                                        meas_dir,silent=TRUE)
    if (length(pathogen_SSonly)>1){
      SSonly_index <- which(colMeans(is.na(datacase_SSonly$BCX))<.9)
    }else{
      SSonly_index <- which(mean(is.na(datacase_SSonly$BCX))<.9)
    }
    cat("==Pathogens with ONLY blood culture measure:==","\n",
        pathogen_SSonly[SSonly_index],"\n")
    JSSonly       <- length(SSonly_index)
  }

  # get case/control indices who have complete observations on
  # NPPCR & BCX(cases), or NPPCR(ctrls):

  if (flag_BCX_in_datacase){
    complete_case_index <- which(rowMeans(is.na(datacase$NPPCR))==0 &
                                   rowMeans(is.na(datacase$BCX[,SS_index,drop=FALSE]))==0)
    case_index_complete_NP_incomp_BCX <- which(rowMeans(is.na(datacase$NPPCR))==0 &
                                               rowMeans(is.na(datacase$BCX[,SS_index,drop=FALSE]))>0)
    case_index_complete_BCX_incomp_NP <- which(rowMeans(is.na(datacase$NPPCR))>0 &
                                                 rowMeans(is.na(datacase$BCX[,SS_index,drop=FALSE]))==0)
    case_index_incomplete             <- which(rowMeans(is.na(datacase$NPPCR))> 0 &
                                                 rowMeans(is.na(datacase$BCX[,SS_index,drop=FALSE]))>0)
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
  }else{
    complete_case_index <- which(rowMeans(is.na(datacase$NPPCR))==0)
    case_index_incomplete             <- which(rowMeans(is.na(datacase$NPPCR))> 0)
    # print incomplete indices:

    if (length(case_index_incomplete)>0){
      cat("==Case indices with incomplete NP:==","\n",
          case_index_incomplete,"\n")
      print(data.frame(datacase$patid,datacase$NPPCR)[case_index_incomplete,])
    }
  }
  complete_ctrl_index <- which(rowMeans(is.na(datactrl$NPPCR))==0)

  

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
  M_BCX_ctrl <- as.data.frame(matrix(NA,nrow=Nu,ncol=length(pathogen_BrS_anyorder)))
  colnames(M_BCX_ctrl) <- paste(pathogen_BrS_anyorder,"BCX",sep="_")
  M_BCX     <- rbind(datacase$BCX[complete_case_index,],M_BCX_ctrl)
  if (!is.null(pathogen_SSonly)){
    M_BCX_SSonly_ctrl <- as.data.frame(matrix(NA,nrow=Nu,ncol=length(pathogen_SSonly)))
    colnames(M_BCX_SSonly_ctrl) <- paste(pathogen_SSonly,"BCX",sep="_")
    if (length(pathogen_SSonly)>1){
      M_BCX_SSonly      <- rbind(datacase_SSonly$BCX[complete_case_index,],M_BCX_SSonly_ctrl)
    }else{
      M_BCX_SSonly_case <- as.matrix(datacase_SSonly$BCX[complete_case_index],ncol=1)
      colnames(M_BCX_SSonly_case) <- colnames(M_BCX_SSonly_ctrl)
      M_BCX_SSonly      <- rbind(M_BCX_SSonly_case,M_BCX_SSonly_ctrl)
    }
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
    #Rdate.case <- as.Date(datcase_X_extra_subset$ENRLDATE, "%d%B%Y")
    #Rdate.ctrl <- as.Date(datctrl_X_extra_subset$ENRLDATE, "%d%B%Y")
    Rdate.case <- lubridate::parse_date_time(datcase_X_extra_subset$ENRLDATE, date_formats)
    Rdate.ctrl <- lubridate::parse_date_time(datctrl_X_extra_subset$ENRLDATE, date_formats)

    uniq.month.case <- unique(paste(lubridate::month(Rdate.case),lubridate::year(Rdate.case),sep="-"))
    uniq.month.ctrl <- unique(paste(lubridate::month(Rdate.ctrl),lubridate::year(Rdate.ctrl),sep="-"))

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
    if (flag_BCX_in_datacase){
      SS_index_all <- sapply(paste(c(pathogen_BrS_anyorder[SS_index],pathogen_SSonly),"BCX",sep="_"),
                             grep,names(datobs))
    }else{
      SS_index_all <- sapply(paste(c(pathogen_SSonly),"BCX",sep="_"),
                             grep,names(datobs))
    }

    print(table(rowSums(datobs[1:Nd,SS_index_all,drop=FALSE])))
    BCX_more_than_one_index <- which(rowSums(datobs[1:Nd,SS_index_all,drop=FALSE])>1)
    if (length(BCX_more_than_one_index)>0){
      cat("==Removed case(s) with >1 positive in BCX:==","\n")
      print(datobs$patid[BCX_more_than_one_index])
      print(datobs[BCX_more_than_one_index,SS_index_all])


      Mobs   <- list(MBS = datobs[-BCX_more_than_one_index, grep("_NPPCR",names(datobs))[MSS_avail_order_index]],
                     MSS = datobs[-BCX_more_than_one_index, c(paste(pathogen_BrS_anyorder,"BCX",sep="_")[MSS_avail_order_index],
                                                              paste(pathogen_SSonly,"BCX",sep="_"))],
                     MGS = NA)
      Y      <- datobs$Y[-BCX_more_than_one_index]
      X      <- datobs[-BCX_more_than_one_index,X_extra]
      Nd     <- Nd - length(BCX_more_than_one_index)
    } else {
      Mobs   <- list(MBS = datobs[,grep("_NPPCR",names(datobs))[MSS_avail_order_index]],
                     MSS = datobs[,c(paste(pathogen_BrS_anyorder,"BCX",sep="_")[MSS_avail_order_index],
                                     paste(pathogen_SSonly,"BCX",sep="_"))],
                     MGS = NA)
      Y      <- datobs$Y
      X      <- datobs[,X_extra]
    }
  } else{# if no SSonly pathogen is supplied:
    # if BCX is available on some pathogens:
    if (flag_BCX_in_datacase){
        cat("==Total positives in silver-stand:==","\n")
        SS_index_all <- sapply(paste(c(pathogen_BrS_anyorder[SS_index]),"BCX",sep="_"),grep,names(datobs))
    
        print(table(rowSums(datobs[1:Nd,SS_index_all,drop=FALSE])))
        BCX_more_than_one_index <- which(rowSums(datobs[1:Nd,SS_index_all,drop=FALSE])>1)
        if (length(BCX_more_than_one_index)>0){
          cat("==Removed case(s) with >1 positive in BCX:==","\n")
          print(datobs$patid[BCX_more_than_one_index])
          print(datobs[BCX_more_than_one_index,SS_index_all])
    
    
          Mobs   <- list(MBS = datobs[-BCX_more_than_one_index, grep("_NPPCR",names(datobs))[MSS_avail_order_index]],
                         MSS = datobs[-BCX_more_than_one_index, c(paste(pathogen_BrS_anyorder,"BCX",sep="_")[MSS_avail_order_index])],
                         MGS = NA)
          Y      <- datobs$Y[-BCX_more_than_one_index]
          X      <- datobs[-BCX_more_than_one_index,X_extra]
          Nd     <- Nd - length(BCX_more_than_one_index)
        } else {
          Mobs   <- list(MBS = datobs[,grep("_NPPCR",names(datobs))[MSS_avail_order_index]],
                         MSS = datobs[,c(paste(pathogen_BrS_anyorder,"BCX",sep="_")[MSS_avail_order_index])],
                         MGS = NA)
          Y      <- datobs$Y
          X      <- datobs[,X_extra]
        }
    }else{
      Mobs   <- list(MBS = datobs[,grep("_NPPCR",names(datobs))[MSS_avail_order_index]],
                     MSS = NA,
                     MGS = NA)
      Y      <- datobs$Y
      X      <- datobs[,X_extra]
    }
  }

  # order the whole pathogen list (excluing SS-only pathogens),
  # so that those have SS data comes first
  pathogen_BrS_ordered_by_MSS <- pathogen_BrS_anyorder[MSS_avail_order_index]

  # get the list of pathogen categories for each of the above ordered pathogens:
  path_cat_ind <- sapply(1:JBrS,function(i)
         which(pathogen_cat_lookup$X==pathogen_BrS_ordered_by_MSS[i]))

  if (!is.null(pathogen_SSonly)){
    pathogen_SSonly_cat <- pathogen_cat_lookup[sapply(1:JSSonly,function(i)
                         which(pathogen_cat_lookup$X==pathogen_SSonly[i])),]
    res <- list(Mobs = Mobs,
                   Y = Y,
                   X = X,
                   JSS  = JSS,
                   pathogen_BrS_ordered_by_MSS = pathogen_BrS_ordered_by_MSS,
                   pathogen_BrS_cat         = pathogen_cat_lookup[path_cat_ind,],
                   JSSonly = JSSonly,
                   pathogen_SSonly_cat = pathogen_SSonly_cat)
  } else{
    res <- list(Mobs = Mobs,
                 Y = Y,
                 X = X,
                 JSS  = JSS,
                 JSSonly = 0,
                 pathogen_BrS_ordered_by_MSS = pathogen_BrS_ordered_by_MSS,
                 pathogen_BrS_cat         = pathogen_cat_lookup[path_cat_ind,])
  }
  if (!is.null(extra_meas_nm)){
    res$extra_Mobs <- extra_Mobs
  }
  res
}
