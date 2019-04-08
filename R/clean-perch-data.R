#' Clean PERCH data
#'
#' \code{clean_perch_data} transforms a raw data table (row for subjects, column
#' for variables - named as \code{\{pathogen name\}_\{specimen\}\{test\}} for lab tests
#' or other covariates) into a list. It is designed for PERCH data format.
#'
#' @param clean_options The list of options for cleaning PERCH data.
#' Its elements are defined as follows:
#' 
#' \itemize{
#' \item{\code{raw_meas_dir}}{: The file path to the raw data;}
#' \item{\code{case_def}}{: Variable name in raw data for \strong{case} definition;}
#' \item{\code{case_def_val}}{: The value for \strong{case} definition;}
#' \item{\code{ctrl_def}}{: Variable name in raw data for \strong{control} definition;}
#' \item{\code{ctrl_def_val}}{: The value for \strong{control} definition;}
#' \item{\code{X_strat}}{: A vector of variable names for stratifying the data
#' to perform SEPARATE analyses;}
#' \item{\code{X_strat_val}}{: A list of values for \code{X_strat}. The output
#' data only have individuals with \code{identical(X_strat,X_strat_val)==TRUE}.
#' To perform analysis on a single site, say \code{"02GAM"}, use \code{X_strat="newSITE"} and
#' \code{X_strat_val=list("02GAM")};}
#' \item{\code{BrS_objects}}{: A list of BrS objects built by \code{\link{make_meas_object}};}
#' \item{\code{SS_objects}}{: A list of SS objects built by \code{\link{make_meas_object}};}
#' \item{\code{X_extra}}{: A vector of covariate names for regression
#' or visualization;}
#' \item{\code{patho_taxo_dir}}{: The file path to the pathogen category or taxonomy
#' information (.csv). The information should be as complete as possible for a particular
#' analysis. If not, the pathogen without taxonomy information could not be
#' assigned to bacterial or viral groups (see \code{\link{plot_group_etiology}});
#' See \code{\link{assign_taxo_cause_list}} that requires this taxonomy information.}.
#'}
#'
#'
#' @return A List: \code{list(Mobs,Y,X)}
#' \itemize{
#' \item \code{Mobs} A list of bronze- (\code{MBS}), silver- (\code{MSS}),
#' and gold-standard (\code{MGS}, if available) measurements. See the formats
#'  of these measurements in \code{\link{extract_data_raw}}.
#' \item \code{Y} 1 for case; 0 for control;
#' \item \code{X} Data frame of covariates for cases and controls. The covariate
#' names are specified in \code{X_extra};
#' }
#' This function does not re-order pathogens that only have silver-standard data.
#'
#' @seealso \link{make_meas_object} for wrapping information about a particular
#' type of measurement; \link{extract_data_raw} for reading raw data table
#' and organizing them into \code{data_nplcm} format. Also see \link{clean_combine_subsites}
#' for combining subsites and \link[lubridate]{parse_date_time} for parsing date.
#' 
#' @family data tidying functions
#' @export

clean_perch_data <- function(clean_options) {
  #
  # cleaning options:
  #
  raw_meas_dir <- clean_options$raw_meas_dir
  case_def     <- clean_options$case_def
  case_def_val <- clean_options$case_def_val
  ctrl_def     <- clean_options$ctrl_def
  ctrl_def_val <- clean_options$ctrl_def_val
  X_strat      <- clean_options$X_strat
  X_strat_val  <- clean_options$X_strat_val
  BrS_objects  <- clean_options$BrS_objects
  SS_objects   <- clean_options$SS_objects
  X_extra      <- clean_options$X_extra
  
  # clean dates: use specified date format; if not specified, try the formats
  # specified in the "else" sub-clause:
  if (!is.null(clean_options$date_formats)) {
    date_formats      <- clean_options$date_formats
  }else{
    date_formats      <- c("%d%B%Y","%d%B%y","%d%b%Y")
  }
  
  #
  # preparing data: 1. combine subsites; 2. remove "X_" from column names:
  #
  
  # combine two sub-sites: 06NTH and 07STH  --> THA; 08MBA and 09DBA  --> BAN:
  PERCH_data_with_newSITE <- clean_combine_subsites(
    raw_meas_dir,
    subsites_list = list(c("06NTH","07STH"),
                         c("08MBA","09DBA")),
    newsites_vec  = c("06THA","07BAN")
  )
  cleanName  <-
    delete_start_with("X_",names(PERCH_data_with_newSITE))
  colnames(PERCH_data_with_newSITE) <- cleanName
  
  # write cleaned data into working directory:
  #meas_dir <-
  #  file.path(dirname(raw_meas_dir),paste0("prepared_",basename(raw_meas_dir)))
  #write.csv(PERCH_data_with_newSITE, meas_dir, row.names = FALSE)
  #rm(PERCH_data_with_newSITE)
  
  #
  # import prepared data:
  #
  
  # stratify on variables:
  X_strat_nm_tmp_case  <- c(X_strat, case_def)
  X_strat_val_tmp_case <- c(X_strat_val, case_def_val)
  
  X_strat_nm_tmp_ctrl  <- c(X_strat, ctrl_def)
  X_strat_val_tmp_ctrl <- c(X_strat_val, ctrl_def_val)
  
  
  # get data for cases and controls
  datacase <-
    extract_data_raw(
      PERCH_data_with_newSITE,X_strat_nm_tmp_case,X_strat_val_tmp_case,
      c(BrS_objects,SS_objects),X_extra
    )
  datactrl <-
    extract_data_raw(
      PERCH_data_with_newSITE,X_strat_nm_tmp_ctrl,X_strat_val_tmp_ctrl,
      c(BrS_objects,SS_objects),X_extra
    )
  
  # clean some covariates:
  if ("ENRLDATE" %in% X_extra) {
    # if ENRLDATE is in the dataset, transform the date format
    #Rdate.case <- as.Date(datcase_X_extra_subset$ENRLDATE, "%d%B%Y")
    #Rdate.ctrl <- as.Date(datctrl_X_extra_subset$ENRLDATE, "%d%B%Y")
    Rdate.case <-
      lubridate::parse_date_time(datacase$X$ENRLDATE, date_formats)
    Rdate.ctrl <-
      lubridate::parse_date_time(datactrl$X$ENRLDATE, date_formats)
    
    uniq.month.case <-
      unique(paste(
        lubridate::month(Rdate.case),lubridate::year(Rdate.case),sep = "-"
      ))
    uniq.month.ctrl <-
      unique(paste(
        lubridate::month(Rdate.ctrl),lubridate::year(Rdate.ctrl),sep = "-"
      ))
    
    #symm.diff.dates <- as.set(uniq.month.case)%D%as.set(uniq.month.ctrl)
    #if (length(symm.diff.dates)!=0){
    #  cat("Cases and controls have different enrollment months:","\n")
    #  print(symm.diff.dates)
    #}
    datacase$X$ENRLDATE <- Rdate.case
    datactrl$X$ENRLDATE <- Rdate.ctrl
  }
  if ("patid" %in% X_extra) {
    datacase$X$patid <- as.character(datacase$X$patid)
    datactrl$X$patid <- as.character(datactrl$X$patid)
  }
  if ("newSITE" %in% X_extra) {
    datacase$X$newSITE <- as.character(datacase$X$newSITE)
    datactrl$X$newSITE <- as.character(datactrl$X$newSITE)
  }
  
  # combine case and control measurements:
  all_MBS <- list()
  for (i in seq_along(datacase$Mobs$MBS)) {
    all_MBS[[i]] <- rbind(datacase$Mobs$MBS[[i]],datactrl$Mobs$MBS[[i]])
  }
  names(all_MBS) <- names(datacase$Mobs$MBS)
  
  all_MSS <- list()
  for (i in seq_along(datacase$Mobs$MSS)) {
    all_MSS[[i]] <- rbind(datacase$Mobs$MSS[[i]],datactrl$Mobs$MSS[[i]])
  }
  names(all_MSS) <- names(datacase$Mobs$MSS)
  if (length(all_MSS) == 0) {
    all_MSS <- NULL
  }
  
  all_X  <- rbind(datacase$X,datactrl$X)
  all_Y  <- c(rep(1,nrow(datacase$X)),rep(0,nrow(datactrl$X)))
  
  data_all <- list(Mobs = list(MBS = all_MBS,MSS = all_MSS),
                   X   = all_X,
                   Y   = all_Y)
}

#' Import Raw PERCH Data
#'  
#' \code{extract_data_raw} imports and converts the raw data to analyzable format
#'
#' @param dat_prepared  The data set prepared in \code{clean_perch_data}.
#' @param strat_nm        The vector of covariate names to separately extract data.
#'    For example, in PERCH data cleaning, \code{X = c("newSITE","CASECONT")}.
#' @param strat_val     The list of covariate values to stratify data.
#'    Each element corresponds to elements in \code{X}. For example, in PERCH
#'    data cleaning, \code{Xval = list("02GAM","1")}.
#' @param meas_object A list of bronze-standard or silver-standard measurement
#' objects made by function \code{\link{make_meas_object}}.
#' @param extra_covariates The vector of covariate name for regression purposes.
#'   The default is NULL, which means not reading in any covariate.
#'
#' @return A list of data.
#' \itemize{
#' \item{Mobs}{
#'    \itemize{
#'    \item{MBS} A list of Bronze-Standard (BrS) measurements.
#'    The names of the list take the form of \code{specimen}_\code{test}. 
#'    Each element of the list is a data frame. The rows of the data frame 
#'    are for subjects; the columns are for measured pathogens.
#'    \item{MSS} A list of Silver-Standard (SS) measurements. 
#'    The formats are the same as \code{MBS} above.
#'    \item{MGS} A list of Gold-Standard (GS) measurements. 
#'    It equals \code{NULL} if no GS data exist.
#'    }
#' }
#' \item{X} A data frame with columns specified by \code{extra_covariates}.
#' }
#' 
#' @family raw data importing functions
#' 
#' @seealso \code{\link{clean_perch_data}}
#' 
#' @export

extract_data_raw <- function(dat_prepared,strat_nm,strat_val,
                             meas_object,
                             extra_covariates = NULL) {
  #dat_prepared <-
  #  utils::read.csv(meas_dir,header = TRUE,stringsAsFactors = FALSE)
  cleanName      <- colnames(dat_prepared)
  ind_this_strat  <- 1:nrow(dat_prepared)
  for (j in 1:length(strat_nm)) {
    ind_this_strat = ind_this_strat[which(dat_prepared[ind_this_strat,strat_nm[j]] ==
                                            strat_val[[j]] &
                                            !is.na(dat_prepared[ind_this_strat,strat_nm[j]]))]
  }
  dat <- dat_prepared[ind_this_strat,]
  
  #
  # get lab measurements:
  #
  Mobs <- list(MBS = list(),MSS = list(), MGS = list())
  count <- c(MBS = 0,MSS = 0,MGS = 0)
  
  for (m in seq_along(meas_object)) {
    object_tmp <- meas_object[[m]]
    meas_tmp  <- read_meas_object(object_tmp,dat)
    # increment number of that quality:
    count[meas_tmp$position] <- count[meas_tmp$position] + 1
    # add to the Mobs:
    curr_pos <- count[meas_tmp$position]
    Mobs[[meas_tmp$position]][[curr_pos]] <- meas_tmp$meas
    names(Mobs[[meas_tmp$position]])[curr_pos] <-
      object_tmp$nm_spec_test
  }
  
  Mobs[which(count == 0)] <- NULL
  
  
  # get extra covariates used for regression purposes:
  X <- list()
  if (!is.null(extra_covariates)) {
    for (i in seq_along(extra_covariates)) {
      if (!extra_covariates[i] %in% colnames(dat)) {
        stop("==[baker]",extra_covariates[i]," is not in the data set. Delete this covariate.==","\n")
      } else {
        X[[i]] = dat[,extra_covariates[i]]
      }
    }
  }
  
  names(X) = extra_covariates
  X <- as.data.frame(X)
  X_stratified <- dat[,strat_nm]
  colnames(X_stratified) <- strat_nm
  cbind(X_stratified,X)
  
  make_list(Mobs,X)
}

#' Combine subsites in raw PERCH data set
#'
#' In the Actual PERCH data set, a study site may have multiple subsites.
#' \code{clean_combine_subsites} combines all the study subjects from the same site.
#'
#' @param raw_meas_dir The file path to the raw data file (.csv)
#' @param subsites_list The list of subsite group names. Each group is a vector of
#'   subsites to be combined
#' @param newsites_vec A vector of new site names. It has the same length as
#' \code{"subsites_list"}
#' @return A data frame with combined sites
#'
#' @export
clean_combine_subsites <-
  function(raw_meas_dir,subsites_list,newsites_vec) {
    if (length(subsites_list) != length(newsites_vec)) {
      stop(
        "==[baker] The length of new site names is not equal to the number of subsite groups!
        Make them equal.=="
      )
    } else{
      tmp.dat         <- utils::read.csv(raw_meas_dir)
      tmp.dat$newSITE <- as.character(tmp.dat$SITE)
      for (i in 1:length(newsites_vec)) {
        tmp.dat$newSITE[which(tmp.dat$SITE %in% subsites_list[[i]])] = newsites_vec[i]
      }
      return(tmp.dat)
    }
  }


#' Make measurement slice
#' 
#' Wrap the information about a particular type of measurement, e.g., NPPCR.
#'
#'@param patho A vector of pathogen names
#'@param specimen Specimen name
#'@param test Test name
#'@param quality Quality category: any of "BrS", "SS" or "GS".
#'@param cause_list The vector of potential latent status
#'@param sep_char a character string that separate the pathogen names and the 
#'specimen-test pair; Default to \code{"_"}
#'
#'@return A list with measurement information
#'\itemize{
#'\item{\code{quality}} same as argument
#'\item{\code{patho}} same as argument
#'\item{\code{name_in_data}} the names used in the raw data to locate these measurements
#'\item{\code{template}} a mapping from \code{patho} to \code{cause_list}.
#'\code{NROW = length(cause_list)+1};
#'\code{NCOL = length(patho)}. This value is crucial in model fitting to determine
#'which measurements are informative of a particular category of latent status.
#'\item{\code{specimen}} same as argument
#'\item{\code{test}} same as argument
#'\item{\code{nm_spec_test}} paste \code{specimen} and \code{test} together
#'}
#'
#'@family data standardization functions
#'
#'@seealso \code{\link{make_template}}
#'@export

make_meas_object <-
  function(patho,specimen,test,quality,cause_list,sep_char="_") {
    nm_spec_test <- paste0(specimen,test)
    name_in_data <- paste(patho,nm_spec_test,sep = sep_char)
    template <- make_template(patho,cause_list)
    make_list(quality,patho,name_in_data,template,specimen,test,nm_spec_test)
  }

#' Read measurement slices
#' 
#' @param object Outputs from \code{\link{make_meas_object}} 
#' @param data Raw data with column names 
#' \code{\{pathogen name\}_\{specimen\}\{test\}}
#' 
#' @return A list with two elements: \code{meas}-data frame with measurements;
#' \code{position}-see \code{\link{lookup_quality}}
#' 
#' @family raw data importing functions
#' @export
#' 
read_meas_object <- function(object,data) {
  position <- lookup_quality(object$quality)
  exist <- sapply(object$name_in_data,grep,colnames(data))
  not_exist_index <- which(unlist(lapply(exist,length))==0)
  if (length(not_exist_index)>0){
    stop(paste0("==[baker]",paste(object$name_in_data[not_exist_index],collapse=", ")," is not in data!=="))  
  }
  meas <- data[,object$name_in_data,drop = FALSE]
  colnames(meas) <- object$patho
  make_list(meas,position)
}


#' assign taxonomy information to cause list
#' 
#' @param cause_list See \code{cause_list} in \code{model_options$likelihood} (see 
#' parameter settings in \code{\link{nplcm}})
#' @param dir_taxo File path to the .csv file storing two columns of information:
#' \code{pathogen}, and \code{pathogen_type}.
#' @return A vector of factors of taxonomy information (currently for "B","F","V"). 
#' The names of the vector are pathogen names.
#' 
#' @export
assign_taxo_cause_list <- function(cause_list, dir_taxo=NULL){
  patho_taxo <- utils::read.csv(dir_taxo,header=TRUE,stringsAsFactors = FALSE)
  patho_taxo$pathogen_type <- factor(patho_taxo$pathogen_type,levels=c("B","F","V"))
  res_taxo   <- patho_taxo$pathogen_type[match(cause_list,patho_taxo[,1])]
  names(res_taxo) <- cause_list
  res_taxo
}



