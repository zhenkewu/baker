#' Combine subsites in raw PERCH data set
#'
#' In the Actual PERCH data set, a study site may have multiple subsites. 
#' \code{clean_combine_subsites} combines all the study subjects from the same site.
#'
#' @param RawMeasDir The file path to the raw data file (.csv)
#' @param subsites_list The list of subsite group names. Each group is a vector of
#'   subsites to be combined
#' @param newsites_vec A vector of new site names. It has the same length as
#' \code{"subsites_list"}
#' @return A data frame with combined sites
#'
#' @export

clean_combine_subsites <- function(RawMeasDir,subsites_list,newsites_vec){
    if (length(subsites_list)!=length(newsites_vec)){
      stop("The length of new site names is not equal to the number of subsite groups!
           Make them equal.")
    } else{
      tmp.dat         <- read.csv(RawMeasDir)
      tmp.dat$newSITE <- as.character(tmp.dat$SITE)
      for (i in 1:length(newsites_vec)){
        tmp.dat$newSITE[which(tmp.dat$SITE %in% subsites_list[[i]])] = newsites_vec[i]
      }
      return(tmp.dat)
    }
  }
