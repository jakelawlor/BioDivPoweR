#' Find target detection rate
#'
#' This function finds the target detection rate (coverage) that indicates that we have sampled "all" of the species in the community. Because coverage is asymptotic at 100, here, we find the point at which coverage indicates we have fewer than half of one species left to find.
#'
#' @param matrix
#'
#' @returns value of coverage at which we've found sufficient species in the community.
#' @export
#'
#' @examples if(FALSE){find_target_detection_rate(pilot_single_trt)}
#' @noRd
find_target_detection_rate <- function(matrix){

  # first, identify numeric columns in the pilot dataset
  # (excluding site and/or category columns)
  numcols <- unlist(lapply(matrix, is.numeric), use.names = FALSE)

  # next, find richness, or the total number of species (columns) that are
  # present in at least one site (colSum >0)
  richness_raw <- sum(colSums(matrix[,numcols]) > 0)

  # find target detection rate to detect all species
  # this is the point where coverage is <.5 species away from detecting the whole community
  target_detection_rate <- (1 / (1 + .5 / richness_raw))

  return(target_detection_rate)
}

