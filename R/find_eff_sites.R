#' Find effective number of sites needed to reach target coverage
#'
#' @param coverage_vec a vector of coverage values
#' @param target coverage value to achieve with n sites
#'
#' @returns number of sites necessary to achieve a target coverage value
#' @export
#' @keywords internal
#'
find_eff_sites <- function(coverage_vec, target = c(1)){

  # if coverage_vals is not supplied,
  # function returns number of sites needed for effectively 100% coverage


  # if coverage_vals is supplied,
  # function returns an array of number of sites needed for specified coverage values
  vec <- lapply(target,
                function(x){dplyr::first(which(coverage_vec >= x))}) %>% unlist()

  if(any(is.na(vec))){
    print('Warning: coverage does not exceed (some) target. Output vector includes NAs')
  }

  return(vec)

}
