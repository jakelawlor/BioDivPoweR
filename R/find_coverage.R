#' Find coverage of ecological community surveys at 1-k samples
#'
#'
#'
#' @param occupancy vector of average occupancy of all species in a community (>0)
#' @param k max number of samples to test. Coverage should asymptote at a large number of samples.
#' @keywords internal
#' @returns vector of coverage achieved in k samples of an ecological community
#' @export
find_coverage <- function(occupancy, k = 1000){

  sapply(1:k,
         function(x){mean(1 - (1-occupancy)^x)})


}
