#' Color palette for magenta-to-cyan colors (mimicing "cool" palette in matlab)
#'
#' @param x
#'
#' @keywords internal
#' @returns color pallette for ggplot
#' @export
#'
#' @examples cool_matlab()
cool_matlab <- function(x = 256) {
  seq_vals <- seq(0, 1, length.out = x)
  colors <- rgb(seq_vals, 1 - seq_vals, 1)
  return(colors)
}
