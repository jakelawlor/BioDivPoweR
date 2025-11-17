#' Extract indices to keep in each bin
#'
#' Find index of bootstraps to extract within bins with count >= min_exp_n
#'
#' @param .p3 histogram of richness differences between bootstrapped communities
#' @param .summary summary df of richness differences
#' @param .min_exp_n threshold for keeping histogram bins
#' @return list of numeric IDs indicating indices of bootstraps to keep in each qualifying effect size bin
#' @noRd
.extract_ids_to_keep <- function(.p3, .summary, .min_exp_n){

  p3_breaks <- .p3 %>% ggplot2::layer_data() %>%
    dplyr::select(xmin,xmax,count,x) %>%
    dplyr::mutate(x = round(x,3)) %>%
    dplyr::mutate(iter = 1:nrow(.)) %>%
    dplyr::filter(count >= .min_exp_n,
                  x != 0) %>%
    split(f = .$iter)

  # make a list of min_exp_n pair IDS to save inside each eff_size_bin.
  ids_in_bins <-
    # find ids within each bin
    purrr::map(
    .x = p3_breaks,
    .f = ~which(.summary$log2_rich_diff > .x$xmin &
                  .summary$log2_rich_diff <= .x$xmax &
                  .summary$log2_rich_diff != 0)[1:.min_exp_n] # keep 1:min_exp_n
  ) %>%
    # set names of bins to the center x value of each bin
    purrr::set_names(purrr::map(p3_breaks,~unique(.x$x)))
}



#' Extract bootstraps to keep
#'
#' subset rarefied bootstrapped communities to only indexes within qualifying bins
#'
#' @param .ids_to_keep list of indices to keep within qualifying bins
#' @param .rarefied_boots full rarefied bootstraps list to subset
#' @return heirarchical list of bootstrapped communities to keep. List of all qualifying bins. Each element is a sublist of two community types. Each sub-sublist is min_exp_n matrixes of bootstrapped rarefied communities
#' @noRd
.keep_selected_boots <- function(.ids_to_keep, .rarefied_boots, .eff_coverage){

  # select the correct rarified bootstraps
  rarefied_boots_keep <- purrr::map(
    .x = .ids_to_keep,
    .f = ~list(comm1 = .rarefied_boots[[1]][.x],
               comm2 = .rarefied_boots[[2]][.x]) %>%
      purrr::set_names(names(.rarefied_boots))
  )

  # cut all coverage vectors (originally 1:full_n_sites) to 1:rarefied_n_sites
  rarefied_coverages <- purrr::map2(
    .x = .eff_coverage,
    .y = .rarefied_boots,
    .f = ~purrr::map2(
      .x = .x,
      .y = .y,
      .f = ~.x[1:nrow(.y)]
    )
  )

  # select the correct coverage vectors,
  # but note these are not rarefied so will need to be cut down to the new sample sizes
  rarefied_coverages_keep <- purrr::map(
    .x = .ids_to_keep,
    .f = ~list(coverage_1 = rarefied_coverages[[1]][.x],
               coverage_2 = rarefied_coverages[[2]][.x]) %>%
      purrr::set_names(names(.eff_coverage))
  )


  boots_tibble <- purrr::map(
    .x = rarefied_boots_keep,
    .f = ~tibble::as_tibble(.x)
  ) %>%
    dplyr::bind_rows(.id = "eff_size") %>%
    purrr::set_names(c(names(.)[1],paste0(names(.)[2:3],".matrix")))

  coverage_tibble <- purrr::map(
    .x = rarefied_coverages_keep,
    .f = ~tibble::as_tibble(.x)
  ) %>%
    dplyr::bind_rows(.id = "eff_size") %>%
    purrr::set_names(c(names(.)[1],paste0(names(.)[2:3],".coverage")))


  return(
    dplyr::bind_cols(boots_tibble,
                     coverage_tibble %>%
                       dplyr::select(-eff_size))
    )

}

