#' Create bootstrap arrays
#'
#' Create the initial bootstrap arrays that will later be rarefied
#'
#' @param data pilot dataset - either single treatment or two treatment
#' @param method specification for "single" or "two" treatment workflow
#' @param category_col name of the column that lists the treatment variable. Imported from wrapper function
#' @param n_boots number of community pairs to simulate. Imported from wrapper funciton
#' @return list of length 2, each element is an array with dimensions (max_sites, max_spp, n_boots)
#' @noRd
.create_boot_arrays <- function(.pilot, .method, .category_col = NULL, .n_boots = n_boots){

  # get array dimensions ----------------------------------------------------
  # first, the number of sites (array rows),
  # where n = total_sites for single treatment,
  # and total_sites within each treatment for two-treatment
  n_community_sites <- switch(
    .method,
    "single" = list("community 1" = nrow(.pilot),
                    "community 2" = nrow(.pilot)),
    "two" =  purrr::map(.x = sort(unique(.pilot[[.category_col]])) %>% purrr::set_names(.),
                        .f = ~sum(.pilot[[.category_col]] == .x))
  )


  # then, maximum number of species (array columns)
  max_n_spp <-  switch(
    .method,
    "single" = list(
      "comm1" = sum(sapply(.pilot, is.numeric) | sapply(.pilot, is.logical)),
      "comm2" = sum(sapply(.pilot, is.numeric) | sapply(.pilot, is.logical))),
    "two" = purrr::map(
      .x =  split(.pilot, .pilot[[.category_col]]),
      .f = ~sum(colSums(.x[,-c(1:2)]) > 0)
    )
  )


  # make empty arrays -------------------------------------------------------
  boot_arrays <- purrr::map2(
    .x = c(n_community_sites),
    .y = max_n_spp,
    .f = ~array(
      NA, # fill with zeros
      # set dimensions
      dim = c(
        .x, # number of sites in original community
        .y, # number of species in original community
        #(this will be the same for the 2 groups)
        .n_boots # number of iterations (default 5000)
      )
    )
  )


  # isolate source communities -------------------------------------------------------------
  # isolate communities in each treatment.
  # this will create a comm_1 and comm_2 matrix, which will be a species x site
  # matrices of only the "treatment 1" or "treatment 2" communities, or if
  # method = "single", they will be identical to the pilot matrix
  comm_1 <- switch(.method,
                   "single" = as.matrix(.pilot[,sapply(.pilot, is.numeric)]),
                   "two" = as.matrix(.pilot[.pilot[[.category_col]] == names(n_community_sites)[1], -c(1, 2)]))
  comm_1 <- comm_1[,colSums(comm_1) > 0] # remove columns where species are not present

  comm_2 <- switch(.method,
                   "single" = as.matrix(.pilot[,sapply(.pilot, is.numeric)]),
                   "two" = as.matrix(.pilot[.pilot[[.category_col]] == names(n_community_sites)[2], -c(1, 2)]))
  comm_2 <- comm_2[,colSums(comm_2) > 0] # remove columns where species are not present


  # fill boot arrays --------------------------------------------------------
  # randomly sample pilot communities (comm_1 and comm_2) n_boots times, with replacement
  # for iterations 1-n_boots
  for(i in 1:.n_boots){

    # fill the "page" of the community 1 n_community_sites[[1]] randomly selected sites
    # from the original category 1 sites
    boot_arrays[[1]][,,i] <- comm_1[sample(nrow(comm_1),n_community_sites[[1]], replace = T),]

    # fill in each "page" in community 2 with n_community_sites[[2]] randomly sampled sites
    # from the original category 2 sites.
    boot_arrays[[2]][,,i] <- comm_2[sample(nrow(comm_2),n_community_sites[[2]], replace = T),]

  }
  rm(comm_1, comm_2)


  # convert from arrays to lists, so each bootstrap is its own list element
  # we're doing this because not all species will be inside each bootstrapped
  # community, and we'll need to remove columns for zero-occurrence species
  boot_arrays <- purrr::map(
    .x = boot_arrays,
    .f = ~lapply(seq(dim(.x)[3]), function(page) .x[,,page])
  )

  # filter to only species present at least once
  # so now each bootstrap will have a different number of columns.
  boot_arrays <- purrr::map(
    .x = boot_arrays,
    .f = ~purrr::map(
      .x = .x,
      .f = ~.x[,colSums(.x) > 0]
    )
  )

  return(boot_arrays)
}



#' find coverage
#'
#' @param .full_boots list of n_boots community pairs from `.create_boot_arrays`
#' @param method specification for "single" or "two" treatment workflow
#' @return list of length 2, each element is an array with dimensions (max_sites, n_species - custom per boot, n_boots)
#' @noRd
.find_rarefied_coverage <- function(.full_boots){

  # a. find average occupancies of all species in each boot
  avg_occupancies <- purrr::map(
    .x = .full_boots,
    .f = ~purrr::map(
      .x = .x,
      .f = ~colMeans(.x)
    )
  )

  # b. find achieved coverage in 1:n_true_samples samples in simulated communities.
  eff_coverage <- purrr::map2(
    .x = avg_occupancies, # occupancy vectors
    .y = .full_boots,
    .f = ~purrr::map(1:length(avg_occupancies[[1]]),
                     .f = function(y) find_coverage(
                       occupancy = .x[[y]],
                       k = nrow(.y[[y]])
                     ))
  )
  # returns a list of two (comm1 and comm2), each n_boots items long,
  # each containing the coverage at 1:max_sites samples

  return(eff_coverage)

}


#' Rarefy bootstrapped communities
#'
#' @param .full_boots list of n_boots community pairs from `.create_boot_arrays`
#' @param .eff_coverage list of coverage vectors from 1:max_sites for each bootstrapped community
#' @return list of length 2, each element is an array with dimensions (max_sites, n_species - custom per boot, n_boots)
#' @noRd
.rarefy_boots <- function(.full_boots, .eff_coverage) {

  n_sites_equal_coverage <- purrr::map2(
    .x = .eff_coverage[[1]], # list of coverage values for community 1 at 1-max_sites
    .y = .eff_coverage[[2]], # list of coverage values for community 2 at 1-max_sites
    .f = ~c(
      # find first number of comm1 samples where coverage is
      # greater than the lower max of the two coverage values
      "comm1" = dplyr::first(which(.x >= min(max(.x), max(.y)))),
      # find first number of comm2 samples where coverage is
      # greater than the lower max of the two coverage values
      "comm2" = dplyr::first(which(.y >= min(max(.x), max(.y))))
    ) %>%
      purrr::set_names(names(.full_boots))
  ) %>%
    purrr::transpose()

  # e. keep only number of samples in each list to match the coverages
  # to the closest value. So we'll sample the full bootstrap (n_sites rows)
  # without replacement, n_sites_equal_coverage times.
  boots_rare <- purrr::map2(
    .x = .full_boots,
    .y = n_sites_equal_coverage,
    .f = ~ purrr::map2(.x, .y,
                       ~ .x[
                         sample(nrow(.x), .y, replace = F),
                         ,
                         drop = FALSE])
    # note that drop=F prevents matrices from becoming vectors if there is only one row
  )

  # f. filter to only present species (so remove columns with a sum of 0)
  boots_rare <- purrr::map(
    .x = boots_rare,
    .f = ~purrr::map(
      .x = .x,
      .f = ~.x[,colSums(.x) > 0, drop = FALSE]
    )
  )

  return(boots_rare)

}


#' Summarize richness change in bootstrapped communities
#'
#' @param .rarefied_boots list of equal-coverage bootstrapped arrays from `.rarefy_boots`
#' @return a tibble of n_boots rows, listing the species richness in rarefied community pairs in all bootstraps, and the log2 ratio richness difference between them
#' @noRd
.summarize_boot_diff <- function(.rarefied_boots) {

  summary <-   # find richness change
    purrr::map2(
      .x = .rarefied_boots[[1]], # community 1 bootstraps
      .y = .rarefied_boots[[2]], # community 2 bootstraps
      .f = ~c("comm_1_rich" = ncol(.x),
              "comm_2_rich" = ncol(.y))
    ) %>%
    dplyr::bind_rows(.id = "boot") %>%
    dplyr::mutate(log_rich_change = log2(comm_2_rich / comm_1_rich)) %>%
    purrr::set_names(c("boot",paste0(names(.rarefied_boots),".rich"),"log2_rich_diff"))

  return(summary)

}


