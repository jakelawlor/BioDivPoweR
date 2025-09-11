#' Bootstrap Biodiversity Samples from Two-Treatment Pilot
#'
#' `bootstrap_two_trt` creates `n_boots` simulated communities from a pilot biodiversity survey, pulled from randomly drawn sites from the `pilot` survey, with replacement. Bootstrapped pairs are compared for differences in richness, and differences in richness (log2 ratio) are sorted into a histogram with `n_eff_size_bins` bins. When possible, `min_exp_n` experiments (pairs of simulated communities, rarefied to equal sample coverage) are retained and kept as the input for the following function.
#'
#' @param pilot species-by-site matrix from a pilot biodiversity survey with two treatments, in which the first column lists site names or codes, the second lists the site category or treatment, and all following columns are named with species names and contain binary (0-1) values. See `data("pilot_two_trts")` for an example.
#' @param category_col the column name in `pilot` that specifies the treatment of sampling sites.
#' @param n_boots The number of simulated community pairs to bootstrap, which will be compared to each other for differences in species richness. Increasing the number of bootstrapped communities should result in more effect sizes (differences in richness between simulated pairs) qualifying for the power analysis (a histogram of richness differences from `n_boots` communities will be sorted into `n_eff_size_bins`, and bins that surpass the `min_exp_n` threshold will be retained for following steps).
#' @param n_eff_size_bins n_eff_size_bins The total number of bins to separate the histogram of  richness differences detected between simulated communities. Increasing `n_eff_size_bins` can offer higher resolution in the following steps (more bins to qualify), but can result in "toothy" histograms in which sequential bins are not all filled. Decreasing `n_eff_size_bins` should fix "toothiness".
#' @param min_exp_n The threshold number of simulated pairs to qualify a given effect size bin for retention. Higher `min_exp_n` values will offer greater resolution in following steps, but will decrease the number of effect sizes that qualify.
#' @param seed Random seed. Defaults to 1 so runs of the same data will provide the same answers, but since the simulations all rely on random draws, changing the seed will result in different answers.
#'
#'
#' @returns a tibble of `min_exp_n` simulated communities, rarified to equal coverage, within all effect size bins that qualify. Additionally, prints a histogram of simulated community richness values, highlighting the bins that are kept for the next step.
#' @export
#'
#' @examples if(FALSE){boostrap_two_trt(pilot2, category_col = "veg_type")}
bootstrap_two_trt <- function(pilot,
                              category_col = NULL,
                              n_boots = 5000,
                              n_eff_size_bins = 40,
                              min_exp_n = 40,
                              seed = 1){

  stopifnot("please specify the name of the treatment column in pilot dataset" = !is.null(category_col),
            "Pilot data need two unique treatments" = length(unique(pilot[[category_col]])) == 2 )

  # set random seed
  set.seed(seed)

  # create a vector of the two types of communities
  community_types <- unique(pilot[[category_col]]) %>% purrr::set_names(.)


  # make empty array to store results
  # cols = n_species in pilot dataset
  # rows = n_sites in pilot dataset within that site
  # dim3 = n communities to compare (2 for richness comparisons)
  # dim4 = n trials (set by user, default to 5000)

  # find number of sites in each group in the original pilot data
  n_community_sites <- purrr::map(.x = community_types,
                                  .f = ~sum(pilot[[category_col]] == .x))
  # find number of sites in each group in the original pilot data
  n_community_sites <- purrr::map(.x = community_types,
                                  .f = ~sum(pilot[[category_col]] == .x))

  # find maximum number of species
  max_n_spp <-   sum(sapply(pilot, is.numeric) | sapply(pilot, is.logical))


  # make empty boot array ---------------------------------------------------
  # make a list of two arrays, one for community 1, one for community 2
  # note that to make it easy, here, the max number of species will be the same
  boot_arrays <- purrr::map(
    .x = c(n_community_sites),
    .f = ~array(
      NA, # fill with zeros
      # set dimensions
      dim = c(
        .x, # number of sites in original community
        max_n_spp, # number of species in original community
        #(this will be the same for the 2 groups)
        n_boots # number of iterations (default 5000)
      )
    )
  )



  # fill in the arrays ----------------------------------------------
  cat("Bootstrapping empirical communities", n_boots,"times... \n")

  # extract all the sites in each community. This is what we'll sample
  # rows from to make the simulated bootstrap communities
  comm_1 <- as.matrix(pilot[pilot[[category_col]] == community_types[1], -c(1, 2)])
  comm_2 <- as.matrix(pilot[pilot[[category_col]] == community_types[2], -c(1, 2)])

  # NOTE that there are different numbers of species in each community,
  # but the vectors are both sized as if they both have max_n_spp.
  sum(colSums(comm_1)>0)
  sum(colSums(comm_2)>0)


  # for iterations 1-5000 (or n)
  for(i in 1:n_boots){

    # fill the "page" of the community 1 n_group_sites[[1]] randomly selected sites
    # from the original category 1 sites
    boot_arrays[[1]][,,i] <- comm_1[sample(nrow(comm_1),n_community_sites[[1]], replace = T),]

    # fill in each "page" in community 2 with n_group_sites[[2]] randomly sampled sites
    # from the original category 2 sites.
    boot_arrays[[2]][,,i] <- comm_2[sample(nrow(comm_2),n_community_sites[[2]], replace = T),]

  }

  # rarefy pairs to equal coverage -----------------------------------------
  # simulated communities can only be compared to each other if they are
  # rarefied to equal sampling coverage, so for each pair, we'll reduce the
  # one with higher coverage to match the coverage of the lower one.
  cat("Rarefying", n_boots, "bootstraps to equal coverage... \n")

  # a. first find average occupancies of species in each community
  avg_occupancies <- purrr::map(
    .x = boot_arrays,
    .f = ~apply(.x,
                3, # apply across the third dimensions (number of bootstraps)
                function(page){colMeans(page)[colSums(page) > 0]}), # keep the means for species whose sum is above zero (e.g., they're present at least once)
    .progress = T
  )

  # b. find achieved coverage in 1:n_true_samples samples in simulated communities.
  eff_coverage <- purrr::map2(
    .x = avg_occupancies,
    .y = n_community_sites,
    .f = ~purrr::map(.x,
                     .f = function(y) BioDivSampler:::find_coverage(y, k = .y)),
    .progress = T
  )
  # returns a list of two (comm1 and comm2), each n_boots items long, each containing the coverage at 1:max_sites samples

  # c. find number of sites in each community necessary to equalize achieved coverage
  # find the first number of sites in 1:k to sample in each community
  n_sites_equal_coverage <- purrr::map2(
    .x = eff_coverage[[1]],
    .y = eff_coverage[[2]],
    .f = ~c(
      # find first number of comm1 samples where coverage is
      # greater than the lower max of the two coverage values
      "comm1" = dplyr::first(which(.x >= min(max(.x), max(.y)))),
      # find first number of comm2 samples where coverage is
      # greater than the lower max of the two coverage values
      "comm2" = dplyr::first(which(.y >= min(max(.x), max(.y))))
    ) %>%
      purrr::set_names(community_types)
  )


  # d. rarefy communities to number of samples that gives equal coverage
  # first, change original arrays into lists
  rarified_boots <- purrr::map(
    .x = boot_arrays,
    .f = ~lapply(seq(dim(.x)[3]), function(page) .x[,,page])
  )

  # e. keep only number of samples in each list to match the coverages
  # here, we're taking just rows 1:n_sites. Since sites are already randomized,
  # we don't need to pull random samples.
  rarified_boots[[1]] <- purrr::map2(
    .x = rarified_boots[[1]],
    .y = n_sites_equal_coverage,
    .f = ~.x[seq(1:.y[community_types[[1]]]),]
  )

  rarified_boots[[2]] <- purrr::map2(
    .x = rarified_boots[[2]],
    .y = n_sites_equal_coverage,
    .f = ~.x[seq(1:.y[community_types[[2]]]),]
  )

  # find richness in rarified boots
  summary <- purrr::map2(
    .x = rarified_boots[[1]],
    .y = rarified_boots[[2]],
    .f = ~c(
      "comm_1_rich" = sum( colSums(.x) > 0),
      "comm_2_rich" = sum( colSums(.y) > 0))) |>
    dplyr::bind_rows(.id = "boot") |>
    dplyr::mutate(log_rich_change = log2(comm_1_rich/ comm_2_rich)) %>%
    purrr::set_names(c("boot",paste0(community_types,".rich"),"log2_rich_change"))


  # make histograms ---------------------------------------------------------
  cat("Sorting differences to",n_eff_size_bins,"effect size bins... \n")

  p1 <- summary |>
    tidyr::pivot_longer(cols = contains(".rich"),
                        names_to = "community",
                        values_to = "richness",
                        names_pattern = "(.*).rich") %>%

   ggplot2::ggplot() +
    ggplot2::geom_histogram(
      ggplot2::aes(x = richness,
                   fill = community,
                   color = community),
      alpha = .5,
      center = 0,
      linewidth = .1,
      binwidth = ceiling(
        (max(summary[paste0(community_types,".rich")]) - min(summary[paste0(community_types,".rich")]))
        /20),
      position = "identity"
    ) +
    ggplot2::scale_fill_manual(values = c("grey35","grey80")) +
    ggplot2::scale_color_manual(values = c("grey20","grey20")) +
    ggplot2::theme(legend.position = "inside",
                   legend.position.inside = c(.02,.98),
                   legend.justification = c(0,1)) +
    ggplot2::labs(x = "Simulated Community Richness",
                  y = "Count",
                  fill = NULL,
                  color = NULL)

  print(p1)

  p1.2 <- summary |>
    dplyr::mutate(raw_change =
                    get(paste0(community_types[[1]],".rich")) -
                    get(paste0(community_types[[2]],".rich"))
                    ) %>%
    ggplot2::ggplot() +
    ggplot2::geom_histogram(ggplot2::aes(x = raw_change),
                            alpha = .5,
                            color = "grey20",
                            fill = "grey80",
                            center = 0,
                            linewidth = .1,
                            binwidth = 1,
                            position = "identity") +
    ggplot2::labs(x = "Raw richness difference",
                  y = "Count")

  print(p1.2)


  p2 <- ggplot2::ggplot(data = summary,
                        ggplot2::aes(x = log2_rich_change)) +
    ggplot2::geom_histogram(
      color = "grey20",
      fill = "grey90",
      linewidth = .1,
      bins = n_eff_size_bins,
      show.legend = F,
      center = 0
    ) +
    ggplot2::geom_hline(yintercept = min_exp_n,
                        linetype = "dashed",
                        linewidth = .4) +
    ggplot2::labs(x = "Log2 Fold-Difference in Richness",
                  y = "Count")

  # add a green box around the sections of the histograms to keep
  p2 <- p2 +
    ggplot2::geom_rect(data = p2 |> ggplot2::layer_data() |>
                         dplyr::filter(count > min_exp_n) |>
                         dplyr::mutate(eff_size_num = (round(x,3)))
                       ,
                       inherit.aes = F,
                       ggplot2::aes(xmin = xmin,
                                    xmax = xmax,
                                    ymin = 0,
                                    ymax = min_exp_n,
                                    fill = abs(eff_size_num)),
                       color = "grey20",
                       linewidth = .1, show.legend = F) +
    ggplot2::scale_fill_gradientn(colors = BioDivSampler:::cool_matlab())

  print(p2)
  # NOTE: Am i still supposed to absolute value here?

  # output the richness simulation dataset

  p2_breaks <- p2 %>% ggplot2::layer_data() %>%
    dplyr::select(xmin,xmax,count,x) %>%
    dplyr::mutate(x = round(x,3)) %>%
    dplyr::mutate(iter = 1:nrow(.)) %>%
    dplyr::filter(count > min_exp_n) %>%
    split(f = .$iter)

  # make a list of min_exp_n pair IDS to save inside each eff_size_bin.
  # note to save the RARIFIED bootstraps, not the FULL bootstrap
  ids_in_bins <- purrr::map(
    .x = p2_breaks,
    .f = ~which(summary$log2_rich_change > .x$xmin &
                  summary$log2_rich_change <= .x$xmax)[1:min_exp_n]
  ) %>%
    purrr::set_names(purrr::map(p2_breaks,~unique(.x$x)))


  # save selected experiments ------------------------------------------
  cat("Saving", min_exp_n, "boostraps from each qualifying bin... \n")

  # select the correct rarified bootstraps
  rarefied_boots_keep <- purrr::map(
    .x = ids_in_bins,
    .f = ~list(comm1 = rarified_boots[[1]][.x],
               comm2 = rarified_boots[[2]][.x])
  )

  # now subset the effective coverage (originally from
  # 1-to-total_site_num) sites to now be from 1 to effective site num
  # (the number at which the two communities have equal coverage).
  eff_coverage_rare <- purrr::map2(
    .x = eff_coverage,              # list of comm1, comm2
    .y = community_types,    # their names
    .f = function(comm_list, comm_name) {
      purrr::map2(
        .x = comm_list,             # each of the 5000 coverage vectors
        .y = n_sites_equal_coverage,
        .f = ~ .x[seq(.y[comm_name])]
      )
    }
  )

  # keep the coverage curves for only the selected sites, at only the lengths of the rarefied communities
  coverage_curves_keep <- purrr::map(
    .x = ids_in_bins,
    .f = ~list(comm1 = eff_coverage_rare[[1]][.x],
               comm2 = eff_coverage_rare[[2]][.x]) %>%
      purrr::set_names(community_types)

  )


  # convert to nested tibble to return
  boots_tibble <- purrr::map(
    .x = rarefied_boots_keep,
    .f = ~tibble::as_tibble(.x)
  ) %>%
    dplyr::bind_rows(.id = "eff_size")
  names(boots_tibble) <- c("eff_size",paste0(community_types,"_matrix"))

  coverage_tibble <- purrr::map(
    .x = coverage_curves_keep,
    .f = ~tibble::as_tibble(.x)
  ) %>% dplyr::bind_rows()
  names(coverage_tibble) <- c(paste0(community_types,"_coverage"))

  full_tibble <- dplyr::bind_cols(boots_tibble, coverage_tibble)

  return(full_tibble)



}
