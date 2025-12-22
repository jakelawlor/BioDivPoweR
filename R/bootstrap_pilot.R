#' Bootstrap Biodiversity Samples from Species-by-Site Matrix
#'
#' Creates `n_boots` simulated communities from a pilot biodiversity survey, pulled from randomly drawn sites from the `pilot` survey, with replacement. Bootstrapped pairs are compared for differences in richness, and differences in richness (log2 ratio) are sorted into a histogram with `n_eff_size_bins` bins. When possible, `min_exp_n` experiments (pairs of simulated communities, rarefied to equal sample coverage) are retained and kept as the input for the following function, `subsample_boots()`.
#'
#' @param pilot species-by-site matrix from a pilot biodiversity survey, either as a single treatment (matrix of multiple samples within one study area), or with two treatments (e.g., high vs. low vegetation sites, restored vs. unrestored), with a column that specifies the treatment type. Single-treatment analysis assesses the power to detect richness change within one study area over time, and two-treatment assesses the power to detect change between two site types.
#' @param method power analysis type, `"single"` for sites with no treatment levels, `"two"` for a two-treatment analysis.
#' @param category_col the column name in `pilot` that specifies the treatment value of sampling sites, if two-treatment analysis is selected.
#' @param n_boots The number of simulated community pairs to bootstrap, which will be compared to each other for differences in species richness. Increasing the number of bootstrapped communities should result in more effect sizes (differences in richness between simulated pairs) qualifying for the power analysis (a histogram of richness differences from `n_boots` communities will be sorted into `n_eff_size_bins`, and bins that surpass the `min_exp_n` threshold will be retained for following steps).
#' @param n_eff_size_bins n_eff_size_bins The total number of bins to separate the histogram of  richness differences detected between simulated communities. Increasing `n_eff_size_bins` can offer higher resolution in the following steps (more bins to qualify), but can result in "toothy" histograms in which sequential bins are not all filled. Decreasing `n_eff_size_bins` should fix "toothiness".
#' @param min_exp_n The threshold number of simulated pairs to qualify a given effect size bin for retention. Higher `min_exp_n` values will offer greater resolution in following steps, but will decrease the number of effect sizes that qualify.
#' @param seed Random seed. Defaults to 1 so runs of the same data will provide the same answers, but since the simulations all rely on random draws, changing the seed will result in different answers.
#'
#' @returns a tibble of `min_exp_n` simulated communities, rarefied to equal coverage, within all effect size bins that qualify. Additionally, prints a histogram of simulated community richness values, highlighting the bins that are kept for the next step.
#' @export
#'
#' @examples if(FALSE){boostrap_pilot(data("pilot_single_trt"))}
bootstrap_pilot <- function(pilot,
                            method = "single",
                            category_col = NULL,
                            n_boots = 5000,
                            n_eff_size_bins = 40,
                            min_exp_n = 40,
                            seed = NULL){

  if(method == "two" & is.null(category_col))
    stop("Please specify the name of the treatment column in pilot dataset")
  if(method == "two")
    if(length(unique(pilot[[category_col]])) != 2)
      stop("Two-treatment pilot must have 2 levels")
  if(method == "single" & sum(sapply(pilot,is.character)) > 1)
    stop("Single-treatment analyses should have only one character column")


  # set random seed
  set.seed(seed)

  # create bootstrap arrays -------------------------------------------------
  cat("Bootstrapping empirical communities", n_boots,"times... \n")
  boot_arrays <- .create_boot_arrays(pilot,
                                     method = method,
                                     category_col = category_col,
                                     n_boots = n_boots)

  # rarefy pairs to equal coverage -----------------------------------------
  # cut one of each of n_boots community pairs to the number of sites
  # that gives the two equal maximum coverage.
  cat("Rarefying", n_boots, "bootstraps to equal coverage... \n")

  # find coverage in bootstrapped communities
  eff_coverage <- .find_rarefied_coverage(boot_arrays)

  rarefied_boots <- .rarefy_boots(boot_arrays, eff_coverage)
  rm(boot_arrays) #remove the full-length boots

  # summarize richness differences between rarefied bootstrapped community pairs
  summary <- .summarize_boot_diff(rarefied_boots)

  # remove true zeros:
  zeros <- which(summary$log2_rich_diff == 0)
  rarefied_boots[[1]] <- rarefied_boots[[1]][-zeros]
  rarefied_boots[[2]] <- rarefied_boots[[2]][-zeros]
  eff_coverage[[1]] <- eff_coverage[[1]][-zeros]
  eff_coverage[[2]] <- eff_coverage[[2]][-zeros]
  summary <- summary[summary$log2_rich_diff != 0,]

  # make bootstrap plots ----------------------------------------------------
  cat("Sorting differences to",n_eff_size_bins,"effect size bins... \n")

  # plot 1: histogram of raw rarefied richness values in community 1 and community 2
  p1 <- .plot_boot_values(summary)
  print(p1)

  # plot 2: histogram of raw difference in richness between
  #         equal-coverage bootstrapped community pairs
  p2 <- .plot_boots_raw_diff(summary)
  print(p2)

  # plot 3: log2 richness differences between bootstrapped pairs
  #         and color histogram bins with counts >= min_exp_n
  p3 <- .plot_boots_log2_diff(summary, min_exp_n = min_exp_n, n_eff_size_bins = n_eff_size_bins)
  print(p3)

  # export boots from qualifying bins -----------------------------------------
  cat("Saving", min_exp_n, "boostraps from each qualifying bin... \n")

  # extract bootstrap IDs to export
  ids_to_keep <- .extract_ids_to_keep(p3, summary, min_exp_n)

  # pull keep id bootstraps and coverage curves filtered to 1:n_sites_rarefied
  rarefied_boots_tibble <- .keep_selected_boots(ids_to_keep, rarefied_boots, eff_coverage)
  rm(rarefied_boots)

  # return the tibble
  return(rarefied_boots_tibble)

}
