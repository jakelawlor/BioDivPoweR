#' Subsample bootstrapped communities to assess power for detecting change across community coverage
#'
#' Downsample bootstrapped communities to lower sample sizes representing a gradient of coverage values, and assess the proportion of correct-direction detections within each coverage interval.
#'
#' @param data tibble of bootstrapped communities across multiple effect sizes. Output from `bootstrap_pilot()`.
#' @param pilot original pilot study. Identical to the pilot input in `bootstrap_pilot()`
#' @param power statistical power with which the user wants to detect richness change (defaults to 80).
#' @param target_eff_size if supplied, the level of richness difference the user wants to detect (in log2 ratio units). Adding a target value will change outputs of the function, specifying sampling needs to reach the richness difference target.
#' @param cost_per_sample cost per unit sample (linear only), which will add an axis to output plots for total cost per unit power.
#' @param seed random seed. Defaults to 1 so repeat runs will be identical, but since simulations rely on random draws, changing the seed will result in different answers.
#'
#' @returns List of outputs: (1) Plot of convergence of detections across community coverage in multiple effect size bins. (2) Plot of mean detection across community coverage within multiple effect size bins. (3) Plot of proportion of correct detections across community coverage within multiple effect size bins. (4) Plot of functional relationship between sample size and power to detect richness differences. (5) Tibble of input and output values.
#' @export
#'
#' @examples if(FALSE){subsample_boots(boots, pilot)}
subsample_boots <- function(boots,
                           pilot,
                           method = "single",
                           power = c(80),
                           target_eff_size = NULL,
                           cost_per_sample = NULL,
                           seed = NULL){

  # various checks
  if(method == "single" & sum(sapply(pilot, is.character)) > 1)
    stop(call.=F,"Singe-treatment analysis selected, but pilot contains multiple character columns. Please check.")
  if(method == "two" & sum(sapply(pilot, is.character)) == 1)
    stop(call.=F,"Two-treatment analysis selected, but pilot appears to be missing a category column. Please check.")

  # set random seed
  set.seed(seed)

  # set waypoints across the coverage scale (x-axis)
  coverage_seq <- c(1,seq(from = 250, to = 10000, by = 250) %>% sqrt())/100

  # get community type names
  community_types <-
    colnames(boots)[stringr::str_detect(colnames(boots),".matrix")] %>%
    stringr::str_remove(.,".matrix")


  # subsample across coverages ---------------------------------------------
  subsamples <- .downsample_bootstraps(boots, coverage_seq, community_types)
  rm(boots)


  # summarize downsamples ---------------------------------------------------
  cat("Calculating richness differences between subsamples...\n")
  # this reformats data in a long df, suitable for plotting
  df <- .summarize_subsamples_diff(subsamples, coverage_seq, community_types)
  rm(subsamples)

  ## calculate n sites in pilot for each coverage step ==============================
  pilot_coverage <- .find_pilot_coverage(pilot, method, coverage_seq)
  sec_axis_labels <- .calc_second_axis(pilot_coverage, coverage_seq, cost_per_sample, method, community_types)


  df <- df %>% dplyr::left_join(sec_axis_labels, dplyr::join_by(coverage, coverage_rank))
  rm(sec_axis_labels)


  # plot --------------------------------------------------------------------

  # funnel plot of effect size across coverage values
  p0 <- .plot_subsample_funnel(df, coverage_seq)

  # line plot of mean detected richness difference across coverages
  p1 <- .plot_subsample_mean(df, coverage_seq)

  # summarize proportion correct --------------------------------------------
  cat("Calculating proportion correct detection across subsamples...\n")
  # find proportion correct detections (same sign, or == 0 in zero bin)
  prop_correct <- .census_proportion_correct(df)

  p2 <- .plot_prop_correct(prop_correct, power, coverage_seq)

  # analyze trends ----------------------------------------------------------
  cat("Modeling detectable richness as a function of sample coverage...\n")

  # find minimum detectable effect size
  # (smallest absolute effect size in which prop_correct was >0 power at every coverage interval).
  min_detectable <- .find_min_detectable_effect(prop_correct, power)

  # model minimum detectable effect size as a function of coverage
  power_mod <- .model_min_detectable(min_detectable)

  # augment power model to full coverage gradient
  power_au <- .augment_power_mods(power_mod, min_detectable)



  # final plot --------------------------------------------------------------

  ## make base plot ==============================
  # regression of minimum detectable effect across powers
  p3 <- .plot_power_base(min_detectable, power_au, coverage_seq)


  ## add achieved coverage ==============================
  # since we know the number of samples in the original pilot community,
  # we can back-calculate coverage, then predict the achieved
  # minimum detectable effect size from the pilot study that we conducted
  pilot_minimum_detectable <- .calc_minEff_at_sampleSize(pilot, pilot_coverage, coverage_seq, power_mod, method)
  p3 <- .add_minEff_at_sampleSize(p3, pilot_minimum_detectable)

  # add true effect size =================================
  # add a point along the y axis showing the probably presumed effect size
  # that could exist in our communities:
  # richness-to-rarified-richness mean for two communities
  if(method == "two"){
    true <- .calc_true_effect(pilot, pilot_coverage)
    p3 <- .add_true_effect(true, p3)
  }
  # or, for single community, just find raw sample size, richness, and coverage values of pilot
  if(method == "single"){
    true <- .calc_raw_vals_single(pilot, pilot_coverage)
  }

  # add samples for target effect size ========
  target_ann <- NULL
  if(!is.null(target_eff_size)){
    target_ann <- .calc_site_to_reach_target(target_eff_size, power_au, pilot_coverage, coverage_seq, cost_per_sample, method)
    p3 <- .add_target_to_plot(p3, target_ann)
  }


  # # save answers as df ------------------------------------------------------
  # make a dataset of all relevant takeaways:

  # get relevant values as dataframe
  out.df <- .get_out_df(.pilot = pilot,
                        .true  = true,
                        .pilot_minimum_detectable = pilot_minimum_detectable,
                        .target_ann = target_ann,
                        .cost_per_sample = cost_per_sample,
                        .target_eff_size = target_eff_size,
                        .method = method)

  # add to plot
  out <- list(p0, p1, p2, p3, out.df)

  return(out)

}







