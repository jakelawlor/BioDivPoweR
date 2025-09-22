#' Subsample bootstrap matrices
#'
#' downsample bootstrap community pairs to test for differences at coverages in coverage_seq
#'
#' @param .boots tibble of rarefied bootstraps. output of `bootstrap_pilot()` function
#' @param .coverage_seq custom sequence of coverage values
#' @param .community_types names of test communities, gathered programmatically at start of `subsample_boots`
#' @return .boots tibble with 4 extra columns: {group1}.n_to_check, {group2}.n_to_check, {group1}.richness, {group2}.richness
#' @noRd
.downsample_bootstraps <- function(.boots, .coverage_seq, .community_types){


  # make names for new columns to save ourselves coding below
  old_cols <- colnames(.boots)[c(-1)]
  new_cols <- c(paste0(.community_types,".n_to_check"), paste0(.community_types,".richness"))

  # find n sites to sample ------
  # first, find the number of sites that correspond to each of the target coverages.
  # here, we apply find_eff_sites() to each bootstrap community's coverage vector,
  # but instead of .01-1, the targets are scaled to the max coverage value for the communty
  # so that 100% cover will be the full number of sites.
  boots2 <- .boots %>%
    dplyr::mutate(
      !!(new_cols[1]) := purrr::map(
        .x = get(old_cols[3]),
        .f = ~BioDivSampler:::find_eff_sites(coverage_vec = .x,
                                             target = .coverage_seq*max(.x)
        )
      ),
      !!(new_cols[2]) := purrr::map(
        .x = get(old_cols[4]),
        .f = ~BioDivSampler:::find_eff_sites(coverage_vec = .x,
                                             target = .coverage_seq*max(.x))
      ),
    )
  # this adds two columns: {group1}.n_to_check and {group2}.n_to_check
  # which list the numbers of samples required to hit every value in coverage_seq
  rm(.boots) # remove previous


  # sample richness at each n_to_check
  boots3 <- boots2 %>%
    dplyr::mutate(
      !!(new_cols[3]) := purrr::map2(
        .x = get(old_cols[1]), # community 1 matrix
        .y = get(new_cols[1]), # community 1 n to check
        .f = ~ {
          # explicitly name the arguments for clarity
          comm_matrix <- .x
          sample_sizes <- .y

          purrr::map_dbl(sample_sizes, ~ {
            n <- .x # current sample size
            # downsample the matrix to 'n' rows - without replacement
            sampled_matrix <- comm_matrix[sample(nrow(comm_matrix), n, replace = FALSE),, drop  = F]

            # find sum of species present more than zero times
            sum(colSums(sampled_matrix) > 0)

          })
        }),
      !!(new_cols[4]) := purrr::map2(
        .x = get(old_cols[2]), # community 1 matrix
        .y = get(new_cols[2]), # community 1 n to check
        .f = ~ {
          # explicitly name the arguments for clarity
          comm_matrix <- .x
          sample_sizes <- .y

          purrr::map_dbl(sample_sizes, ~ {
            n <- .x # current sample size
            # downsample the matrix to 'n' rows - without replacement
            sampled_matrix <- comm_matrix[sample(nrow(comm_matrix), n, replace = FALSE),, drop  = F]

            # find sum of species present more than zero times
            sum(colSums(sampled_matrix) > 0)

          })
        })
    )
  # this adds two columns to the tibble:
  # {group1}.richness, and {group2}.richness,
  # which list the richness observed at samples
  # of n_to_check sites of each community.
  rm(boots2) # remove previous

  return(boots3)


}

#' Summarize sub-sampled bootstraps
#'
#' Find richness differences between downsampled communities at increasing coverage intervals
#'
#' @param .subsamples tibble of downsampled bootstraps. output of `.downsample_bootstraps()` function
#' @param .coverage_seq custom sequence of coverage values
#' @return dataframe of richness effect sizes from coverage .01-1 for `min_exp_n` community pairs in all qualifying effect size bins
#' @noRd
.summarize_subsamples_diff <- function(.subsamples, .coverage_seq, .community_types){


  # find differences in richness between community 1 and community 2
  # at sample sizes that represent each coverage step.
  df <- .subsamples %>%
    # keep only effect size, and richness vectors at the coverage waypoints
    dplyr::select(eff_size, dplyr::contains(".richness")) %>%
    dplyr::group_by(eff_size) %>%
    # add trial (which will be the number of repeats within each effect size bin)
    dplyr::mutate(trial = c(1:dplyr::n())) %>%
    # unnest list columns
    tidyr::unnest(dplyr::contains(".richness")) %>%
    dplyr::group_by(eff_size,trial) %>%
    dplyr::mutate(coverage = .coverage_seq) %>%
    dplyr::mutate(coverage_rank = 0:(dplyr::n()-1)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(log2diff = log2(
      get(paste0(.community_types[[2]],".richness")) / get(paste0(.community_types[[1]],".richness"))
    )) %>%
    # note that in some cases, log2diff comes out as infinite if
    # numerator community is zero, and comes out as NA if denominator
    # community richness is zero. Manually fix those here.
    dplyr::mutate(log2diff2 = dplyr::case_when(log2diff == -Inf ~ min(log2diff[is.finite(log2diff)], na.rm=T),
                                               log2diff == Inf ~ max(log2diff[is.finite(log2diff)], na.rm=T),
                                               is.na(log2diff) ~ 0,
                                               TRUE ~ log2diff)) %>%
    # note that in rare cases, log2diff is Inf when comm2_rich == 0,
    # not sure what we should do with this in the future, but let's
    # manually change for now.
    dplyr::mutate(eff_size_num = as.numeric(eff_size))
  rm(.subsamples)

  return(df)

}

#' Census proportion correct
#'
#' Find proportion of `min_exp_n` (or `min_exp_n`*2 for absolute valued bins) in all qualifying effect size bins across coverage values
#'
#' @param .df long-format tibble of detected changes in every qualifying experiment across coverage values
#' @return summarized df displaying number and proportion of correct-direction detections of all experiments within effect size bins and coverage intervals
#' @noRd
.census_proportion_correct <- function(.df){

  prop_correct <-  .df %>%
    # first calculate correct sign, except in zero effect size bin
    dplyr::mutate(correct_sign =
                    dplyr::case_when(
                      # when eff_size is not zero, calculate whether detected change is same
                      # sign as true change
                      eff_size_num != 0 ~ ifelse(sign(log2diff2) == sign(eff_size_num),
                                                 T,F),
                      # when eff size is zero, calculate whether log2 def is also zero
                      eff_size_num == 0 ~ ifelse(abs(log2diff2) == 0,
                                                 T,F))) %>%
    dplyr::mutate(abs_eff_size = abs(eff_size_num)) %>%
    dplyr::group_by(abs_eff_size, coverage_rank, coverage) %>%
    dplyr::summarize(n_correct  = sum(correct_sign),
                     n_total = dplyr::n(),
                     .groups = "drop") %>%
    dplyr::mutate(prop_correct = (n_correct/n_total)*100) %>%
    dplyr::ungroup()

  return(prop_correct)

}



.get_out_df <- function(.pilot,
                        .pilot_minimum_detectable,
                        .target_ann = NULL,
                        .cost_per_sample,
                        .target_eff_size){

  # first, arrange the achieved values
  out.achieved <- .pilot_minimum_detectable %>%
    dplyr::mutate(group = "achieved") %>%
    dplyr::relocate(group) %>%
    dplyr::rename(min_detectble_effect = achieved_min_eff_size,
                  coverage = pilot_achieved_cover,
                  sample_size = pilot_ss) %>%
    dplyr::select(-pilot_achieved_rank)


  pilot_ss <- switch(
    class(.pilot[[2]]),
    "numeric" = data.frame("sample_size.total" = nrow(.pilot)),
    "character" = dplyr::count(.pilot, across(2)) %>%
    tidyr::pivot_wider(values_from = n,
                       names_from = 1,
                       names_prefix = "sample_size.") %>%
    dplyr::mutate(sample_size.total = sum(dplyr::c_across(starts_with("sample_size.")), na.rm = TRUE)))

  out.achieved <- out.achieved %>%
    dplyr::select(-sample_size) %>%
    dplyr::bind_cols(pilot_ss)

  out.df <- out.achieved

  if(!is.null(.target_eff_size)){
    out.target <- .target_ann %>%
      dplyr::mutate(group = "target", .before = power ) %>%
      dplyr::rename(min_detectble_effect = target_eff_size,
                    coverage = coverage_value
      ) %>%
      dplyr::select(-coverage_rank, -ann)
    out.df <- out.achieved %>% dplyr::bind_rows(out.target)
  }


  # add cost if applicable
  if(!is.null(.cost_per_sample)){
    out.df <- out.df %>%
      dplyr::mutate(total_cost = sample_size.total*.cost_per_sample)
  }

  return(out.df)
}
