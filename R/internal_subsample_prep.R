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
        .f = ~find_eff_sites(coverage_vec = .x,
                                             target = .coverage_seq*max(.x)
        )
      ),
      !!(new_cols[2]) := purrr::map(
        .x = get(old_cols[4]),
        .f = ~find_eff_sites(coverage_vec = .x,
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
    dplyr::mutate(coverage_rank = 1:(dplyr::n())) %>%
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
.census_proportion_correct <- function(.df, .analysis_type){


  # find 1/2 bin width (center of smallest bin / 2)
  if(.analysis_type != "sign"){
    vec <- unique(.df$eff_size_num)
    half <- min(vec[vec > 0])/2 # find half of smallest positive bin
    half2 <- abs(max(vec[vec<0])/2) # find half of negative smallest bin
    stopifnot(half == half2) # make sure they're the same
  }

# .df %>% dplyr::group_by(eff_size_num) %>%
#   dplyr::select(eff_size_num, coverage, log2diff2, coverage_rank) %>%
#   dplyr::filter(coverage_rank > 36,
#                 eff_size_num > 0) %>%
#   dplyr::mutate(
#     eff_size_num_high = eff_size_num + half,
#     eff_size_num_low = eff_size_num - half
#   ) %>%
#   dplyr::slice(1:4) %>%
#   dplyr::mutate(correct = dplyr::case_when(
#     eff_size_num < 0 ~ ifelse(log2diff2 < 0 & log2diff2 > (eff_size_num - half),T,F),
#     eff_size_num > 0 ~ ifelse(log2diff2 > 0 & log2diff2 > (eff_size_num + half),T,F),
#     eff_size_num == 0 ~ ifelse(abs(log2diff2) == 0,T,F))
#   )

 #df %>% filter(eff_size_num ==vec[1]) %>%
 #  filter(coverage_rank %in% c(1,10,20,30,40)) %>%
 #  dplyr::mutate(correct = dplyr::case_when(
 #    eff_size_num < 0 ~ ifelse(log2diff2 < 0 & log2diff2 >= (eff_size_num - half),T,F),
 #    eff_size_num > 0 ~ ifelse(log2diff2 > 0 & log2diff2 <= (eff_size_num + half),T,F),
 #    eff_size_num == 0 ~ ifelse(log2diff2 <= half & log2diff2 >= half*-1,T,F))
 #  ) %>%
 #  ggplot() +
 #  geom_vline(aes(xintercept = eff_size_num - half ), linewidth = .2 )+
 #  geom_vline(aes(xintercept = eff_size_num + half ), linewidth = .2 )+
 #  geom_vline(xintercept = 0, linetype = "dotted")+
 #  geom_histogram(aes(x = log2diff2,
 #                     fill = correct),
 #                 boundary = 0,
 #                 binwidth = half) +
 #  geom_rect(data = . %>% group_by(coverage_rank) %>% slice(1),
 #            aes(xmin = eff_size_num - half,
 #                xmax = eff_size_num + half,
 #                ymin = -Inf,
 #                ymax = Inf),
 #            fill = "grey80",
 #            alpha = .5)+
 #
 #  geom_vline(aes(xintercept = eff_size_num - half ), linewidth = .1 )+
 #  geom_vline(aes(xintercept = eff_size_num + half ), linewidth = .1 )+
 #
 #  facet_wrap(~coverage_rank) +
 #  coord_cartesian(xlim = c(-.5,.5))

  # find proportion correct based on analysis type
  prop_correct1 <- switch(
    .analysis_type,

    # first, if analysis is type "sign": --------
    "sign" = .df %>%
      # first calculate correct sign, except in zero effect size bin
      dplyr::mutate(correct = ifelse(sign(log2diff2) == sign(eff_size_num),
                                                          T,F)
                      #dplyr::case_when(
                      #  # when eff_size is not zero, calculate whether detected change is same
                      #  # sign as true change
                      #  eff_size_num != 0 ~ ifelse(sign(log2diff2) == sign(eff_size_num),
                      #                             T,F),
                      #  # when eff size is zero, calculate whether log2 def is also zero
                      #  eff_size_num == 0 ~ ifelse(abs(log2diff2) == 0,
                      #                             T,F))
                      ),

    # option 2: at least
    # detections will be considered correct when they:
    # - are the same sign as true difference
    # - are less than or equal to the magnitude of the true difference,
    # so that when you detect a difference,
    # you know the difference is at least that big.
    "at_least" = .df %>%
      dplyr::mutate(correct = dplyr::case_when(
        eff_size_num < 0 ~ ifelse(log2diff2 < 0 & log2diff2 >= (eff_size_num - half),T,F),
        eff_size_num > 0 ~ ifelse(log2diff2 > 0 & log2diff2 <= (eff_size_num + half),T,F),
        eff_size_num == 0 ~ ifelse(log2diff2 <= half & log2diff2 >= half*-1,T,F))
      ),

    # OPTION 3:
    # in this option, detections will be considered correct if they:
    # - are inside the bounds of their effect size bin
    "inside" = .df %>%
      dplyr::mutate(correct = ifelse(log2diff > eff_size_num - half &
                                       log2diff < eff_size_num + half, T, F))
  )

 # prop_correct1 %>% mutate(abs_eff_size = abs(eff_size_num)) %>%
 #   filter(abs_eff_size == .166) %>%
 #   filter(coverage_rank > 30) %>%
 #   dplyr::select(eff_size_num, coverage, log2diff2, coverage_rank, correct) %>%
 #   dplyr::mutate(
 #     eff_size_num_high = eff_size_num + half,
 #     eff_size_num_low = eff_size_num - half
 #   ) %>%
 #   print(n = Inf)

  prop_correct2 <- prop_correct1 %>%

    dplyr::mutate(abs_eff_size = abs(eff_size_num)) %>%
    dplyr::group_by(abs_eff_size, coverage_rank, coverage,sec_axis_labels) %>%
    dplyr::summarize(n_correct  = sum(correct),
                     n_total = dplyr::n(),
                     .groups = "drop",
                     across(starts_with("pilot_n_sites"), unique, .names = "{.col}"),
                     include = unique(include)) %>%
    dplyr::mutate(prop_correct = (n_correct/n_total)*100) %>%
    dplyr::ungroup()


  return(prop_correct2)

}



.get_out_df <- function(.pilot,
                        .true,
                        .pilot_minimum_detectable,
                        .target_ann = NULL,
                        .cost_per_sample,
                        .target_eff_size,
                        .method){


  # first, arrange the achieved values
  out.achieved <- .pilot_minimum_detectable %>%
    dplyr::mutate(group = "achieved") %>%
    dplyr::relocate(group) %>%
    dplyr::rename(min_detectable_effect = achieved_min_eff_size,
                  coverage = pilot_achieved_cover,
                  sample_size = pilot_ss) %>%
    dplyr::select(-pilot_achieved_rank)

  if(!is.null(.true)){

    out.achieved <- switch(
      .method,
      "single" =  out.achieved %>%
        dplyr::bind_cols(.true %>% dplyr::select(raw_richness)) %>%
        dplyr::rename(sample_size.total = sample_size),
      "two" = out.achieved %>%
      dplyr::select(-sample_size) %>%
      dplyr::bind_cols(.true %>%
                         dplyr::select(-min_coverage)) %>%
      dplyr::mutate(raw_sample_size.total =
                      sum(dplyr::across(dplyr::contains("raw_sample_size."))),
                    .before = dplyr::contains("rarefied_richness")) %>%
      dplyr::mutate(rarefied_sample_size.total =
                      sum(dplyr::across(dplyr::contains("rarefied_sample_size."))),
                    .before = true_abs_effect)
    )

  }

#  pilot_ss <- switch(
#    class(.pilot[[2]]),
#    "numeric" = data.frame("sample_size.total" = nrow(.pilot)),
#    "character" = dplyr::count(.pilot, across(2)) %>%
#    tidyr::pivot_wider(values_from = n,
#                       names_from = 1,
#                       names_prefix = "sample_size.") %>%
#    dplyr::mutate(sample_size.total = sum(dplyr::c_across(starts_with("sample_size.")), na.rm = TRUE)))
#
#  out.achieved <- out.achieved %>%
#    dplyr::select(-sample_size) %>%
#    dplyr::bind_cols(pilot_ss)
#
  out.df <- out.achieved

  # if there is a target effect size, calculate "target" values
  if(!is.null(.target_eff_size)){
    out.target <- .target_ann %>%
      dplyr::mutate(group = "target", .before = power ) %>%
      dplyr::rename(min_detectable_effect = target_eff_size,
                    coverage = coverage_value
      ) %>%
      dplyr::select(-coverage_rank, -ann, -dplyr::contains("cost"))

    # rename to match other columns (with treatment names) if it's a two-treatment
    if(.method == "two"){

      colnames(out.target)[stringr::str_detect(colnames(out.target),"sample_size")] <-
        stringr::str_replace(    colnames(out.target)[stringr::str_detect(colnames(out.target),"sample_size")],
                                 "sample_size","rarefied_sample_size")

    }

  out.df <- out.achieved %>% dplyr::bind_rows(out.target)
  }


  # add cost if applicable
  if(!is.null(.cost_per_sample)){
    out.df <- switch(
      .method,
      "single" = out.df %>% dplyr::mutate(total_cost = sample_size.total*.cost_per_sample),
      "two" = out.df %>%
      dplyr::mutate(total_cost.raw = raw_sample_size.total*.cost_per_sample,
                    total_cost.rarefied = rarefied_sample_size.total*.cost_per_sample)
    )
    }

  return(out.df)
}
