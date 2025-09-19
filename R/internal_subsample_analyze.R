

#' Find minimum detectable effect size
#'
#' Identify the smallest effect size for which the proportion correct detections is >= target power for every coverage interval
#'
#' @param .prop_correct dataframe of proportion of correct detections for every effect size bin at every coverage interval, output of `.census_proportion_correct()`.
#' @return list of length `length(power)`. Each element is a tibble of coverage ranks and the associated smallest absolute effect size for which detectability was >= power.
#' @noRd
.find_min_detectable_effect <- function(.prop_correct, .power){

  # plot the minimum effect size that can be detected with a given power
  # at each step of sample coverage along coverage_seq
  min_detectable <- purrr::map(
    .x = .power %>% purrr::set_names(.power),
    .f = ~ .prop_correct %>%
      # make binary column indicating whether the proportion correct is above power threshold
      dplyr::mutate(threshold = prop_correct >= .x) %>%
      # group by coverage value
      dplyr::group_by(coverage, coverage_rank) %>%
      # filter to only values more correct than the threshold
      dplyr::filter(threshold == T) %>%
      # arrange by abs_eff_size to find the smallest effect size
      # which we can detect at the given coverage value
      dplyr::arrange(abs_eff_size) %>%
      # take the first
      dplyr::slice(1) %>% dplyr::ungroup() %>%
      dplyr::mutate(power = .x) %>%
      # remove the coverage = 100 point beause it necessarily detects correctly.
      dplyr::filter(coverage != 1)
  )

  return(min_detectable)

}


#' Model detectable effect size over coverage
#'
#' Linear model the (absolute) minimum detectable effect size at every interval of coverage.
#'
#' @param .min_detectable list of tibbles of minimum detectable effect sizes at every coverage interval. Output of `.find_min_detectable_effect()`.
#' @return list of linear model objects for every power.
#' @noRd
.model_min_detectable <- function(.min_detectable){

  # make a linear model of minimum detectable effect size as a function of coverage
  # note that here, we are using the linearly scaled "coverage rank" (0-41), which
  # actually represents coverage on a non-linear scale (coverage_seq)
  power_mod <- purrr::map(
    .x = .min_detectable,
    .f = ~lm(data = .x,
             formula = abs_eff_size ~ coverage_rank)
  )

}


#' Augment effect size models
#'
#' Predict linear model outputs over all intervals of coverage.
#'
#' @param .power_mod list of linear models of detectable effect size across coverage intervals at every power. Output of `.model_min_detectable()`.
#' @return list of length `length(power)`, where each item is a tibble of predicted model outputs at coverage value 0-40
#' @noRd
.augment_power_mods <- function(.power_mod, .min_detectable){

  power_au <- purrr::map2(
    .x = .power_mod,
    .y = .min_detectable,
    .f = ~broom::augment(.x,
                         newdata = .y %>% tidyr::complete(
                           data.frame(coverage_rank = c(0:40)) # complete to 0-40
                         ),
                         interval = "confidence")
  ) %>%
    dplyr::bind_rows(.id = "power")

  return(power_au)

}
