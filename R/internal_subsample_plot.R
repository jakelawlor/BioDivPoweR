# Helper functions for plotting within subsample_boots function

#' Plot subsample funnels
#'
#' Make funnel plots of detected richness differences in bootstrapped communities, subsampled at a gradient of coverage values.
#'
#' @param df dataframe of richness differences between boots at all values in `coverage_seq`. Output of `.summarize_subsamples_diff()`
#' @return ggplot of effect size as coverage increases for `min_exp_n` experiments in all qualifying effect size bins
#' @noRd
.plot_subsample_funnel <- function(.df, .power, .coverage_seq){

  # plot funnels within effect size bins
  p0 <-   .df %>%
    dplyr::group_by(eff_size, coverage_rank,eff_size_num) %>%
    dplyr::summarize(quant10 = quantile(log2diff2, .1),
                     quant50 = quantile(log2diff2, .5),
                     quant90 = quantile(log2diff2, .9),
                     .groups = "drop") %>%
    dplyr::mutate(eff_size = forcats::fct_reorder(eff_size, eff_size_num)) %>%
    ggplot2::ggplot(ggplot2::aes(x = coverage_rank,
                                 y = quant50)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = quant10,
                                      ymax = quant90),
                         fill = "grey80",
                         color = "grey20") +
    ggplot2::geom_line() +
    ggplot2::geom_hline(yintercept = 0, linewidth = .4) +
    ggplot2::geom_hline(ggplot2::aes(yintercept = as.numeric(as.character(eff_size))),
                        color = "darkred", linetype = "dashed") +
    ggplot2::facet_wrap(~eff_size) +
    ggplot2::coord_cartesian(ylim = c(-1,1)) +
    ggplot2::scale_x_continuous(
      breaks = seq(0,40,by=4),
      labels = round(c(.coverage_seq*100),0)[c(T,F,F,F)]
    ) +
    ggplot2::labs(x = "Coverage (%)",
                  y = "Detected Richness Change",
                  title = "Convergence across coverage within effect size bins") +
    ggplot2::theme_bw()

  return(p0)

}

#' Plot subsample means across coverage intervals
#'
#' Make line plot of mean effect sizes of all experiments within effect size bins across coverage intervals.
#'
#' @param .df dataframe of richness differences between boots at all values in `coverage_seq`. Output of `.summarize_subsamples_diff()`
#' @return ggplot of mean effect size. should converge to effect size bin values.
#' @noRd
.plot_subsample_mean <- function(.df, .coverage_seq){

  # plot mean richness change across coverage within bins
  p1 <- .df %>%
    dplyr::group_by(eff_size_num, coverage, coverage_rank) %>%
    dplyr::summarize(mean_rich = mean(log2diff2),
                     .groups = "drop") %>%
    dplyr::ungroup() %>%
    ggplot2::ggplot(ggplot2::aes(x = coverage_rank,
                                 y = mean_rich,
                                 group = eff_size_num,
                                 color = abs(eff_size_num))) +
    ggplot2::geom_line(linewidth = 1) +
    ggplot2::geom_point(data = . %>% dplyr::filter(coverage == 1),
                        shape = 21, fill ="white",
                        stroke = 1) +
    ggplot2::scale_color_gradientn(colors = BioDivSampler:::cool_matlab()) +
    ggplot2::scale_x_continuous(
      breaks = seq(0,40,by=2),
      labels = round(c(.coverage_seq*100),0)[c(T,F)]
    ) +
    ggplot2::labs(x = "Coverage (%)",
                  y = "Detected Richness Difference\n(log2 fold-difference)",
                  color = "Absolute Effect Size",
                  title = "Mean across coverage within effect size bins") +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position = "inside",
                   legend.position.inside = c(.05,.05),
                   legend.justification = c(0,0),
                   legend.key.width = ggplot2::unit(25,"pt"),
                   legend.key.height = ggplot2::unit(6,"pt"),

                   legend.direction = "horizontal",
                   legend.background = ggplot2::element_rect(color = "black",
                                                             linewidth = .1)) +
    ggplot2::guides(color = ggplot2::guide_colorbar(theme = ggplot2::theme(legend.title.position = "top"))) +
    ggplot2::theme(panel.grid = ggplot2::element_blank())

  return(p1)

}


#' Plot subsample means across coverage intervals
#'
#' Make line plot of mean effect sizes of all experiments within effect size bins across coverage intervals.
#'
#' @return ggplot of mean effect size. should converge to effect size bin values.
#' @noRd
.plot_prop_correct <- function(.prop_correct, .power = power, .coverage_seq){

  p2 <- .prop_correct %>%
    ggplot2::ggplot(ggplot2::aes(x = coverage_rank,
                                 y = prop_correct,
                                 group = abs_eff_size,
                                 color = abs_eff_size)) +
    ggplot2::geom_line(linewidth = 1) +
    ggplot2::scale_color_gradientn(colors = BioDivSampler:::cool_matlab(),
    ) +
    ggplot2::scale_x_continuous(
      breaks = seq(0,40,by=2),
      labels = round(.coverage_seq*100)[c(T,F)]
    ) +
    ggplot2::scale_y_continuous(labels = scales::percent_format(scale = 1)) +
    ggplot2::theme_bw() +
    ggplot2::labs(x = "Coverage (%)",
                  y = "Proportion Correct Direction",
                  color = "Absolute\nEffect Size") +
    ggplot2::theme(legend.position = "right",
                   legend.key.width = ggplot2::unit(12,"pt"),
                   legend.key.height = ggplot2::unit(40,"pt"),
                   panel.grid = ggplot2::element_blank()) +
    ggplot2::guides(color = ggplot2::guide_colorbar(theme = ggplot2::theme(legend.title.position = "top"))) +
    ggplot2::geom_hline(data = data.frame(power = .power),
                        ggplot2::aes(yintercept = .power),
                        linetype = "dashed")


  return(p2)

}



#' Plot coverage by detectable effect - base only
#'
#' Plot regression line and points for detectable effect across coverage values.
#'
#' @param .min_detectable list of tibbles of minimum detectable effect sizes at every coverage interval. Output of `.find_min_detectable_effect()`.
#' @param .power_au list of length `length(power)`, where each item is a tibble of predicted model outputs at coverage value 0-40. Outout of `..augment_power_mods()`
#' @return ggplot of trend of increasing detectable effect across coverage.
#' @noRd
.plot_power_base <- function(.min_detectable, .power_au, .coverage_seq ){


  .min_detectable <- .min_detectable %>%
    dplyr::bind_rows(.id = "power")

  # plot minimum detectable effect size across coverage
  p3 <- .min_detectable %>%
    ggplot2::ggplot(ggplot2::aes(x = coverage_rank,
                                 y = abs_eff_size,
                                 group = power)) +
    ggplot2::geom_point(size = 2, alpha = .8,
                        ggplot2::aes(fill = abs_eff_size,
                                     color = power),
                        shape = 21,
                        stroke = 1) +
    #  ggplot2::geom_ribbon(data = power_au,
    #              inherit.aes=F,
    #              ggplot2::aes(x = coverage_rank,
    #                  ymin = .lower,
    #                  ymax = .upper,
    #                  group = power),
    #              alpha = .1) +

    ggplot2::geom_line(data = .power_au,
                       ggplot2::aes(x = coverage_rank,
                                    y = .fitted,
                                    group = power,
                                    color = power)) +
    ggplot2::coord_cartesian(
      ylim = c(
        # higher limit is maximum of either observed effect size or auugmented
        # effect size from model ouput
        pmax(max(abs(.power_au$.fitted)),max(abs(.min_detectable$abs_eff_size))),
        # minimum is -.05
        0),
      xlim = c(0,40),
      clip = "off"
    ) +
    ggplot2::scale_y_reverse(expand = ggplot2::expansion(mult = c(.1, 0))) +
    ggplot2::scale_fill_gradientn(colors = BioDivSampler:::cool_matlab(),
                                  limits = c(0,max(abs(.min_detectable$abs_eff_size)))) +
    ggplot2::scale_color_grey(end = .55,
                              start = 0) +
    ggplot2::scale_x_continuous(
      breaks = seq(0,40,by=2),
      labels = round(.coverage_seq*100)[c(T,F)]
    ) +
    ggplot2::theme_bw() +
    ggplot2::labs(y = "Minimum Detectable Fold Change in Richness",
                  x = "Coverage (%)",
                  fill = "Absolute\nEffect Size",
                  color = "Power") +
    ggplot2::theme(panel.grid = ggplot2::element_blank()) +
    ggplot2::geom_hline(yintercept = 0,
                        linewidth = .2)

  return(p3)

}



#' Calculate N sites for second axis
#'
#' Use coverage in the pilot to estimate n samples to reach coverage intervals in pilot system.
#'
#' @param .pilot pilot dataset
#' @param .method analysis method: "single" or "two" treatment.
#' @param .coverage_seq coverage intervals to test. imported internally.
#' @param .community_types labels for community types. imported internally.
#' @return character vector of axis 41 axis labels, corresponding to number of sites in empirical community (or communities) to reach coverage intervals.
#' @noRd
.rescale_pilot <- function(.pilot, .method, .coverage_seq, .community_types ){


  # a. find average occupancy in pilot survey
  # note that this is just mean occurrence for single-treatment,
  # but mean treatment split by site category if two-treatment.
  pilot_occupancy <- switch(
    .method,
    "single" =  colMeans(.pilot[,-c(1)]),
    "two" =  purrr::map(
      .x = .pilot %>% split(f = .pilot[[2]]),
      .f = ~colMeans(.x[,-c(1:2)])[colSums(.x[,-c(1:2)]) > 0],
      .progress = T
    )
  )

  # b. find coverage at 1:10000 samples in pilot survey
  # again, this will be one value if single-treatment,
  # and a list of two if two-treatment.
  pilot_coverage <- switch(
    .method,
    "single" = find_coverage(pilot_occupancy),
    "two" = purrr::map(
      .x = pilot_occupancy,
      .f = ~find_coverage(occupancy = .x))
  )

  # c. find target detection rate, or the coverage at which we have <.5 of a species left to find.
  # we'll call that the coverage target at which we've found all the species.
  # again, this will be done separately for the two treatments in method == two
  pilot_target_detection_rate <- switch(
    .method,
    "single" = find_target_detection_rate(.pilot),
    "two" =  purrr::map(
      .x =  .pilot %>% split(f = .pilot[[2]]),
      .f = ~find_target_detection_rate(.x))
  )

  # d. rescale pilot coverage to zero-target detection rate, instead of 0-100
  pilot_coverage_rescale <- switch(
    .method,
    "single" = pilot_coverage / pilot_target_detection_rate,
    "two" = purrr::map2(
      .x = pilot_coverage,
      .y = pilot_target_detection_rate,
      .f = ~.x/.y
    )
  )
  # this establishes a finite number of samples at which we achieve 100% coverage
  # which will be necessarily >= than the number of samples in the pilot.

  return(pilot_coverage_rescale)

}



#' Create labels for the second X-axis
#'
#' Use rescaled coverage to identify the number of samples in pilot communit(y/ies)
#' that it would take to reach each coverage interval.
#' If applicable, add cost.
#'
#' @param .pilot_coverage_rescale pilot coverage at 1-1000 samples, rescaled by half a species so it reaches 1.
#' @param .coverage_seq sequence of coverage values in our nonlinear scale
#' @param .cost_per_sample cost per unit sample
#' @param .method single or double analysis pipeline.
#' @return character vector of axis 41 axis labels, corresponding to number of sites in empirical community (or communities) to reach coverage intervals.
#' @noRd
.calc_second_axis <- function(.pilot_coverage_rescale, .coverage_seq, .cost_per_sample, .method, .community_types){
  # remove unneeded variables

  # e. find number of samples in the empirical community (with coverage rescaled)
  # to reach the coverage values in coverage_seq.
  pilot_n_sites <- switch(
    .method,
    "single" = find_eff_sites(.pilot_coverage_rescale, .coverage_seq),
    "two" = purrr::map(
      .x = .pilot_coverage_rescale,
      .f = ~find_eff_sites(.x, .coverage_seq)
    )
  )
  # returns a vector (or two vectors if method = 2) of number of sites to sample
  # to reach each coverage interval in coverage_seq.
  # NOTE: could throw a warning here or automatically recalculate if there are any NAs.

  # f. turn into second axis labels
  sec_axis_labels <- switch(
    .method,
    "single" = as.character(pilot_n_sites),
    "two" = {
      # add the category to the first label
      int <- Map(function(vec, label) {
        vec[1] <- paste0(label, ": ", vec[1])
        vec
      }, pilot_n_sites, .community_types)

      # linebreak
      paste0(int[[1]],"\n", int[[2]])
    }
  )


  # add cost, if needed ---------------------------------------------------------
  if(!is.null(.cost_per_sample)){

    sec_axis_labels <- switch(
      .method,
      "single" = {

        # calculate cost at each interval
        cost_labs <- (pilot_n_sites * .cost_per_sample)

        # add dollar signs, and original labels underneath
        sec_axis_labels <- paste0(scales::dollar(round(cost_labs,0)),"\n", pilot_n_sites)

        # make only one of every 4 costs show for spacing reasons
        sec_axis_labels[c(F,T,T,T)] <- gsub(".*\n","",sec_axis_labels[c(F,T,T,T)])

        sec_axis_labels
      },

      "two" = {

        # calculate cost labels for sum of sample numbers in both communities
        cost_labs <- ((pilot_n_sites[[1]] + pilot_n_sites[[2]]) * .cost_per_sample)

        # add dollar signs, and original labels underneath, only to one in every 4
        sec_axis_labels[seq(3,41,by = 4)] <- paste0(scales::dollar(round(cost_labs[seq(3,41,by = 4)],0)),"\n", sec_axis_labels[seq(3,41,by = 4)])

      }
    )

  }

  return(sec_axis_labels)

}



# find the achieved power at pilot_sample_size
# here, we start with pilot sample number
# - convert to coverage (lower of two communities, if applicable)
# - linearly approximate coverage to coverage rank (0-40)
# - predict regression y at x = approx rank
.calc_minEff_at_sampleSize <- function(.pilot, .pilot_coverage_rescale, .coverage_seq, .power_mod, .method){

  # 0. find sample sizes from 1 or both pilot communities to test
  ss <- switch(
    .method,
    "single" = nrow(.pilot),
    "two" = {
     .pilot %>% dplyr::count(across(2)) %>% tibble::deframe() %>% as.list
    }
  )

  # 1. find rescaled coverage value at pilot ss
  # if 1 treatment, find coverage at ss,
  # if 2 trts, find respective and choose the lower
  ss_coverage <- switch(
    .method,
    "single" = .pilot_coverage_rescale[ss],
    "two" = {
      vals <- purrr::map2(
        .x = .pilot_coverage_rescale,
        .y = .pilot %>% dplyr::count(across(2)) %>% tibble::deframe() %>% as.list,
        .f = ~.x[.y] )
      vals[which.min(vals)]}
  )


  # 2. convert coverage to linear scale
  # coverage from the pilot is between 0 and 1, but we need it to be
  # between 0 and 40, which is our "ranked" linear coverage scale.
  # use linear approximation from the points 0-40 and the values in coverage_seq.
  pilot_coverage_rank <- approx(x= .coverage_seq,
                                y=c(0:40),
                                xout = ss_coverage)$y
  # NOTE ASK EDEN: here, I'm using the rescaled coverage scale (assuming that
  # full coverage is the target detection rate. Is that right?)


  # step 2. predict achieved minimum detectable effect size
  # at the coverage rank (0-40) from the pilot study
  pilot_minimum_detectable <- purrr::map(
    .x = .power_mod,
    .f = ~predict(.x,
                  list(coverage_rank = pilot_coverage_rank))
  ) %>%
    dplyr::bind_rows(.id = "power") %>%
    dplyr::rename("achieved_min_eff_size" = `1`)

  pilot_minimum_detectable$pilot_ss <- switch(.method,
                                              "single" = ss,
                                              "two" = paste0(ss[[names(ss_coverage)]]," (",names(ss_coverage),")"))
  pilot_minimum_detectable$pilot_achieved_cover <- ss_coverage
  pilot_minimum_detectable$pilot_achieved_rank <- pilot_coverage_rank

  return(pilot_minimum_detectable)

}


# add lines to plot for achieved minimum detectable effect and sample size

.add_minEff_at_sampleSize <- function(.p3, .pilot_minimum_detectable){

  # add to plot
  p3 <- .p3 +
    # add vertical line at n samples and text
    ggplot2::geom_segment(data = .pilot_minimum_detectable,
                          stat = "unique",
                          ggplot2::aes(x =pilot_achieved_rank,
                                       xend = pilot_achieved_rank,
                                       y = -Inf,
                                       yend = Inf),
                          linetype = "longdash") +
    ggplot2::geom_text(data = .pilot_minimum_detectable,
                       ggplot2::aes(x =pilot_achieved_rank,
                                    y = 0.01,
                                    label = paste0("A = ",pilot_ss)),
                       hjust = -.1, vjust = 1, size = 3.5, fontface = "bold") +

    # add horizontal segment at effect size and text
    ggplot2::geom_segment(data = .pilot_minimum_detectable,
                          stat = "unique",
                          ggplot2::aes(x =-Inf,
                                       xend = pilot_achieved_rank,
                                       y = achieved_min_eff_size,
                                       yend = achieved_min_eff_size),
                          linetype = "longdash") +
    ggplot2::geom_text(data = .pilot_minimum_detectable,
                       ggplot2::aes(x =0,
                                    y = achieved_min_eff_size,
                                    label = paste0("A = ",round(achieved_min_eff_size,2)),
                                    color = power),
                       hjust = 0, vjust = -.5, size = 3.5, fontface = "bold",
                       show.legend = F)


  return(p3)

}


# COME BACK TO THIS
#  .add_true_effect <- function(.boots, .method){
#
#    # calculate the true effect size as a reference point along the y axis.
#
#    true <- switch(
#      .method,
#
#      # if one community, this will be one standard deviation of bootstrapped changes
#      "single" = )
#
#  }


# here we basically do the opposite as .calc_minEff_at_sampleSize
# - start with target effect size
# - approximate rank along regression line
# - convert rank to sample size
.calc_site_to_reach_target <- function(.target_eff_size, .power_au, .pilot_coverage_rescale, .coverage_seq, .cost_per_sample = cost_per_sample, .method){


  # get the fitted model values (coverage rank) that correspond to the target effect size
    coverage_rank_to_reach_target <- purrr::map(
      .x = .power_au %>% split(f = .$power),
      # interpolate the relationship between model fit and coverage rank
      # to find the coverage rank where the model fit would be target_eff_size
      .f = ~approx(.x$.fitted,
                   .x$coverage_rank,
                   .target_eff_size) %>%
        tibble::as_tibble() %>%
        dplyr::rename(target_eff_size = x,
                      coverage_rank = y)
    ) %>%
      dplyr::bind_rows(.id = "power")


    # now we have the target coverage rank.
    # Interpolate coverage ranks against coverage sequence to get target coverage value
    coverage_value_to_reach_target <- purrr::map(
      .x = coverage_rank_to_reach_target %>% split(.$power),
      .f = ~approx(c(0:40),
                   .coverage_seq,
                   .x$coverage_rank) %>%
        tibble::as_tibble() %>%
        dplyr::rename(coverage_rank = x,
                      coverage_value = y)
    ) %>%
      dplyr::bind_rows(.id = "power")


    # re-establish number of sites to reach each rescaled coverage value
    pilot_n_sites <- switch(
      .method,
      "single" = find_eff_sites(.pilot_coverage_rescale, .coverage_seq),
      "two" = purrr::map(
        .x = .pilot_coverage_rescale,
        .f = ~find_eff_sites(.x, .coverage_seq)
      )
    )

    # now use coverage value to find sample size to reach target effect size
    # this one is tricky because sometimes we have to split by multiple powers,
    # and for method == "two", we have to split by multiple pilot_n_sites
    ss_to_reach_target <-
      # first, split by power group =====
    purrr::map(
      coverage_value_to_reach_target %>% split(f = .$power),

      function(covgroup) {

        # then, split by commuity and approximate with respective pilot_n_sites
        purrr::map(.x = switch(.method, "two" = pilot_n_sites,"single" = list(pilot_n_sites)),
                   .f = ~approx(c(0:40),
                                .x,
                                covgroup$coverage_rank) %>%
                     tibble::as_tibble() %>%
                     dplyr::rename(coverage_rank = x,
                                   n_samples = y) %>%
                     dplyr::mutate(n_samples = round(n_samples,0)))
        # then merge together by community
      } %>% dplyr::bind_rows(.id = "comm_group"))  %>%
      # then merge together by power
      dplyr::bind_rows(.id = "power")



    # merge all together
    ann <- coverage_rank_to_reach_target %>%
      dplyr::left_join(coverage_value_to_reach_target,
                       by = dplyr::join_by(power,coverage_rank)) %>%
      dplyr::left_join(ss_to_reach_target, by = dplyr::join_by(power,coverage_rank))

    # add group annotation if method is two
    if(.method == "single") {ann <- ann %>% dplyr::mutate(ann = n_samples)}
    if(.method == "two") {ann <- ann %>% dplyr::mutate(ann = paste(n_samples, comm_group))}


    ann_sum <- ann %>%
      dplyr::group_by(power, target_eff_size) %>%
      dplyr::summarize(target_eff_size = unique(target_eff_size),
                       coverage_rank = unique(coverage_rank),
                       coverage_value = unique(coverage_value),
                       n_samples = sum(n_samples),
                       ann = paste0(ann,collapse = "\n"),
                       .groups = "drop")

    # if cost is supplied, add cost to annotation
    if(!is.null(.cost_per_sample)) {

      ann_sum <- ann_sum %>%
        dplyr::mutate(total_cost = scales::dollar(n_samples * .cost_per_sample)) %>%
        dplyr::mutate(ann = paste0(ann, "\n(",total_cost,")"))


    } # end if cost

    return(ann_sum)
}


.add_target_to_plot <- function(.p3, .target_ann){

  .p3 +
    # add horizontal segment from regression line to minimum detectable effect
    ggplot2::geom_segment(data = .target_ann,
                          ggplot2::aes(x = -Inf,
                                       xend = coverage_rank,
                                       y = target_eff_size,
                                       yend = target_eff_size,
                                       color = power),
                          linetype = "dotdash") +

  #  ggplot2::geom_text(data = .target_ann,
  #                     ggplot2::aes(x = 4,
  #                                  y = target_eff_size,
  #                                  label = paste0("T = ",round(target_eff_size,2))),
  #                     hjust = 0,
  #                     vjust = -.5,
  #                     size = 3.5,
  #                     fontface = "bold",
  #                     show.legend = F) +
#
    # add vertical line at ss for target eff size
    ggplot2::geom_segment(data = .target_ann,
                          ggplot2::aes(x = coverage_rank,
                                       xend = coverage_rank,
                                       y = -Inf,
                                       yend = target_eff_size,
                                       color = power),
                          linetype = "dotdash") +
   # ggplot2::geom_text(data = .target_ann,
   #                    ggplot2::aes(x = coverage_rank,
   #                                 y = 0.03,
   #                                 label = paste0("T = ",ann),
   #                                 color = power),
   #                    show.legend = F,
   #                    hjust = -.1,
   #                    vjust = 1,
   #                    size = 3.5,
   #                    lineheight = .85,
   #                    fontface = "bold") +
    # add label in corner listing both values
    ggplot2::geom_label(data = .target_ann,
                       ggplot2::aes(x = coverage_rank,
                                    y = target_eff_size,
                                    label = paste0("T = ",round(target_eff_size,2),"\n",
                                                   ann),
                                    color = power),
                       hjust = 0,
                       vjust = 1,
                       size = 3.5,
                       lineheight = .85,
                       fontface = "bold",
                       alpha = .8,
                       show.legend = F)

}

