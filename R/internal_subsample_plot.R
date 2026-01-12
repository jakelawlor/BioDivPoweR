# Helper functions for plotting within subsample_boots function

#' Plot subsample funnels
#'
#' Make funnel plots of detected richness differences in bootstrapped communities, subsampled at a gradient of coverage values.
#'
#' @param df dataframe of richness differences between boots at all values in `coverage_seq`. Output of `.summarize_subsamples_diff()`
#' @return ggplot of effect size as coverage increases for `min_exp_n` experiments in all qualifying effect size bins
#' @noRd
.plot_subsample_funnel <- function(.df, .coverage_seq){

  # plot funnels within effect size bins
  p0 <- .df %>%
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
      breaks = seq(1,41,by=4),
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
    ggplot2::scale_color_gradientn(colors = cool_matlab()) +
    ggplot2::scale_x_continuous(
      breaks = seq(1,41,by=2),
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
    ggplot2::theme(panel.grid = ggplot2::element_blank()) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed")

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
    ggplot2::scale_color_gradientn(colors = cool_matlab(),
    ) +
    ggplot2::scale_x_continuous(
      breaks = unique(.prop_correct$coverage_rank)[c(T,F)],
      labels = round(unique(.coverage_seq*100))[c(T,F)]
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



#' Plot coverage by detectable effect - base plot only
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
  #labels <- .min_detectable$pilot_n_sites[.min_detectable$power == min(.min_detectable$power)]

  # plot minimum detectable effect size across coverage
  p3 <- .min_detectable %>%
    ggplot2::ggplot(ggplot2::aes(x = coverage_rank,
                                 y = abs_eff_size,
                                 group = power)) +
    ggplot2::geom_point(data = . %>% dplyr::filter(include %in% c(T,F)),
                        size = 2, alpha = .8,
                        ggplot2::aes(fill = abs_eff_size,
                                     color = power,
                                     shape = include),
                        #shape = 21,
                        stroke = 1) +
    ggplot2::scale_shape_manual(values = c(`FALSE` = 4, `TRUE` = 21)) +
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
      xlim = c(1,40),
      clip = "off"
    ) +
    ggplot2::scale_y_reverse(expand = ggplot2::expansion(mult = c(.1, 0))) +
    ggplot2::scale_fill_gradientn(colors = cool_matlab(),
                                  limits = c(0,max(abs(.min_detectable$abs_eff_size)))) +
    ggplot2::scale_color_grey(end = .55,
                              start = 0) +
    ggplot2::scale_x_continuous(
      breaks = .min_detectable$coverage_rank[1:40][c(T,F)],
      labels = round(unique(.coverage_seq)[1:40]*100)[c(T,F)],
      sec.axis = ggplot2::sec_axis(~.,
                                   breaks = .min_detectable$coverage_rank[1:40][c(T,F)],
                                   labels = .min_detectable$sec_axis_labels[1:40][c(T,F)],
                                   name = "Sample n in pilot community")

    ) +
    ggplot2::theme_bw() +
    ggplot2::labs(y = "Minimum Detectable Fold Change in Richness",
                  x = "Coverage (%)",
                  fill = "Absolute\nEffect Size",
                  color = "Power") +
    ggplot2::theme(panel.grid = ggplot2::element_blank()) +
    ggplot2::geom_hline(yintercept = 0,
                        linewidth = .2) +
    #ggplot2::theme(axis.text.x.top = ggplot2::element_text(hjust = c(.95,rep(.5, times = 19)))) +
    ggplot2::guides(shape = ggplot2::guide_none())

  return(p3)

}



#' Find coverage in pilot community
#'
#' find coverage in the pilot community at 1-k samples.
#'
#' @param .pilot pilot dataset
#' @param .method analysis method: "single" or "two" treatment.
#' @param .coverage_seq coverage intervals to test. imported internally.
#' @param .community_types labels for community types. imported internally.
#' @return character vector of axis 41 axis labels, corresponding to number of sites in empirical community (or communities) to reach coverage intervals.
#' @noRd
.find_pilot_coverage <- function(.pilot, .method, .coverage_seq ){


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

  # if the highest value isn't enough, do it again for 5000
  if(.method == "single")
    if(max((pilot_coverage)) < .coverage_seq[40]){
      pilot_coverage <- find_coverage(pilot_occupancy, 5000)}
  if(.method == "two")
    if(any(purrr::map(pilot_coverage, max) < .coverage_seq[40])){
      pilot_coverage <- purrr::map(.x = pilot_occupancy, .f = ~find_coverage(.x, 5000))}


  return(pilot_coverage)


}


#' Create labels for the second X-axis
#'
#' Use pilot coverage to identify the number of samples in pilot communit(y/ies)
#' that it would take to reach each coverage interval.
#' If applicable, add cost.
#'
#' @param .pilot_coverage pilot coverage at 1-1000 samples
#' @param .coverage_seq sequence of coverage values in our nonlinear scale
#' @param .cost_per_sample cost per unit sample, if provided
#' @param .method single or double analysis pipeline.
#' @return character vector of axis 41 axis labels, corresponding to number of sites in empirical community (or communities) to reach coverage intervals.
#' @noRd
.calc_second_axis <- function(.pilot_coverage, .coverage_seq, .cost_per_sample, .method, .community_types){

  # e. find number of samples in the empirical community (with coverage rescaled)
  # to reach the coverage values in coverage_seq.
  pilot_n_sites <- switch(
    .method,
    "single" = find_eff_sites(.pilot_coverage, .coverage_seq[1:40]),
    "two" = purrr::map(
      .x = .pilot_coverage,
      .f = ~c(find_eff_sites(.x, .coverage_seq[1:40]))
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
        sec_axis_labels

      }
    )

  }

  sec_axis_df <- data.frame(coverage_rank = c(1:40),
                            coverage = .coverage_seq[1:40])

  sec_axis_df <- switch(
    .method,

    # if single treatment:
    "single" = sec_axis_df %>%
      dplyr::mutate(
        pilot_n_sites = pilot_n_sites,
        sec_axis_labels = sec_axis_labels
      ) %>%
      dplyr::mutate(
        include = ifelse(pilot_n_sites < 6,F,T)
      ),

    # if two treatment:
    "two" = sec_axis_df %>%
      dplyr::mutate(
        !!paste0("pilot_n_sites.",.community_types[[1]]) := pilot_n_sites[[paste0(.community_types[[1]])]],
        !!paste0("pilot_n_sites.",.community_types[[2]]) := pilot_n_sites[[paste0(.community_types[[2]])]],
        sec_axis_labels = sec_axis_labels
      ) %>%
      dplyr::mutate(
        include = ifelse(

            get(paste0("pilot_n_sites.",.community_types[[1]])) < 6 |
            get(paste0("pilot_n_sites.",.community_types[[2]])) < 6,
          F,T)
      )

  ) # end switch



  return(sec_axis_df)

}


#' Calculate achieved detectable effect size at pilot sample size
#'
#' Use linear interpolations to backcalculate the achieved minimum detectable
#' effect from the regression that we established for detectable effect as a
#' function of sample coverage
#'
#' @param .pilot species-by-site matrix from pilot study
#' @param .pilot_coverage community coverage of pilot at 1-1000 samples
#' @param .coverage_seq sequence of coverage values in our nonlinear scale
#' @param .power_mod functional relationship between coverage and detectable effect
#' @param .method single or double analysis pipeline.
#' @return tibble of achieved minimum effect size at pilot sample size for each value power
#' @noRd
.calc_minEff_at_sampleSize <- function(.pilot, .pilot_coverage, .coverage_seq, .power_mod, .method){

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
    "single" = .pilot_coverage[ss],
    "two" = {
      vals <- purrr::map2(
        .x = .pilot_coverage,
        .y = .pilot %>% dplyr::count(across(2)) %>% tibble::deframe() %>% as.list,
        .f = ~.x[.y] )
      vals[which.min(vals)]}
  )


  # 2. convert coverage to linear scale
  # coverage from the pilot is between 0 and 1, but we need it to be
  # between 0 and 40, which is our "ranked" linear coverage scale.
  # use linear approximation from the points 0-40 and the values in coverage_seq.
  pilot_coverage_rank <- approx(x= .coverage_seq,
                                y=c(1:41),
                                xout = ss_coverage)$y

  # step 2. predict achieved minimum detectable effect size
  # at the coverage rank (1-41) from the pilot study
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
  pilot_minimum_detectable$pilot_achieved_cover <- unlist(ss_coverage)
  pilot_minimum_detectable$pilot_achieved_rank <- pilot_coverage_rank

  return(pilot_minimum_detectable)

}


#' Annotate achieved minimum detectable effect
#'
#' Add lines to base plot of achieved minimum effect size at the pilot sample size at each power
#'
#' @param .p3 plot of coverage-to-detectable effect relationship
#' @param .pilot_minimum_detectable tibble of coverage value and detectable effect size at each power
#' @return p3 with lines and annotations for "achieved" detectable effect
#' @noRd
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





#' Calculate sites for target power
#'
#' Approximate number of samples needed to achieve target minimum detectable effect size
#'
#' @param .target_eff_size effect size that user wants to reliably detect
#' @param .power_au augmented coverage-to-detectable-effect model
#' @param .pilot_coverage community coverage of pilot community at 1-1000 samples
#' @param .coverage_seq vector of coverage values represented on the x axis (custom scale)
#' @param .cost_per_sample if applicable, cost per sample
#' @param .method single- or two-treatment analysis type
#' @return tibble of number of sites needed to reach target detectable effect size
#' @noRd
.calc_site_to_reach_target <- function(.target_eff_size = target_eff_size,
                                       .power_au = power_au,
                                       .pilot_coverage = pilot_coverage,
                                       .coverage_seq = coverage_seq,
                                       .cost_per_sample = cost_per_sample,
                                       .method = method){


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
      .f = ~approx(c(1:41),
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
      "single" = find_eff_sites(.pilot_coverage, .coverage_seq[1:40]),
      "two" = purrr::map(
        .x = .pilot_coverage,
        .f = ~find_eff_sites(.x, .coverage_seq[1:40])
      )
    )

    # now use coverage value to find sample size to reach target effect size
    # this one is tricky because sometimes we have to split by multiple powers,
    # and for method == "two", we have to split by multiple pilot_n_sites
    ss_to_reach_target <-
      # first, split by power group - - - - - - -
    purrr::map(
      coverage_value_to_reach_target %>% split(f = .$power),

      function(powergroup) {

        # then, split by commuity and approximate with respective pilot_n_sites
        purrr::map(.x = switch(.method, "two" = pilot_n_sites,"single" = list(pilot_n_sites)),
                   .f = ~approx(c(1:40),
                                .x,
                                powergroup$coverage_rank) %>%
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
    if(.method == "single") {ann <- ann %>% dplyr::mutate(ann = as.character(n_samples))}
    if(.method == "two") {ann <- ann %>% dplyr::mutate(ann = paste(n_samples, comm_group))}

    # if any of the target values are NA, that probably means the requested power is too low

    if(any(is.na(ann[,c("coverage_rank", "coverage_value")]))){
      pow <- ann$power[is.na(ann$coverage_rank)]
      eff <- ann$target_eff_size[is.na(ann$coverage_rank)]
      lab <- paste0("effect size ", eff," at power ",pow)
      lab <- paste0(lab, collapse = ", ")
      warning(paste0("Requested effect size and/or power are outside augmented model range. They will be labeled on plot and returned as sample size zero.",
                     "(",lab,")")
              , call. = F)
      ann <- ann %>% tidyr::replace_na(
        list(coverage_rank = 0,
             coverage_value  = 0,
             n_samples = 0,
             ann = "outside range"))

      }


    ann_sum <- ann %>%
      tidyr::pivot_wider(
        names_from = comm_group,
        values_from = n_samples,
        names_prefix = "sample_size."
      ) %>%
      dplyr::group_by(power, target_eff_size) %>%
      dplyr::summarize(
        coverage_rank  = unique(coverage_rank),
        coverage_value = unique(coverage_value),
        sample_size.total  = sum(dplyr::c_across(starts_with("sample_size.")), na.rm = TRUE),
        ann            = paste0(ann, collapse = "\n"),
        dplyr::across(starts_with("sample_size."), function(x) max(x, na.rm=T)),   # keep the new columns
        .groups = "drop"
      )

    # remove the n_samles.community 1 if single-treatment analysis (redundant)
    if(.method == "single"){
      ann_sum$sample_size.1 <- NULL
      ann_sum$sample_size.2 <- NULL
    }


    # if cost is supplied, add cost to annotation
    if(!is.null(.cost_per_sample)) {
      ann_sum <- ann_sum %>%
        dplyr::mutate(total_cost = scales::dollar(sample_size.total * .cost_per_sample)) %>%
        dplyr::mutate(ann = paste0(ann, "\n",total_cost))
    } # end if cost

    return(ann_sum)
}


#' Annotate samples to reach target on base plot
#'
#' Add lines to base plot of number of samples needed to reach target detectable effect size
#'
#' @param .p3 plot of coverage-to-detectable effect relationship
#' @param .target_ann tibble of coverage values needed to detect target effect size
#' @return p3 with lines and annotations for "target" detectable effect
#' @noRd
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

    # add vertical line at ss for target eff size
    ggplot2::geom_segment(data = .target_ann,
                          ggplot2::aes(x = coverage_rank,
                                       xend = coverage_rank,
                                       y = -Inf,
                                       yend = target_eff_size,
                                       color = power),
                          linetype = "dotdash") +

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
                       alpha = 1,
                       show.legend = F)

}



#' Calculate true effect
#'
#' For two-treatment comparison, calculate the true difference in richness between the two input communities, which will be equal to richness in the undersampled community - rarefied richness in 50 bootstrapped iterations of the oversampled community
#'
#' @param .pilot species-by-site matrix of pilot sample
#' @param .pilot_coverage vector of coverage values at 1-1000 samples in pilot community
#' @return tibble of values for pilot community
#' @noRd
.calc_true_effect <- function(.pilot, .pilot_coverage){

  # for two-treatment, first find n_sites for equal coverage
  pilot_split <- .pilot %>% split(f = .[2])

  raw_richness <- purrr::map(
    .x = pilot_split,
    .f = ~sum(colSums(.x[,-c(1:2)]) > 0)
  ) %>% dplyr::bind_rows() %>%
    purrr::set_names(paste0("raw_richness.",names(.)))

  raw_sample_size <- purrr::map(.x = pilot_split,
                                .f = ~nrow(.x)) %>%
    dplyr::bind_rows() %>%
    purrr::set_names(paste0("raw_sample_size.",names(.)))


  # find the minimum coverage of the two communities at their original samples size
  min_cov <- min(.pilot_coverage[[1]][raw_sample_size[[1]]],
                 .pilot_coverage[[2]][raw_sample_size[[2]]])

  # find the number of sites in original communities at which coverage is equal
  # to the smaller of the two max coverages
  pilot_n_equal_cov <- c(
    dplyr::first(which(.pilot_coverage[[1]] >= min_cov)),
    dplyr::first(which(.pilot_coverage[[2]] >= min_cov))
  ) %>%
    purrr::set_names(names(.pilot_coverage))

  # identify the community that's equal-coverage n_samples is the same as its total n_samples.
  # that one is the relatively udersampled community, for which richness will be measured once
  # the community whose n_equal_coverage is not the same as its full n is the relatively
  # oversampled community, for which richness will be measured by bootstrapping 50 times and
  # finding mean richness in each bootstrap.
  equal <- purrr::map2(
    lapply(pilot_split, nrow),
    pilot_n_equal_cov,
    ~.x == .y
  )

  # sometimes, neither community is undersampled,
  # if they both reach the same target coverage in the same number of sites.
  bothequal <- all(purrr::map(.x = equal,.f = ~isTRUE(.x)) %>% unlist())


  if(bothequal == F){
    undersampled <- names(equal[equal == T]) %>% purrr::set_names(.)
    oversampled <- names(equal[equal != T]) %>% purrr::set_names(.)
  } else if (bothequal){
    undersampled <- names(equal[1]) %>% purrr::set_names(.)
    oversampled <- names(equal[2]) %>% purrr::set_names(.)
  }

  # undersampled_rich
  # (unedited richness in the non-rarefied community)
  # Note that this only happens if the communities are not both equally-sampled.
  # if they are equally sampled, we'll take rarefied richness of both of them
  rich_under <- switch(
    paste(bothequal),
    "FALSE" = sum(colSums(pilot_split[[undersampled]][,-c(1:2)],)>0) %>% purrr::set_names(undersampled),
    "TRUE" = purrr::map(
      .x = c(1:50),
      .f = ~sum(
        colSums(pilot_split[[undersampled]][sample(nrow(pilot_split[[undersampled]]),
                                                  pilot_n_equal_cov[[undersampled]],
                                                  replace = F),][,-c(1:2)]) > 0)
    ) %>% unlist()
  )


  # oversampled_rich
  # rarefied richness in 50 iterations of the rarefied oversampled community
  rich_rare_over <- purrr::map(
    .x = c(1:50),
    .f = ~sum(
      colSums(pilot_split[[oversampled]][sample(nrow(pilot_split[[oversampled]]),
                                                pilot_n_equal_cov[[oversampled]],
                                                replace = F),][,-c(1:2)]) > 0)
  ) %>% unlist()

  rich_rare_over_mean <- rich_rare_over %>% mean() %>%  purrr::set_names(oversampled)
  rich_rare_over_sd <- rich_rare_over %>% sd() %>%  purrr::set_names(oversampled)

  richness <- c( rich_rare_over_mean, rich_under)[names(.pilot_coverage)]

  true_eff <- log2(richness[1]/richness[2])

  richness %>% purrr::set_names(paste0(names(.),".rarefied_rich")) %>% c()

  richness.all <-
    raw_richness %>%
    dplyr::bind_cols(raw_sample_size) %>%
    dplyr::bind_cols(min_cov %>% purrr::set_names("min_coverage") %>% dplyr::bind_rows()) %>%
    dplyr::bind_cols(richness %>% purrr::set_names(paste0("rarefied_richness.",names(.))) %>% dplyr::bind_rows()) %>%
    dplyr::bind_cols(pilot_n_equal_cov %>% purrr::set_names(paste0("rarefied_sample_size.",names(.))) %>% dplyr::bind_rows()) %>%
    dplyr::bind_cols(true_eff %>% purrr::set_names("true_abs_effect") %>% dplyr::bind_rows() %>% abs())


  return(richness.all)

}



#' find pilot values; single
#'
#' Calculate raw sample size, raw community coverage, and raw richness of pilot sample
#'
#' @param .pilot species-by-site matrix of pilot sample
#' @param .pilot_coverage vector of coverage values at 1-1000 samples in pilot community
#' @return tibble of values for pilot community
#' @noRd
.calc_raw_vals_single <- function(.pilot, .pilot_coverage){


  sample_size <- nrow(.pilot)

  coverage <- .pilot_coverage[sample_size]

  richness <- sum(colSums(.pilot[,-1]) > 0)

  richness.all <- data.frame(
    raw_sample_size = sample_size,
    coverage = coverage,
    raw_richness = richness
  )

  return(richness.all)

}


#' Add true effect to plot
#'
#' Draw a line on the y axis of p3 to display the true effect from pilot sample (for two-treatment analysis only)
#'
#' @param .true_effect tibble of true richness difference between two community types
#' @param .p3 ggplot of coverage-to-detectable-effect relationship
#' @return p3 with a line marking "true" effect on the y axis
#' @noRd
.add_true_effect <- function(.true_effect, .p3){


  true_eff <- .true_effect$true_abs_effect


  p3 <- .p3 +
    ggplot2::annotate(geom = "rug",
                      y = true_eff,
                      linewidth = 1,
                      length = ggplot2::unit(6,"pt"),
                      color = "black") +
    ggplot2::annotate(geom = "text",
                      y = true_eff,
                      x = 0,
                      label = "True\ndifference",
                      hjust = 0,
                      lineheight = .85,
                      size = 3.5)

  return(p3)

}
