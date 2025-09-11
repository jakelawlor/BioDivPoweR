subsample_two_trt <- function(data,
                              pilot,
                              power = c(80),
                              target_eff_size = NULL,
                              cost_per_sample = NULL,
                              seed = 1){


  # set random seed
  set.seed(seed)

  cat("Subsampling bootstrapped communities...\n")

  # create a vector of the two types of communities
  community_types <- stringr::str_remove(
    colnames(data)[stringr::str_detect(colnames(data),"matrix")],
    "_matrix") %>% purrr::set_names(.)


  # rescale coverage ----------------------------------------------
  # in boots, we have the true coverage at every sample in each site.
  # Here, we'll assume the max sample number results in 100% coverage.

  # set waypoints across the coverage scale (x-axis)
  coverage_seq <- c(1,seq(from = 250, to = 10000, by = 250) %>% sqrt())/100

  # first, find the number of sites that correspond to each of the target coverages.
  data2 <- data %>%
    dplyr::mutate(
      !!(paste0("n_to_check_",community_types[[1]])) := purrr::map(
        .x = get(paste0(community_types[[1]],"_coverage")),
        # here, we'll rescale each of these coverage_seq values so
        # that instead of coverage_seq going 0-1, it will go 0-max(coverage)
        .f = ~BioDivSampler:::find_eff_sites(coverage_vec = .x,
                                             target = coverage_seq*max(.x)
                                             # note that the target here, we're multiplying the coverage sequence
                                             # 10, 15, 22, ... 100, by the max achieved coverage (~93% instead of 100),
                                             # so we're rescaling the coverage waypoints to the maximum coverage in the site
        )
      ),
      !!(paste0("n_to_check_",community_types[[2]])) := purrr::map(
        .x = get(paste0(community_types[[2]],"_coverage")),
        # here, we'll rescale each of these coverage_seq values so
        # that instead of coverage_seq going 0-1, it will go 0-max(coverage)
        .f = ~BioDivSampler:::find_eff_sites(coverage_vec = .x,
                                             target = coverage_seq*max(.x))
      ),
    )
  # this adds two columns to the boots dataset:
  # n_to_check_{comm_name_1}, and n_to_check_{comm_name_2}, which list the number
  # of sites needed to reach the coverage values across the waypoints
  # that we created for coverage scaling (1, 15, 22, etc..)

  rm(data)


  # sample richness at each n_to_check
  data3 <- data2 %>%
    dplyr::mutate(
      !!(paste0(community_types[[1]],"_rich")) := purrr::map2(
        .x = get(paste0(community_types[[1]],"_matrix")),
        .y = get(paste0("n_to_check_",community_types[[1]])),
        .f = ~ {
          # explicitly name the arguments for clarity
          comm_matrix <- .x
          sample_sizes <- .y

          purrr::map_dbl(sample_sizes, ~ {
            n <- .x # current sample size
            # downsample the matrix to 'n' rows - without replacement
            sampled_matrix <- comm_matrix[sample(nrow(comm_matrix), n, replace = FALSE), ]

            # calculate richness as the number of columns with colSum > 0
            # note that if n is 1, colSums does work, so we're doing it this way.
            if(n == 1){
              sum(sampled_matrix > 0)
            } else {
              sum(colSums(sampled_matrix) > 0)
            }
          })
        }),
      !!(paste0(community_types[[2]],"_rich")) := purrr::map2(
        .x = get(paste0(community_types[[2]],"_matrix")),
        .y = get(paste0("n_to_check_",community_types[[2]])),
        .f = ~ {
          # explicitly name the arguments for clarity
          comm_matrix <- .x
          sample_sizes <- .y

          purrr::map_dbl(sample_sizes, ~ {
            n <- .x #  current sample size
            # downsample the matrix to 'n' rows
            sampled_matrix <- comm_matrix[sample(nrow(comm_matrix), n,, replace = FALSE), ]

            # Calculate richness as the number of columns with colSum > 0
            if(n == 1){
              sum(sampled_matrix > 0)
            } else {
              sum(colSums(sampled_matrix) > 0)
            }
          })
        })
    )
  rm(data2)
  # this adds two columns to the tibble:
  # {comm1}_rich, and {comm2}_rich, which list the richness observed at samples
  # of n_to_check sites of each community.

  # unlist ------------------------------------------------------------------
  cat("Calculating richness differences between subsamples...\n")

  # find differences in richness between community 1 and community 2
  # at sample sizes that represent each coverage step.
  df <- data3 %>%
    # keep only effect size, and richness vectors at the coverage waypoints
    dplyr::select(eff_size, dplyr::contains("_rich")) %>%
    dplyr::group_by(eff_size) %>%
    # add trial (which will be the number of repeats within each effect size bin)
    dplyr::mutate(trial = c(1:dplyr::n())) %>%
    # unnest list columns
    tidyr::unnest(dplyr::contains("_rich")) %>%
    dplyr::group_by(eff_size,trial) %>%
    dplyr::mutate(coverage = coverage_seq) %>%
    dplyr::mutate(coverage_rank = 0:(dplyr::n()-1)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(log2diff = log2(
      get(paste0(community_types[[1]],"_rich")) / get(paste0(community_types[[2]],"_rich"))
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
  rm(data3)



  # plot funnels within effect size bins
  p0 <-   df %>%
    dplyr::group_by(eff_size, coverage_rank,eff_size_num) %>%
    dplyr::summarize(quant10 = quantile(log2diff2, .05),
                     quant50 = quantile(log2diff2, .5),
                     quant90 = quantile(log2diff2, .95),
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
      labels = round(c(coverage_seq*100),0)[c(T,F,F,F)]
    ) +
    ggplot2::labs(x = "Coverage (%)",
                  y = "Detected Richness Change",
                  title = "Convergence across coverage within effect size bins") +
    ggplot2::theme_bw()

  # plot mean richness change across coverage within bins
  p1 <- df %>%
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
      labels = round(c(coverage_seq*100),0)[c(T,F)]
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

  # find the width of the smallest bin
  zero_bin_width_0 <- df %>%
    dplyr::distinct(eff_size_num) %>%
    dplyr::arrange(abs(eff_size_num)) %>%
    dplyr::slice(2) %>% dplyr::pull()
  zero_bin_width <- abs(zero_bin_width_0)/2
  rm(zero_bin_width_0)


  # count correct detections ------------------------------------------------
  cat("Calculating proportion correct detection across subsamples...\n")

  # plot proportion correct sign, note that for the smallest bin which includes
  # zero, correct sign will be calculated instead as log2 richness being any value
  # inside the bin (from -small change to +small change)
  prop_correct <-  df %>%
    # first calculate correct sign, except in zero effect size bin
    dplyr::mutate(correct_sign =
                    dplyr::case_when(
                      # when eff_size is not zero, calculate whether detected change is same
                      # sign as true change
                      eff_size_num != 0 ~ ifelse(sign(log2diff2) == sign(eff_size_num),
                                                 T,F),
                      # when eff size is zero, calculate whether log2 change is within
                      # the limits of the zero effect size bin
                      eff_size_num == 0 ~ ifelse(abs(log2diff2) == 0,
                                                 T,F))) %>%
    dplyr::mutate(abs_eff_size = abs(eff_size_num)) %>%
    dplyr::group_by(abs_eff_size, coverage_rank, coverage) %>%
    # group_by(eff_size, eff_size_num,coverage,coverage_rank) %>%
    dplyr::summarize(n_correct  = sum(correct_sign),
                     n_total = dplyr::n(),
                     .groups = "drop") %>%
    dplyr::mutate(prop_correct = (n_correct/n_total)*100) %>%
    dplyr::ungroup()

  p2 <- prop_correct %>%
    ggplot2::ggplot(ggplot2::aes(x = coverage_rank,
                                 y = prop_correct,
                                 group = abs_eff_size,
                                 color = abs_eff_size)) +
    ggplot2::geom_line(linewidth = 1) +
    ggplot2::scale_color_gradientn(colors = BioDivSampler:::cool_matlab(),
    ) +
    ggplot2::scale_x_continuous(
      breaks = seq(0,40,by=2),
      labels = round(coverage_seq*100)[c(T,F)]
    ) +
    ggplot2::scale_y_continuous(labels = scales::percent_format(scale = 1)) +
    ggplot2::theme_bw() +
    ggplot2::labs(x = "Coverage (%)",
                  y = "Proportion Correct Direction",
                  color = "Absolute\nEffect Size") +
    ggplot2::theme(legend.position = "right") +
    ggplot2::theme(legend.key.width = ggplot2::unit(12,"pt"),
                   legend.key.height = ggplot2::unit(40,"pt"),
                   panel.grid = ggplot2::element_blank()) +
    ggplot2::guides(color = ggplot2::guide_colorbar(theme = ggplot2::theme(legend.title.position = "top"))) +
    ggplot2::geom_hline(data = data.frame(power = power),
                        ggplot2::aes(yintercept = power),
                        linetype = "dashed")




  # find minimum detectable effect size -------------------------------------

  cat("Modeling detectable richness as a function of sample coverage...\n")
  # plot the minimum effect size that can be detected with a given power
  # at each step of sample coverage along coverage_seq
  min_detectable <- purrr::map(
    .x = power %>% purrr::set_names(power),
    .f = ~ prop_correct %>%
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
      # remove the coverage = 100 point
      dplyr::filter(coverage != 1)
  )


  # model
  power_mod <- purrr::map(
    .x = min_detectable,
    .f = ~lm(data = .x,
             formula = abs_eff_size ~ coverage_rank)
  )

  # augment the model to get a regression line across the coverage_seq
  power_au <- purrr::map2(
    .x = power_mod,
    .y = min_detectable,
    .f = ~broom::augment(.x,
                         newdata = .y %>% tidyr::complete(
                           data.frame(coverage_rank = c(0:40))
                         ),
                         interval = "confidence")
  ) %>%
    dplyr::bind_rows(.id = "power")

  # bind detectable effect size list into one
  min_detectable <- min_detectable %>%
    dplyr::bind_rows(.id = "power")

  # plot minimum detectable effect size across coverage
  p3 <- min_detectable %>%
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

    ggplot2::geom_line(data = power_au,
                       ggplot2::aes(x = coverage_rank,
                                    y = .fitted,
                                    group = power,
                                    color = power)) +
    ggplot2::coord_cartesian(
      ylim = c(
        # higher limit is maximum of either observed effect size or auugmented
        # effect size from model ouput
        pmax(max(abs(power_au$.fitted)),max(abs(min_detectable$abs_eff_size))),
        # minimum is -.05
        0),
      xlim = c(0,40),
      clip = "off"
    ) +
    ggplot2::scale_y_reverse(expand = ggplot2::expansion(mult = c(.1, 0))) +
    ggplot2::scale_fill_gradientn(colors = cool_matlab(),
                                  limits = c(0,max(abs(min_detectable$abs_eff_size)))) +
    ggplot2::scale_color_grey(end = .55,
                              start = 0) +
    ggplot2::scale_x_continuous(
      breaks = seq(0,40,by=2),
      labels = round(coverage_seq*100)[c(T,F)]
    ) +
    ggplot2::theme_bw() +
    ggplot2::labs(y = "Minimum Detectable Fold Change in Richness",
                  x = "Coverage (%)",
                  fill = "Absolute\nEffect Size",
                  color = "Power") +
    ggplot2::theme(panel.grid = ggplot2::element_blank()) +
    ggplot2::geom_hline(yintercept = 0,
                        linewidth = .2)



  # calculate second axis - number of empirical samples ---------------------
  # we want to relate the coverage values along the bottom of the axis to the
  # actual number of samples needed in the empirical community. Note that this
  # will differ from the method for single treatment, as we have to get the
  # number of samples separately for each of the two treatments

  # a. find average species occupancy in pilot survey
  pilot_occupancy <- purrr::map(
    .x = pilot %>% split(f = pilot[[2]]),
    .f = ~colMeans(.x[,-c(1:2)])[colSums(.x[,-c(1:2)]) > 0],
    .progress = T
  )
  # here, we have a 2-element list (one for each community type) that includes
  # the average occupancy of each species within the bootstrapped samples,
  # filtered to species that have occupancy > 0, so they will be <= max_n_spp.


  # b. find coverage at 1:1000 samples in pilot survey
  pilot_coverage <-  purrr::map(
    .x = pilot_occupancy,
    .f = ~BioDivSampler:::find_coverage(occupancy = .x)
  )

  # c. find target detection rate,
  # or the coverage needed to detect all species in the community
  # find target detection rate,
  # or the coverage needed to detect all species in the community
  pilot_target_detection_rate <- purrr::map(
    .x =  pilot %>% split(f = pilot[[2]]),
    .f = ~find_target_detection_rate(.x)
  )

  # d. rescale the pilot coverage from 0-target detection rate, instead of 0-100
  pilot_coverage_rescale <- purrr::map2(
    .x = pilot_coverage,
    .y = pilot_target_detection_rate,
    .f = ~.x/.y
  )
  rm(pilot_occupancy, pilot_coverage, pilot_target_detection_rate)

  # e. find number of sites to reach each coverage_seq() value on the rescaled coverage scale
  pilot_n_sites <- purrr::map(
    .x = pilot_coverage_rescale,
    .f = ~BioDivSampler:::find_eff_sites(.x, coverage_seq)
  )
  # now we have the number of sites in the empirical community that would
  # allow us to reach target coverage of coverage_seq()
  sec_axis_labels <- pilot_n_sites

  # add dollar signs, and original labels underneath
  sec_axis_labels <- c(
    paste0(community_types[[1]],":\n",community_types[[2]],":"),
    paste0(pilot_n_sites[[1]],"\n", pilot_n_sites[[2]]))

  # add cost, if applicable, to second axis labels
  if(!is.null(cost_per_sample)){

    # calculate cost
    cost <- ((pilot_n_sites[[1]] + pilot_n_sites[[2]]) * cost_per_sample)

    # add dollar signs, and original labels underneath
    sec_axis_labels[2:42] <- paste0(scales::dollar(round(cost,0)),"\n", sec_axis_labels[2:42])

    # make only one of every 4 costs show.
    sec_axis_labels[2:41][c(T,F,T,T)] <- gsub("^.*?\\n","",sec_axis_labels[2:41][c(T,F,T,T)])

  }


  # add second axis scale to p3
  suppressMessages(
    suppressWarnings(
      p3 <- p3 +
        ggplot2::scale_x_continuous(
          breaks = seq(0,40,by=3),
          labels = round(coverage_seq*100)[c(T,F,F)],
          sec.axis = ggplot2::sec_axis(~.,
                                       breaks = seq(-1,40,by=2),
                                       labels = c(sec_axis_labels[c(T,F)]),
                                       name = "Sample n in pilot community")
        ) +
        ggplot2::theme(axis.text.x.top = ggplot2::element_text(hjust = c(1,rep(.5, times = 20))))
    )
  )


  # add achieved coverage values --------------------------------------------
  # we know how many samples were taken in the pilot community, so we'll add
  # the achieved coverage and achieved minimum detectable effect size at the
  # given power.

  # 0. because the power analysis only applies to the sample treatment with lower
  #    coverage, we need to first pick only the treatment with the lower of the two coverages.
  min_pilot_coverage <- purrr::map2(
    .x = pilot_coverage_rescale,
    .y = pilot %>% dplyr::count(veg_score) %>% tibble::deframe() %>% as.list,
    .f = ~.x[.y]
  )

  low_group <- names(which.min(unlist(min_pilot_coverage)))
  low_coverage_nrow <- nrow(pilot[pilot[[2]] == low_group,])

  # 1. approximate the coverage rank (0-40) from the coverage achieved
  #   in the pilot study, rescaled to 0-100, at the number of samples in the pilot
  pilot_approx_coverage_rank <- purrr::map(
    .x = min_pilot_coverage,
    .f = ~approx(x = coverage_seq,
                 y = c(0:40),
                 xout = .x)$y
  )

  # approximate achieved minium detectable effect size at pilot coverage
  # now that we know the coverage rank (0-40) of the pilot study, we can
  # predict the minimum detectable effect size from the model
  pilot_approx_effect_size <- purrr::map(
    .x = power_mod,
    .f = ~purrr::map(
      pilot_approx_coverage_rank,
      function(blah) predict(.x,
                             list(coverage_rank = blah))

  ) %>% dplyr::bind_rows(.id = "group"))%>%
    dplyr::bind_rows(.id = "power") %>%
    dplyr::rename("achieved_min_eff_size" = `1`)

  # add to plot
  p3 <- p3 +
    # add line segment at true number of samples
    ggplot2::annotate(geom = "segment",
                      x = pilot_approx_coverage_rank[[low_group]],
                      xend = pilot_approx_coverage_rank[[low_group]],
                      y = Inf,
                      yend = -Inf,
                      linetype = "dashed") +
    # add text at true number of samples
    ggplot2::annotate(geom = "text",
                      x = pilot_approx_coverage_rank[[low_group]],
                      y = 0.01,
                      hjust = -.05,
                      vjust = 1,
                      label = paste0("A = ",low_coverage_nrow," (",low_group,")"),
                      size = 3.5,
                      fontface = "bold") +
    # add horizontal segment to minimum detectable effect size
    ggplot2::geom_segment(inherit.aes = F,
                          data = pilot_approx_effect_size[pilot_approx_effect_size$group == low_group,],
                          ggplot2::aes(x = -Inf,
                                       xend = pilot_approx_coverage_rank[[low_group]],
                                       y = achieved_min_eff_size,
                                       yend = achieved_min_eff_size,
                                       color = power),
                          linetype = "dashed") +
    # add y-side annotation
    ggplot2::geom_text(
      data = pilot_approx_effect_size[pilot_approx_effect_size$group == low_group,],
      ggplot2::aes(
        x = 0,
        y = achieved_min_eff_size,
        label = paste0("A = ", round(achieved_min_eff_size,2)),
        color = power
      ),
      hjust = 0,
      vjust = -.5,
      size = 3.5,
      fontface = "bold",
      show.legend = F)



  # add samples to reach target ---------------------------------------------
  # IF the user supplied a target effect size, calculate the number of samples
  # needed to reach the target
  if(!is.null(target_eff_size)){

    coverage_rank_to_reach_target <- purrr::map(
      .x = power_au %>% split(f = .$power),
      # interpolate the relationship between model fit and coverage rank
      # to find the coverage rank where the model fit would be target_eff_size
      .f = ~approx(.x$.fitted,
                   .x$coverage_rank,
                   target_eff_size) %>%
        tibble::as_tibble() %>%
        dplyr::rename(target_eff_size = x,
                      coverage_rank = y)
    ) %>%
      dplyr::bind_rows(.id = "power")


    # now translate that coverage rank (0-40) to a coverage value (0-100)
    coverage_value_to_reach_target <- purrr::map(
      .x = power_au %>% split(f = .$power),
      .f = ~approx(.x$.fitted,
                   coverage_seq,
                   target_eff_size) %>%
        tibble::as_tibble() %>%
        dplyr::rename(target_eff_size = x,
                      coverage_value = y)
    ) %>%
      dplyr::bind_rows(.id = "power")


    #translate that coverage rank to a number of samples in the pilot community
    sample_n_to_reach_target <- purrr::map(
      .x = pilot_n_sites,
      .f = ~approx(unique(power_au$coverage_rank),
                   .x,
                   coverage_rank_to_reach_target$coverage_rank) %>%
        tibble::as_tibble() %>%
        dplyr::rename(coverage_rank = x,
                      n_samples = y) %>%
        dplyr::mutate(n_samples = round(n_samples,0))
    ) %>%
      dplyr::bind_rows(.id = "group")




    ann <- coverage_rank_to_reach_target %>%
      dplyr::left_join(coverage_value_to_reach_target,
                       by = dplyr::join_by(power, target_eff_size)) %>%
      dplyr::left_join(sample_n_to_reach_target,
                       by = dplyr::join_by(coverage_rank)) %>%
      dplyr::mutate(ann = paste(n_samples, group)) %>%
      dplyr::group_by(power, target_eff_size)

    ann_sum <- ann %>%
      dplyr::summarize(target_eff_size = unique(target_eff_size),
                coverage_rank = unique(coverage_rank),
                n_samples = sum(n_samples),
                ann = paste0(ann,collapse = "\n"))

    # if cost is supplied, add cost to annotation
    if(!is.null(cost_per_sample)) {

      ann_sum <- ann_sum %>%
        dplyr::mutate(total_cost = scales::dollar(n_samples * cost_per_sample)) %>%
        dplyr::mutate(ann = paste0(ann, "\n(",total_cost,")"))


    } # end if cost


    p3 <- p3 +
      ggplot2::geom_segment(data = ann_sum,
                            ggplot2::aes(x = coverage_rank,
                                         xend = coverage_rank,
                                         y = -Inf,
                                         yend = target_eff_size,
                                         color = power),
                            linetype = "twodash") +
      # add horizontal segment from regression line to minimum detectable effect
      ggplot2::geom_segment(data = ann_sum,
                            ggplot2::aes(x = -Inf,
                                         xend = coverage_rank,
                                         y = target_eff_size,
                                         yend = target_eff_size,
                                         color = power),
                            linetype = "twodash") +
      # add annotation at number of samples
      ggplot2::geom_text(data = ann_sum,
                         ggplot2::aes(x = coverage_rank,
                                      y = 0.01,
                                      label = paste0("T = ",ann),
                                      color = power),
                         show.legend = F,
                         hjust = -.1,
                         vjust = 1,
                         size = 3.5,
                         nudge_x = -.5,
                         lineheight = .85,
                         fontface = "bold") +
      ggplot2::geom_text(data = ann_sum,
                         ggplot2::aes(x = 0,
                                      y = target_eff_size,
                                      label = paste0("T = ",round(target_eff_size,2))),
                         hjust = 0,
                         vjust = -.5,
                         size = 3.5,
                         fontface = "bold",
                         show.legend = F)

  }


  # # save answers as df ------------------------------------------------------
  out.df <- pilot_approx_effect_size %>%
    dplyr::left_join(
      pilot %>% dplyr::group_by(dplyr::across(2)) %>% dplyr::count(name = "n_samples") %>% dplyr::rename(group = colnames(pilot)[2]),
      by = dplyr::join_by("group")
    ) %>%
    dplyr::rename(treatment = group,
                  min_detectable_eff_size = achieved_min_eff_size) %>%
    dplyr::relocate(treatment ) %>%
    dplyr::mutate(group = "achieved")

  if(!is.null(target_eff_size)){
    out.df <- out.df %>%
      dplyr::bind_rows(
        ann %>%
          dplyr::rename(treatment = group) %>%
          dplyr::rename(min_detectable_eff_size = target_eff_size) %>%
          dplyr::mutate(group = "target") %>%
          dplyr::select(power,
                        treatment,
                        min_detectable_eff_size,
                        group,
                        n_samples)
      ) %>%
      dplyr::relocate(treatment, group) %>%
      dplyr::arrange(group, power)
  }

  # add cost if applicable
  if(!is.null(cost_per_sample)){
    out.df <- out.df %>%
      dplyr::group_by(treatment) %>%
      dplyr::mutate(total_trt_samples = sum(n_samples)) %>%
      dplyr::mutate(sample_cost = total_trt_samples*cost_per_sample)
  }



  out <- list(p0, p1, p2, p3, out.df)
  return(out)

}
