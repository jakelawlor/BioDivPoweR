

#' Plot bootstrapped community raw values
#'
#' make histograms of richness in rarefied bootstrap community 1 and community 2
#'
#' @param summary tibble of summary values from rarefied bootstraps. Output of `.summarize_boot_diff`
#' @return ggplot histogram of richness in bootstrapped, equal-coverage community pairs
#' @noRd
.plot_boot_values <- function(.boot_summary){

  p1 <- .boot_summary |>
    tidyr::pivot_longer(cols = contains(".rich"),
                        names_to = "community",
                        values_to = "richness",
                        names_pattern = "(.*).rich") %>%

    ggplot2::ggplot() +
    ggplot2::theme_bw() +
    ggplot2::geom_histogram(
      ggplot2::aes(x = richness,
                   fill = community,
                   color = community),
      alpha = .5,
      center = 0,
      linewidth = .1,
      binwidth = ceiling(
        (max(.boot_summary[2:3]) - min(.boot_summary[2:3]))
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

  return(p1)

}

#' Plot bootstrapped community raw differences
#'
#' make histogram of raw richness difference between bootstrapped community pairs
#'
#' @param summary tibble of summary values from rarefied bootstraps. Output of `.summarize_boot_diff`
#' @return ggplot histogram of richness difference between bootstrapped, equal-coverage community pairs
#' @noRd
.plot_boots_raw_diff <- function(.boot_summary){

  p2 <- .boot_summary %>%
    # get raw difference - or third col minus second col
    dplyr::mutate(raw_diff = dplyr::pull(.,3) - dplyr::pull(., 2)) %>%
    ggplot2::ggplot() +
    ggplot2::theme_bw() +
    ggplot2::geom_histogram(ggplot2::aes(x = raw_diff,
                                         fill = ifelse(log2_rich_diff == 0,"grey20","grey80")),
                            alpha = .5,
                            color = "grey20",
                            #fill = "grey80",
                            center = 0,
                            linewidth = .1,
                            binwidth = 1,
                            position = "identity") +
    ggplot2::labs(x = "Raw richness difference",
                  y = "Count") +
    ggplot2::scale_fill_identity()


}


#' Plot bootstrapped community log2 differences
#'
#' make histogram of log2 richness difference between bootstrapped community pairs
#'
#' @param summary tibble of summary values from rarefied bootstraps. Output of `.summarize_boot_diff`
#' @return ggplot histogram of log2-fold richness difference between bootstrapped, equal-coverage community pairs
#' @noRd
.plot_boots_log2_diff <- function(.boot_summary, min_exp_n, n_eff_size_bins){

  # first calculate bindiwth
  xrange <- range(.boot_summary$log2_rich_diff, na.rm = TRUE)
  binwidth <- (xrange[2] - xrange[1]) / n_eff_size_bins


  p3 <- ggplot2::ggplot(data = .boot_summary[.boot_summary$log2_rich_diff != 0,],
                        ggplot2::aes(x = log2_rich_diff)) +
    ggplot2::theme_bw() +
    ggplot2::geom_histogram(
      boundary = 0,
      color = "grey20",
      fill = "grey90",
      linewidth = .1,
      binwidth = binwidth,
      show.legend = F
    ) +
    ggplot2::geom_hline(yintercept = min_exp_n,
                        linetype = "dashed",
                        linewidth = .4) +
    ggplot2::labs(x = "Log2 Fold-Difference in Richness",
                  y = "Count")

  # color the histogram bins in which count surpasses min_exp_n
  p3 <- p3 +
    # add colored rectangles
    ggplot2::geom_rect(data = p3 %>% ggplot2::layer_data() %>%
                         dplyr::filter(count >= min_exp_n) %>%
                         dplyr::mutate(eff_size_num = (round(x,3)))# %>%
                        # dplyr::filter(round(x,3) != 0)
                       ,
                       inherit.aes = F,
                       ggplot2::aes(xmin = xmin,
                                    xmax = xmax,
                                    ymin = 0,
                                    ymax = min_exp_n,
                                    fill = abs(eff_size_num)),
                       color = "grey20",
                       linewidth = .1, show.legend = F) +

   ## add grey rectangle
   #ggplot2::geom_rect(data = p3 %>% ggplot2::layer_data() %>%
   #                     dplyr::filter(round(x,3) == 0)
   #                   ,
   #                   inherit.aes = F,
   #                   ggplot2::aes(xmin = xmin,
   #                                xmax = xmax,
   #                                ymin = 0,
   #                                ymax = min_exp_n),
   #                   fill = "grey50",
   #                   color = "grey20",
   #                   linewidth = .1, show.legend = F) +

    ggplot2::scale_fill_gradientn(colors = cool_matlab())

  return(p3)

}


