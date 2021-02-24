
# categorize_strategy -----------------------------------------------------
#' @title Categorize the strategy of brain areas based on the relationship between
#' median count of + cells and median of mean intensity of protein expression.
#'
#' @description  By using a lm or loess model as discriminant, the function outputs a dataframe
#' where each brain area ("my_grouping") is provided with a strategy score ("strategy"). The strategy
#' score is a number between 0 and 1, where 1 indicates that all samples of that brain area were categorized
#' in "count strategy", and 0 that all samples of that brain area were categorized as "intensity strategy"
#'
#' @param region_df region_based dataframe. Each row is a brain area ("my_grouping") per sample ("sample_id"), where
#' corrected cell count ("cells_perthousand") and average maximum intensity of the protein of interest ("intensity")
#' have been summarized. It can be output from summarize_per_region(), where the meta variable "group" has been added.
#' @param group_to_categorize refers to the group to categorize. All groups are used to build the discriminant.
#' @param lm_or_loess specify the type of model to build the discriminant. Can have values c("lm", "loess"). Defaults to lm.
#'
#' @return
#' @export
#'
#' @examples
#' x <- data.frame(
#' batch = rep(c(1,1,2,2), each = 5),
#' group = rep(c("control", "exp", "exp", "control"), each = 5),
#' sample_id = rep(c("a", "b", "c", "d"), each = 5),
#' my_grouping = rep(c("CA1", "CA2", "CA3", "DG", "BLA"), 4),
#' intensity = sample(10000, 20, replace = TRUE),
#' cells_perthousand = abs(rnorm(20))
#' )
#'
#' y <- categorize_strategy(x, "control")

categorize_strategy <- function(region_df, group_to_categorize, lm_or_loess = "lm") {

  # check that function inputs have correct class
  assertthat::assert_that(is.data.frame(region_df))
  # check that there are the correct variable names
  assertthat::has_name(region_df, "sample_id")
  assertthat::has_name(region_df, "my_grouping")
  assertthat::has_name(region_df, "intensity")
  assertthat::has_name(region_df, "cells_perthousand")
  assertthat::has_name(region_df, "group")

  assertthat::assert_that(all(group_to_categorize %in% region_df$group))

  assertthat::are_equal(length(lm_or_loess), 1)
  assertthat::are_equal(lm_or_loess %in% c("lm", "loess"), TRUE)

  # get prediction line
  line <- ggplot2::ggplot() +
    ggplot2::geom_smooth(data = region_df,
                ggplot2::aes(cells_perthousand, intensity), method = lm_or_loess)

  # from ggplot, I take the predicted values of loess
  gg_lm <- ggplot2::ggplot_build(line)$data[[1]][,c("x","y")]
  md <- loess(gg_lm$y ~ gg_lm$x)

  # categorize areas based on the fitted line
  strategy <- region_df %>%

    # select only timepoint 0
    dplyr::filter(group %in% group_to_categorize) %>%

    # get predicted values
    dplyr::ungroup() %>%
    dplyr::mutate(y_predicted = predict(md, cells_perthousand)) %>%
    dplyr::mutate(strategy = ifelse(intensity >= y_predicted,
                             "intensity_strategy", "count_strategy"))


  # calculate consistency categorizations across samples
  strategy_summ <-  strategy %>%

    # count per categorization
    dplyr::group_by(my_grouping, strategy) %>%
    dplyr::count()  %>%

    # change to percentages
    dplyr::ungroup() %>%
    tidyr::pivot_wider(names_from = strategy, values_from = n) %>%
    dplyr::mutate(dplyr::across(dplyr::everything(), ~tidyr::replace_na(.x, 0)),
           total = count_strategy+intensity_strategy,
           strategy = count_strategy/total) %>%
    dplyr::select(-c(total, count_strategy, intensity_strategy))

  return(strategy_summ)



}



# find_changed_strategy --------------------------------------------
#' @title Find brain areas that meet the criteria for analysis of strategy across groups
#'
#' @description Evaluation of the rate of change in count and intensity across brain areas. The criterion is set by
#' two conditions: first condition based on consistency (i.e. change in the same direction across samples), and the
#' second condition based on effect sizes (i.e. minimum effect size for the rate of change to be of interest). The function
#' outputs a dataframe with the comparison performed, the brain area ("my_grouping") and towards which strategy (more intensity
#'  or more count) the brain area met the criteria.
#'
#' @param region_df region_based dataframe. Each row is a brain area ("my_grouping") per sample ("sample_id"), where
#' corrected cell count ("cells_perthousand") and average maximum intensity of the protein of interest ("intensity")
#' have been summarized. These variables have been preprocessed to "cells_perthousand_box_scaled" and "intensity_box_scaled".
#' The samples belong to experimental groups ("group") and were processed in batches ("batch", i.e. one sample per group).
#' This data frame can be the output of preprocess_per_region().
#' @param consistency_n_samples number of samples which require to have consistent results in (at least one)
#' experimental group.
#' @param comparison_group group against which comparison are made.
#' @param min_change_rate minimum change in rate between comparison and experimental groups. Defaults to 1, i.e. any change.
#'
#' @return
#' @export
#'
#' @examples
#' # prepare dataframe
#' x <- data.frame(
#' batch = rep(c(1,1,2,2), each = 5),
#' group = rep(c("control", "exp", "exp", "control"), each = 5),
#' sample_id = rep(c("a", "b", "c", "d"), each = 5),
#' my_grouping = rep(c("CA1", "CA2", "CA3", "DG", "BLA"), 4),
#' intensity = sample(10000, 20, replace = TRUE),
#' cells_perthousand = abs(rnorm(20))
#' )
#' x$intensity_box_scaled <- scale(x$intensity)
#' x$cells_perthousand_box_scaled <- scale(x$cells_perthousand)
#'
#' y <- find_changed_strategy(x, 2, "control")
#' z <- find_changed_strategy(x, 2, "control", 2)
#'

find_changed_strategy <- function(region_df, consistency_n_samples, comparison_group, min_change_rate = 1) {

  # check that function inputs have correct class
  assertthat::assert_that(is.data.frame(region_df))

  # check that there are the correct variable names
  assertthat::has_name(region_df, "sample_id")
  assertthat::has_name(region_df, "batch")
  assertthat::has_name(region_df, "my_grouping")
  assertthat::has_name(region_df, "intensity")
  assertthat::has_name(region_df, "cells_perthousand")
  assertthat::has_name(region_df, "group")
  assertthat::has_name(region_df, "cells_perthousand_box_scaled")
  assertthat::has_name(region_df, "intensity_box_scaled")

  tot_samples <- length(unique(region_df$batch))
  assertthat::is.count(consistency_n_samples)
  assertthat::assert_that(consistency_n_samples <= tot_samples)

  assertthat::are_equal(length(comparison_group), 1)
  assertthat::assert_that(comparison_group %in% unique(region_df$group))

  assertthat::is.number(min_change_rate)
  assertthat::assert_that(min_change_rate != 0)


  rate_of_change <- region_df %>%

    # change comparison_group to my own level (required for later)
    dplyr::mutate(group = ifelse(group == comparison_group, "base", as.character(group))) %>%

    # get median of each ba for each timepoint
    dplyr::group_by(my_grouping, batch, group) %>%
    dplyr::summarize(count = median(cells_perthousand_box_scaled, na.rm = T),
                    intensity = median(intensity_box_scaled, na.rm = T)) %>%

    # wide format
    tidyr::pivot_wider(id_cols = c(my_grouping, batch), names_from = group,
                       values_from = c(count, intensity)) %>%

    # calculate difference from 0
    dplyr::mutate(
      dplyr::across(dplyr::starts_with("count_"), ~{.x - count_base}),
      dplyr::across(dplyr::starts_with("intensity_"), ~{.x - intensity_base})
    ) %>%

    # keep relevant vars
    dplyr::select(-c(count_base, intensity_base)) %>%
    dplyr::select(my_grouping, batch, starts_with("count_"), starts_with("intensity_")) %>%

    # to long format
    tidyr::pivot_longer(cols = -c(my_grouping, batch),
                        names_to = "comparison",
                        values_to = "difference") %>%

    # split comparison var
    tidyr::separate(comparison, into = c("type", "diff")) %>%

    # calculate ration of rate of change
    tidyr::pivot_wider(names_from = type, values_from = difference) %>%
    dplyr::mutate(rate = count / intensity) %>%

    # unique identified
    dplyr::mutate(ba_group = paste(my_grouping, diff, sep = "_"))


  # probability of change > change for count and intensity in same direction
  first_condition_met <- rate_of_change %>%

    # count how often condition met
    dplyr::group_by(ba_group) %>%
    dplyr::summarize(
      count_positive = sum(count > 0, na.rm = TRUE),
      intensity_positive = sum(intensity > 0, na.rm = TRUE)
    ) %>%
    dplyr::ungroup() %>%

    # filter based on cons9stency
    dplyr::mutate(
      extreme_count = ifelse(count_positive >= consistency_n_samples | count_positive <= (tot_samples - consistency_n_samples), TRUE, FALSE),
      extreme_intensity = ifelse(intensity_positive >= consistency_n_samples | intensity_positive <= (tot_samples - consistency_n_samples), TRUE, FALSE),
      extreme_both = ifelse(extreme_count == TRUE & extreme_intensity == TRUE, TRUE, FALSE)
    ) %>%

    # filter condition met
    dplyr::filter(extreme_both == TRUE) %>%

    # get brain areas
    dplyr::pull(ba_group) %>%
    unique()

  if(is.null(first_condition_met)) stop("No brain area met the consistency condition.")

  # keep only first condition met
  region_df_filtered <- region_df %>%
    dplyr::mutate(ba_group = paste(my_grouping, group, sep = "_")) %>%
    dplyr::filter(ba_group %in% first_condition_met) %>%
    droplevels()


  # calculate effect sizes
  effect_sizes <- data.frame()

  for (i in unique(region_df_filtered$ba_group)) {


    filtering_var <- stringr::str_split(i, pattern = "_")[[1]]


    dat <- region_df %>%
      dplyr::filter(my_grouping == filtering_var[[1]],
                    group %in% c(comparison_group, filtering_var[[2]])) %>%
      dplyr::select(my_grouping, batch, group, cells_perthousand, intensity) %>%
      dplyr::mutate(group = ifelse(group == comparison_group, "control", "exp")) %>%
      dplyr::mutate(rate = (cells_perthousand + 10) / (intensity + 10)) %>%
      tidyr::pivot_wider(names_from = group, values_from = c(cells_perthousand, intensity, rate))



    effect_sizes <- effect_sizes %>%
      dplyr::bind_rows(
        data.frame(
          ba_group = i,
          rate_hedges = effsize::cohen.d(dat$rate_exp,
                                         dat$rate_control,
                                         hedges.correction = TRUE)$estimate

        )
      )


  }


  res <- effect_sizes %>%

    # get goruping back
    tidyr::separate(ba_group, into = c("my_grouping", "comparison")) %>%

    # second_condition
    dplyr::filter(rate_hedges <= (1/min_change_rate) | rate_hedges >= min_change_rate) %>%
    dplyr::mutate(strategy_towards = ifelse(rate_hedges <= 1/min_change_rate, "intensity","count")) %>%

    # clean up
    dplyr::select(comparison, my_grouping, strategy_towards) %>%
    dplyr::mutate(comparison = paste(comparison, "vs", comparison_group))

  if(nrow(res) == 0) stop("No brain area met the effect size condition.")

  return(res)

}



