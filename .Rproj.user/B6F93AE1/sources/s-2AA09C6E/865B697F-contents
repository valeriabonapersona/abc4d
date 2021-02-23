
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
#'
#'
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


