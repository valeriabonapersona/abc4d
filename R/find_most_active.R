
# find_most_active --------------------------------------------------
#' @title Find most active brain areas
#'
#' @description Wrapper around dplyr functions to find most active brain areas
#' by looking at the frequency of the top 5% (or other probability specified by the user) of the var distribution.
#' The function returns a dataframe with a summary per group of which brain areas were in the top distribution in how many batches.
#'
#' @param region_df region_based dataframe. Each row is a brain area ("my_grouping") per sample ("sample_id"), where
#' corrected cell count ("cells_perthousand") has been summarized. It contains a variable "batch" that identifies the unit where to
#' perform the calculation. If from a block design, "batch" identifies a unique set control and experimental groups (var "group"), with 1 sample each.
#' It can be output from summarize_per_region() or preprocess_per_region().
#' @param high_prob number between 0 and 1 indication the threshold for being a highly active region. 0.95 corresponds to the top 5%.
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
#' intensity_ave = sample(10000, 20, replace = TRUE),
#' cells_perthousand = abs(rnorm(20))
#' )
#'
#' find_most_active(x)


find_most_active <- function(region_df, high_prob = 0.95) {

  # check that function inputs have correct class
  assertthat::assert_that(is.data.frame(region_df))
  # check that there are the correct variable names
  assertthat::has_name(region_df, "sample_id")
  assertthat::has_name(region_df, "my_grouping")
  assertthat::has_name(region_df, "cells_perthousand")
  assertthat::has_name(region_df, "batch")
  assertthat::has_name(region_df, "group")


  # high_prob must be a number between 0 and 1
  assertthat::is.number(high_prob)
  assertthat::are_equal(all(high_prob >= 0 & high_prob <= 1), TRUE)

    cells_ba_probs <- region_df %>%

      # get upper quantile
      dplyr::group_by(batch) %>%
      dplyr::summarize(high_count = quantile(cells_perthousand, probs = high_prob)) %>%

      # clean up
      dplyr::ungroup() %>%

      # join with cells_ba
      dplyr::right_join(region_df, by = "batch")


    dat <- cells_ba_probs %>%

      # select ba above boundary
      dplyr::filter(cells_perthousand >= high_count) %>%

      # calculate blocks
      dplyr::select(batch, group, my_grouping) %>%
      unique() %>%

      dplyr::group_by(group, my_grouping) %>%
      dplyr::summarize(n_samples = length(unique(batch)))

    return(dat)

}




# prepare_sim_weights ---------------------------------------------------------
#' @title Prepare weights for simulation
#'
#' @description Wrapper around dplyr functions to to prepared weights for simulation study, by considering aba expression
#' levels as well as increase in protein expression due to the experimental manipulation. It returns a dataframe in long format with one
#' brain area "my_grouping" per group ("group") with the expression levels ("mean_expression" and "sd_expression") as well as the
#' median ratio of expression against the against_group ("weight") have been summarized.
#'
#' @param region_df region_based dataframe. Each row is a brain area ("my_grouping") per sample ("sample_id"), where
#' corrected cell count ("cells_perthousand") has been summarized. It contains a variable "batch" that identifies the unit where to
#' perform the calculation. If from a block design, "batch" identifies a unique set control and experimental groups (var "group"), with 1 sample each.
#' It can be output from summarize_per_region() or preprocess_per_region().
#' @param classifications dataframe with one brain area per row ("my_grouping") classified according to categorizations ("parents") found in Allen Brain Atlas
#' expression levels of the mRNA of interest. For an example of how to create this dataframe see X.
#' @param aba_api_summary dataframe where each row is a brain area "acronym" for which the average Allen Brain Atlas expression levels ("mean_expression") and deviation
#' ("sd_expression") have been summarized. It can be the output of from_aba_api_to_df
#' @param against_group group from "group" of region_df against which comparisons will be made
#' @return
#' @export
#'
#' @examples
#'   # create dataframes
#'   x <- data.frame(
#'   batch = rep(c(1,1,2,2), each = 5),
#'   group = rep(c("control", "exp", "exp", "control"), each = 5),
#'   sample_id = rep(c("a", "b", "c", "d"), each = 5),
#'   my_grouping = rep(c("CA1", "CA2", "CA3", "DG", "BLA"), 4),
#'   intensity_ave = sample(10000, 20, replace = TRUE),
#'   cells_perthousand = abs(rnorm(20))
#'   )
#'   y <- data.frame(
#'   my_grouping = c("CA1", "CA2", "CA3", "DG", "BLA"),
#'   parents = c(rep("hippocampus",4), "cortical subplate")
#'   )
#'   z <- data.frame(
#'   acronym = c("hippocampus", "cortical subplate"),
#'   mean_expression = rnorm(2, 10, 1),
#'   sd_expression = abs(rnorm(2))
#'   )
#'
#'   # run
#'   prepare_sim_weights(x,y,z, "control")
#'


prepare_sim_weights <- function(region_df, classifications, aba_api_summary, against_group) {

  # check that function inputs have correct class
  assertthat::assert_that(is.data.frame(region_df))
  # check that there are the correct variable names
  assertthat::has_name(region_df, "sample_id")
  assertthat::has_name(region_df, "my_grouping")
  assertthat::has_name(region_df, "cells_perthousand")
  assertthat::has_name(region_df, "batch")
  assertthat::has_name(region_df, "group")

  # classifications
  assertthat::assert_that(is.data.frame(classifications))
  # check that there are the correct variable names
  assertthat::has_name(classifications, "my_grouping")
  assertthat::has_name(classifications, "parents")

  # classifications
  assertthat::assert_that(is.data.frame(aba_api_summary))
  # check that there are the correct variable names
  assertthat::has_name(aba_api_summary, "acronym")
  assertthat::has_name(aba_api_summary, "mean_expression")
  assertthat::has_name(aba_api_summary, "sd_expression")

  # check against_group
  assertthat::are_equal(length(against_group),1)
  assertthat::are_equal(against_group %in% region_df$group, TRUE)

  # vars
  n_group <- length(unique(region_df$group))
  main_groups <- unique(region_df$group)

  # create dataset with weights
  weights <- classifications %>%

    # get cfos values
    dplyr::left_join(aba_api_summary %>%
                       dplyr::rename(parents = acronym), by = "parents") %>%

    # areas without expression levels get average of the rest
    dplyr::mutate(
      mean_expression = ifelse(is.na(mean_expression), mean(mean_expression, na.rm = TRUE), mean_expression),
      sd_expression = ifelse(is.na(sd_expression), 0, sd_expression)
    )


  # Weights through all groups
  group_weights <- region_df %>%

    # select relevant vars
    dplyr::select(batch, group, sample_id, my_grouping, cells_perthousand) %>%

    # get t0 values for ratios
    dplyr::left_join(region_df %>%
                       dplyr::filter(group == against_group) %>%
                       dplyr::select(batch, my_grouping, cells_perthousand) %>%
                       dplyr::rename(cells_perthousand_baseline = cells_perthousand),
              by = c("batch", "my_grouping")) %>%

    # calculate ratios
    dplyr::mutate(ratio = cells_perthousand / cells_perthousand_baseline) %>%

    # get median ratio per time point across blocks
    dplyr::group_by(group, my_grouping) %>%
    dplyr::summarize(weight = median(ratio, na.rm = TRUE))

  # merge expression and weights
  weights_sim <- do.call("rbind", replicate(n_group, weights, simplify = FALSE))

  weights_sim <- weights_sim %>%

    # get time var and associate different probs
    dplyr::mutate(
      group = rep(main_groups, each = nrow(weights))
    ) %>%

    # get time ratio probabilities
    dplyr::left_join(group_weights, by = c("group", "my_grouping"))

    return(weights_sim)

}


# from_aba_api_to_df ------------------------------------------------------
#' @title Convert data downloaded from Allen Brain Atlas API to a dataframe
#'
#' @description The data must have been previously downloaded from the Allen Brain Atlas API, and saved as .csv.
#' The function converts the data so that the final output is a dataframe with a row for each brain area ("acronym"),
#' its full name ("name") and the mean ("mean_expression") as well as standard deviation ("sd_expression") of the protein of
#' interest across all Allen Brain Atlas successful experiments.
#'
#' @param api_in_csv data downloaded from Allen Brain Atlas api in .csv format
#'
#' @return
#' @export
#'
#' @examples
#' FOR EXAMPLE SEE X.

from_aba_api_to_df <- function(api_in_csv) {

  # check that function inputs have correct class
  assertthat::assert_that(is.data.frame(api_in_csv))


  # find variables from api
  structure_var <-
    names(api_in_csv)[c(which(str_detect(names(api_in_csv), "structure.acronym")),
                      which(str_detect(names(api_in_csv), "structure.name")))]
  expression_var <- names(api_in_csv)[which(str_detect(names(api_in_csv), "expression.energy"))][1]

  # convert data
  cfos <- api_in_csv %>%

    # select relevant columns
    dplyr::select(all_of(c(structure_var, expression_var))) %>%

    # rename columns
    dplyr::rename(
      acronym = structure_var[1],
      name = structure_var[2],
      expression = expression_var[1]
    ) %>%

    # remove rows
    dplyr::drop_na()


  # summary stats per brain area
  cfos_summary_stats <- cfos %>%

    # per brain area
    dplyr::group_by(acronym, name) %>%
    dplyr::summarize(
      mean_expression = mean(expression, na.rm = TRUE),
      sd_expression = sd(expression, na.rm = TRUE)
    )

  return(cfos_summary_stats)


}



