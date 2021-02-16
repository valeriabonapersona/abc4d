
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
#' levels as well as increase in protein expression due to the experimental manipulation.
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


prepare_aba_weights <- function(region_df, high_prob = 0.95) {

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



