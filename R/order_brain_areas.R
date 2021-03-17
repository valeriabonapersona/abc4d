
# find_intercept ----------------------------------------------------------
#' @title Function to find where an interpolated linear curve crosses a y value
#'
#' @description The function interpolates a linear model within two points, and finds at
#' which x the curve has a specified y.
#'
#' @param x number of x variable
#' @param y number of y variable
#' @param inter number of y at which you want to calculate x
#'
#' @return
#' @export
#'
#' @examples
#'find_intercept(x = c(-1.1, 1.1) , y = c(-1.1, 0.9))


find_intercept <- function(x, y, inter = 0) {

  assertthat::is.number(x)
  assertthat::is.number(y)
  assertthat::is.number(inter)

  approximate_function <- approxfun(x, y)
  approximate_function(inter)
}



# order_brain_areas -------------------------------------------------------
#' @title Order time of activation of brain areas based on their expression levels
#'
#' @description By using expression levels in different groups (e.g. time points),
#' brain areas are ordered .
#'
#' @param region_df region_based dataframe. Each row is a brain area ("my_grouping") per sample, with a variable for counts
#' (e.g. "cells_perthousand_box_scaled_ba"). The variable "group" specifies the experimental group, and it is used for the ordering.
#' It can be output from summarize_per_region() or preprocess_per_region().
#' @param order_of_groups string with the same unique elements as region_df$group, but ordered as should be used for the analysis.
#' @param count_var variable specified by the user on which to conduct the analysis.
#'
#' @return
#' @export
#'
#' @examples
#' # baby example
#' x <- data.frame(
#' batch = rep(c("1", "2"), 2),
#' group = rep(rep(c("0", "90"), each = 2), 2),
#' my_grouping = "CA1",
#' cells_perthousand_box_scaled_ba = c(-2, -1, 2, 1, -3, -2, 4, 2)
#' )
#' order_brain_areas(x, order_of_groups = c("0", "90"))



order_brain_areas <- function(region_df, order_of_groups, count_var = "cells_perthousand_box_scaled_ba") {

  # check region_df is a dataframe
  assertthat::assert_that(is.data.frame(region_df))

  # has correct variables
  assertthat::has_name(region_df, "my_grouping")
  assertthat::has_name(region_df, "group")

  # check groups are consistent
  assertthat::are_equal(length(order_of_groups) >= 2, TRUE)
  assertthat::are_equal(all(order_of_groups %in% region_df$group), TRUE)

  # count_var is 1
  # count_var is in names(region_df)
  assertthat::are_equal(length(count_var), 1)
  assertthat::has_name(region_df, count_var)


  # prepare df
  dat <- region_df %>%
    dplyr::select(group, my_grouping, !!count_var)
  names(dat) <- c("group", "my_grouping", "cells")

  # median df
  dat_summary <- dat %>%

    # get medians
    dplyr::group_by(group, my_grouping) %>%
    dplyr::summarize(cells_median = median(cells, na.rm = TRUE)) %>%
    dplyr::ungroup()

  # get intercept for all i+1 groups
  # gets saved in res
  i <- 1
  res <- NULL

  repeat {

    cross <- dat_summary %>%

      # keep only groups of interest for intercept
      dplyr::filter(group %in% order_of_groups[i:(i+1)]) %>%
      dplyr::mutate(group = as.numeric(as.character(group))) %>%

      # get intercept
      dplyr::group_by(my_grouping) %>%
      dplyr::summarize(inter = ifelse(any(cells_median >= 0), find_intercept(cells_median, group), NA)) %>%

      # clean_up
      dplyr::ungroup() %>%
      tidyr::drop_na()

    res <- rbind(res, cross)
    i <- i+1

    if (i == length(order_of_groups)) break
  } # repeat until no more found

  # look if there are duplicates
  duplicated_ba <- unique(res$my_grouping[which(duplicated(res$my_grouping))])

  if (length(duplicated_ba) != 0) {
    warning(paste(paste(duplicated_ba, collapse = ", "), "crossed the 0 more than one time."))
  }

  return(res)

}









