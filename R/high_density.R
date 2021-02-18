
# calculate_density --------------------------------------------------
#' @title Calculate density of active cells
#'
#' @description Wrapper around kde() function from the ks package to calculate the density of
#' a brain area across samples. It saves each sample separately in a specified folder.
#'
#' @param cells_df dataframe with xyz coordinates ("xPos", "yPos", "zPos") of the brain area of interest.
#' Each row corresponds to a different cell that belongs to a sample ("sample_id")
#' @param n_slice integer which specifies in how many parts to slice the estimate. Default = 10.
#' @param path specify path were data for each sample separtely will be saved.
#'
#' @return
#' @export
#'
#' @examples
#' x <- data.frame(
#' sample_id = rep(c(1:3), each = 100),
#' xPos = c(rnorm(100), rpois(100, 2), sample(100),
#' yPos = c(rnorm(100), rpois(100, 2), sample(100),
#' zPos = c(rnorm(100), rpois(100, 2), sample(100)
#' )



calculate_density <- function(cells_df, n_slice = 10, path) {

  # check that function inputs have correct class
  assertthat::assert_that(is.data.frame(cells_df))
  # check that there are the correct variable names
  assertthat::has_name(cells_df, "sample_id")
  assertthat::has_name(cells_df, "xPos")
  assertthat::has_name(cells_df, "yPos")
  assertthat::has_name(cells_df, "zPos")

  # check that n_slice is an integer
  assertthat::is.count(n_slice)

  # housekeep
  all_samples <- unique(cells_df$sample_id)
  my_prob <- seq(0, 1, l = n_slice + 1)

  #run
  x <- purrr::map(all_samples, function(x) {
    # filter data
    y <- cells_df %>%
      dplyr::filter(sample_id == x) %>%
      dplyr::select(xPos, yPos, zPos)

    # calculate density at each point
    k <-  ks::kde(x = y, eval.points = y)

    # divide based on quantiles
    k$quantiles <- quantile(k$estimate, probs = my_prob)
    k$slice <- c(1:n_slice)[cut(k$estimate, breaks = k$quantiles)]

    # create output
    k_res <- k$eval.points %>%

      # create df with info
      dplyr::mutate(
        sample_id = x,
        slice = k$slice
        )

    return(k_res)

  })

  # bind_together ?

  y <- dplyr::bind_rows(x) %>%
    # put sample_id as first column
    dplyr::relocate(sample_id)

  return(y)


}

