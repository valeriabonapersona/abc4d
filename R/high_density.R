
# calculate_density --------------------------------------------------
#' @title Calculate density of active cells
#'
#' @description Wrapper around kde() function from the ks package to calculate the density of
#' a brain area across samples. It saves each sample separately in a specified folder.
#'
#' @param cells_df dataframe with xyz coordinates ("xPos", "yPos", "zPos") of the brain area of interest.
#' Each row corresponds to a different cell that belongs to a sample ("sample_id")
#' @param n_slice integer which specifies in how many parts to slice the estimate. Default = 10.
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
#' y <- calculate_density(x)



calculate_density <- function(cells_df, n_slice = 10) {

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





# plotting_coordinates_2d ----------------------------------------------------
#' @title Plot combinations of xyz coordinates in 2D
#'
#' @description Wrapper around ggplot2 to plot xyz coordinates. Used later in plotting_all_coordinates()
#'
#' @param cells_df dataframe with xyz coordinates ("xPos", "yPos", "zPos") of the brain area of interest.
#' Each row corresponds to a different cell that belongs to a sample ("sample_id")
#'
#' @return


plotting_coordinates_2d <- function(data_main, var_to_colour, cols,
                                    x_coord, y_coord) {

  p <- ggplot2::ggplot() +
    ggplot2::geom_jitter(data = data_main,
                         ggplot2::aes_string(x_coord, y_coord, colour = var_to_colour)) +

    ggplot2::scale_colour_manual(values = cols) +
    ggplot2::theme(legend.position = "bottom")

  return(p)

}


# plotting_all_coordinates_2d --------------------------------------------------
#' @title Plot all combinations of xyz coordinates in 2D
#'
#' @description Wrapper around ggplot2 to plot xyz coordinates. Plots xy, yz, and xz coordinates and
#' merges the different plots together. It returns a list with 3 elements, each of which is a ggplot object.
#' These can be merged with other functions, such as ggpubr
#'
#' @param cells_df dataframe with xyz coordinates ("xPos", "yPos", "zPos") of the brain area of interest.
#' Each row corresponds to a different cell.
#' @param var_to_colour variable to use for colouring xyz coordinates. Default to NULL
#' @param cols colours to use. Must have length equal to the levels of var_to_colour. Default to NULL
#'
#' @return
#' @export
#'
#' @examples
#' x <- data.frame(
#' sample_id = rep(c(1:3), each = 100),
#' xPos = c(rnorm(100), rpois(100, 2), sample(100),
#' yPos = c(rnorm(100), rpois(100, 2), sample(100),
#' zPos = c(rnorm(100), rpois(100, 2), sample(100),
#' group = rep(c("a", "b", "c"), each = 100)
#' )
#' my_cols <- c("black", "green", "red")
#' g <- plotting_all_coordinates_2d(x, "group", my_cols)
#' ggpubr::ggarrange(g[[1]], g[[2]], g[[3]], nrow = 3, common.legend = TRUE, legend = "top")


plotting_all_coordinates_2d <- function(cells_df, var_to_colour = NULL, cols = NULL) {

  # check cells_df is a dataframe with xPos, yPoz, zPos
  assertthat::assert_that(is.data.frame(cells_df))
  assertthat::has_name(cells_df, "xPos")
  assertthat::has_name(cells_df, "yPos")
  assertthat::has_name(cells_df, "zPos")

  # if var_to_colour is not NULL, var to colour must be in cells_df.
  if (!is.null(var_to_colour)) {
    assertthat::assert_that(length(var_to_colour) == 1)
    assertthat::has_name(cells_df, var_to_colour)

    if (!is.null(cols)) {
      assertthat::are_equal(length(cols), length(unique(cells_df[,var_to_colour])))
    }
  }

  xy <- plotting_coordinates_2d(cells_df, var_to_colour, cols, x_coord = "xPos", y_coord = "yPos")
  xz <- plotting_coordinates_2d(cells_df, var_to_colour, cols, x_coord = "xPos", y_coord = "zPos")
  zy <- plotting_coordinates_2d(cells_df, var_to_colour, cols, x_coord = "yPos", y_coord = "zPos")


  g <- list(xy, xz, zy)
  return(g)

}




# voxelize ------------------------------------------------------------
#' @title Voxelize xyz coordinates
#'
#' @description Converts xyz coordinates in pixels to voxels of a specified µm size.
#' The function summarizes the number of samples per group. It outputs a list where the first element ("data")
#' is the data voxelized, and the second element ("summary") is the data summarized, i.e. the number of samples
#' per voxel per group.
#'
#'
#' @param data dataframe with xyz coordinates ("xPos", "yPox", "zPos") for example output of calculate_density().
#' It includes a variable ("sample_id") that is the quantity that will be summarizes over the group variable ("group").
#' @param conv_factor number to be used for the conversion. The conversion factor depends on the imaging settings.
#' For more information, see for example: https://mindresearchfacility.nl/wp-content/uploads/2019/12/Zoom-factor-and-corresponding-NA-and-resolution-by-pixel-size.pdf
#' Defaults to 5.16
#' @param voxel_size in µm. Defaults to 30. If you want to specify in pixels, put conv_factor = 1
#'
#' @return
#' @export
#'
#' @examples
#' x <- data.frame(
#' xPos = rnorm(100, 300, 100), yPos = rnorm(100, 300, 100), zPos = rnorm(100, 300, 100),
#' group = rep(c("control", "exp"), each = 50),
#' sample_id = rep(c(1:10), each = 10)
#' )
#' y <- voxelize(x, voxel_size = 100)

voxelize <- function(data, conv_factor = 5.16, voxel_size = 30) {

  # data is a dataframe with the required variables
  assertthat::assert_that(is.data.frame(data))
  # check that there are the correct variable names
  assertthat::has_name(data, "sample_id")
  assertthat::has_name(data, "group")
  assertthat::has_name(data, "xPos")
  assertthat::has_name(data, "yPos")
  assertthat::has_name(data, "zPos")

  # conv_factor is a number
  assertthat::is.number(conv_factor)
  assertthat::is.number(voxel_size)

  # size conv_factor and voxel_size is 1
  assertthat::are_equal(length(conv_factor), 1)
  assertthat::are_equal(length(voxel_size), 1)
  assertthat::are_equal(voxel_size != 0, FALSE)


  # get conv coefficients
  µm <- sapply(c("xPos", "yPos", "zPos"), function(x) {(max(data[,x]) - min(data[,x])) * conv_factor})

  # cut size
  µm <- round(µm / voxel_size)

  if(all(µm > 1) == FALSE) {warning("Voxel size requested is smaller than resolution.")}

  vox <- data %>%

    # voxelize xyq coordinates
    dplyr::mutate(
      xPos_vox = as.numeric(cut(xPos, µm[[1]])),
      yPos_vox = as.numeric(cut(yPos, µm[[2]])),
      zPos_vox = as.numeric(cut(zPos, µm[[3]]))
    )

  # summarizes n samples
  vox_summ <- vox %>%

    # summarize how many cells per sample in each voxel
    dplyr::group_by(group, xPos_vox, yPos_vox, zPos_vox) %>%
    dplyr::summarize(n_samples = length(unique(sample_id)))

  res <- list(
    data = vox,
    summary = vox_summ
  )

  return(res)


}

