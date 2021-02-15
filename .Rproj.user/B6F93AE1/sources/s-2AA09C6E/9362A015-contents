
# adapt_estimation_atlas --------------------------------------------------
#' @title Adapt estimation atlas from Erö et al. to be used in the analysis
#'
#' @description Wrapper around dplyr functions to select only relevant
#' brain areas of the cell estimation atlas. The atlas estimating the number of cells per brain area was generated
#' in this publication (doi:10.3389/fninf.2018.00084). The atlas follows the Allen Brain
#' Reference Atlas Categorization (mouse brain).
#'
#' @param estimation_atlas atlas from Erö et al. Use atlas as present in the package, or
#' provide a dataframe where each row is a brain area. The dataframe must contain a variable
#' called "Regions" with the names on the brain areas. The other variable(s) are the estimations,
#' i.e "Neurons", "Glia", "Inhibitory" etc.
#' @param adj_aba_atlas dataframe with Allen Brain Atlas tree, with an additional variable
#' called "my_grouping" with the level of categorization of interest. The dataframe contains also
#' the variable "name" specifying the name of the brain areas. For an example of how to create
#' this dataframe, see X.
#'
#' @return
#' @export
#'
#' @examples
#' x <- data.frame(
#' Regions = c("a_1", "a_2", "a_3", "b_1", "b_2"),
#' Cells = c(10000, 2100, 39847, 754, 923)
#' )
#'
#' y <- data.frame(
#' name = c("a_1", "a_2", "a_3", "b_1", "b_2", "c_1", "c_2"),
#' my_grouping = c(rep("a",3), rep("b",2), rep("c", 2))
#' )
#'
#' adapt_estimation_atlas(x,y)


## upload data to package
adapt_estimation_atlas <- function(estimation_atlas, adj_aba_atlas) {

  # check that estimation atlas and adj aba atlas are dataframe
  assertthat::assert_that(is.data.frame(estimation_atlas))
  assertthat::assert_that(is.data.frame(adj_aba_atlas))

  # check that there are the correct variable names
  assertthat::has_name(estimation_atlas, "Regions")
  assertthat::has_name(adj_aba_atlas, "my_grouping")
  assertthat::has_name(adj_aba_atlas, "name")


  # check that Regions and name have the same levels
  assertthat::are_equal(sum(!estimation_atlas$Regions %in% adj_aba_atlas$name, na.rm = TRUE), 0)

  # function
  estimation_atlas %>%

    # rename vars for consistency
    dplyr::rename_all(tolower) %>%
    dplyr::rename(name = regions) %>%

    # merge atlases
    dplyr::full_join(adj_aba_atlas[,c("name", "my_grouping")], by = "name") %>%

    # summarize per grouping
    dplyr::group_by(my_grouping) %>%
    dplyr::summarise(dplyr::across(-name,
                                   ~ sum(.x, na.rm = TRUE))) %>%

    # remove NA of non grouped areas
    dplyr::filter(!is.na(my_grouping))

}


# summarize_per_region --------------------------------------------------
#' @title Transform xyz coordinates data in a summary per brain region
#'
#' @description Wrapper around dplyr functions to create a summary dataframe
#' with information about count and intensity per brain area. In the output:
#' count_perthousand refers to the number of cells expressing the protein of
#' interest per thousand of total cells in that brain area, while intensity_ave
#' refers to the average intensity of active cells. You can also specify the type
#' of cells estimated.
#'
#' @param xyz_coordinates dataframe where each row is a cell. It contains the following variables:
#' 'sample_id' to describe the sample; 'my_grouping' for the brain areas; 'maxInt' for the maximum intensity
#' of the protein per identified cell.
#'
#' @param estimation_atlas estimation atlas already adjusted for the brain areas of interest. Can
#' be output from adapt_estimation_atlas().
#' @param meta datafame with meta-information about samples. It expects a variable named "sample_id".
#' @param cells_type type of cells from estimation atlas to be used for count correction. Look at the variables
#' of estimation_atlas for types. If not specified, uses all types.
#'
#' @return
#' @export
#'
#' @examples
#'
#' x <- data.frame(
#' sample_id = c(rep("a", 100), rep("b", 75), rep("c", 50)),
#' my_grouping = c(rep(c("CA1","CA2","CA3", "DG", "BLA"), each = 20),
#' rep(c("CA1", "CA2", "CA3"), each = 20), rep("BLA", 15), rep(c("DG", "BLA"), each = 25)),
#' maxInt = abs(rnorm(250, 100))
#' )
#'
#' y <- data.frame(
#' my_grouping = c("CA1","CA2","CA3", "DG", "BLA"),
#' cells = sample(100000, 5),
#' glia = sample(10000, 5)
#' )
#'
#' z <- data.frame(
#' sample_id = c("a", "b", "c"),
#' batch = c(1,1,2),
#' group = c("A", "B", "A")
#' )
#'
#'summarize_per_region(x,y,z)


## upload data to package
summarize_per_region <- function(xyz_coordinates, estimation_atlas, meta, cells_type = "cells") {

  # check that function inputs have correct class
  assertthat::assert_that(is.data.frame(xyz_coordinates))
  assertthat::assert_that(is.data.frame(estimation_atlas))
  assertthat::assert_that(is.data.frame(meta))

  # check that there are the correct variable names
  assertthat::has_name(xyz_coordinates, "sample_id")
  assertthat::has_name(xyz_coordinates, "my_grouping")
  assertthat::has_name(xyz_coordinates, "maxInt")

  assertthat::has_name(estimation_atlas, "my_grouping")

  assertthat::has_name(meta, "sample_id")

  # check that Regions and name have the same levels
  if(cells_type != "cells"){
    assertthat::assert_that(length(cells_type) == 1)
    assertthat::assert_that(is.character(cells_type))
  }

  assertthat::has_name(estimation_atlas, cells_type)

  # function
  xyz_coordinates %>%

    # summarize per brain area - hemispheres together
    dplyr::group_by(sample_id, my_grouping) %>%
    dplyr::summarize(
      count = length(maxInt),
      intensity = sum(maxInt)
    ) %>%

    # correct by total cells per brain areas
    dplyr::left_join(estimation_atlas[,c("my_grouping", cells_type)], by = "my_grouping") %>%

    # normalize
    dplyr::mutate(cells_perthousand = count/cells*1000,
           intensity_ave = intensity/count
    ) %>%
    dplyr::rename(atlas_cells = cells) %>%

    # clean up tibble
    dplyr::ungroup() %>%
    droplevels() %>%

    # merge with meta information
    left_join(meta, by = "sample_id")

}


# preprocess_per_region --------------------------------------------------
#' @title Normalize and scale for region-based analyses
#'
#' @description Wrapper around dplyr functions to normalize and standardize
#' the count and intensity variables. It expects the output of summarize_per_region.
#'
#' @param region_df region_based dataframe. Each row is a brain area ("my_grouping") per sample ("sample_id"), where
#' corrected cell count ("cells_perthousand") and average maximum intensity of the protein of interest ("intensity_ave") have been summarized.
#' It can be output from summarize_per_region(). The data will be normalized according to "batch",
#' and it will be scaled per unit ("cells_perthousand_box_scaled", "intensity_ave_box_scaled"),
#' as well as per brain area per unit ("cells_perthousand_box_scaled_ba", "intensity_ave_box_scaled_ba")
#' @param
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
#' preprocess_per_region(x)


preprocess_per_region <- function(region_df) {

  # check that function inputs have correct class
  assertthat::assert_that(is.data.frame(region_df))
  # check that there are the correct variable names
  assertthat::has_name(region_df, "sample_id")
  assertthat::has_name(region_df, "my_grouping")
  assertthat::has_name(region_df, "intensity_ave")
  assertthat::has_name(region_df, "cells_perthousand")
  assertthat::has_name(region_df, "batch")

  region_df %>%

    # box cox transformation
    dplyr::group_by(batch) %>%
    dplyr::mutate(
      cells_perthousand_box = normalize(cells_perthousand),
      intensity_ave_box = normalize(intensity_ave)) %>%

    # normalization by block
    dplyr::mutate(
      cells_perthousand_box_scaled = scale(cells_perthousand_box),
      intensity_ave_box_scaled = scale(intensity_ave_box)

    ) %>%

    dplyr::ungroup() %>%

    # normalization by block by brain area
    dplyr::group_by(batch, my_grouping) %>%

    dplyr::mutate(
      cells_perthousand_box_scaled_ba = scale(cells_perthousand_box),
      intensity_ave_box_scaled_ba = scale(intensity_ave_box)
    ) %>%

    # clean up
    dplyr::ungroup() %>%
    droplevels()


}

