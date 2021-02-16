# pairwise_activation --------------------------------------------------
#' @title Calculate pairwise activation across groups
#'
#' @description
#'
#' @param region_df region_based dataframe. It can be the output of preprocess_per_region(). Each row corresponds
#' to a brain area ("my_grouping") per sample, with (corrected) counts ("cells_perthousand_box_scaled"). The variable "group" refers to the experimental groups. The variable
#' "batch" is used for pairwise comparisons. Of note, for each "batch" there can be only a sample per group.
#' @param against_group group from "group" var against which comparisons will be made
#' @param count_var string with the name of the variable on which to perform the analyses
#' @param heatmap if TRUE, the function returns a list where the first element is the dataframe, and the second
#' element is a visualization with heatmap.
#'
#' @return
#' @export
#'
#' @examples
#'
#' dat <- data.frame(
#' group = rep(c("0", "30", "90"), each = 50),
#' my_grouping = rep(c(1:50), 3),
#' my_counts = rnorm(150)
#' )

#' res <- pairwise_activation(dat)
#' res_with_visualization <- pairwise_activation(dat, heatmap = TRUE, classifications_df = cat, classifications_cat = "parents")
#'
#'


pairwise_activation <- function(region_df, against_group = "0", heatmap = FALSE) {

  # check that function inputs have correct class
  assertthat::assert_that(is.data.frame(region_df))

  # check that against_group and count_var are one element
  assertthat::are_equal(length(against_group),1)

  # check that there are the correct variable names
  assertthat::has_name(region_df, "my_grouping")
  assertthat::has_name(region_df, "batch")
  assertthat::has_name(region_df, "group")

  # check that against_group is in group var
  assertthat::are_equal(against_group %in% region_df$group, TRUE)

  # check that heatmap is either TRUE or FALSE
  assertthat::are_equal(heatmap == TRUE | heatmap == FALSE, TRUE)

  # calculate pvals ---------------------------------------------------------

  # get all t0 brain areas per block
  translating_factor <- abs(min(region_df[,"cells_perthousand_box_scaled"], na.rm = TRUE))

  region_df_transformed <- region_df %>%
    dplyr::mutate(cells_perthousand_box_scaled = cells_perthousand_box_scaled + translating_factor)

  baseline <- region_df_transformed %>%

    # select only baseline group
    dplyr::filter(group == against_group) %>%

    # select vars of interest for later re-merging
    dplyr::select(batch, my_grouping, cells_perthousand_box_scaled) %>%
    dplyr::rename(count_baseline = cells_perthousand_box_scaled)

  # merge
  pairwise_df <- region_df_transformed %>%

    # remove baseline
    dplyr::filter(group != against_group) %>%

    # remerge with t0 dataset
    dplyr::full_join(baseline, by = c("batch", "my_grouping"))

  # t-test approach
  pval_df <- pairwise_df %>%

    # group by time and brain area
    dplyr::group_by(group, my_grouping) %>%

    # get test statistics
    dplyr::summarize(
      p_val = t.test(cells_perthousand_box_scaled, count_baseline,
                     alternative = "greater")$p.value) %>%

    # adjust the pvalue
    dplyr::ungroup() %>%
    dplyr::mutate(p_val_adj = p.adjust(p_val,method="BH")) %>%

    # beautify
    dplyr::mutate(group = paste(group, "vs", against_group)) %>%
    unique()
  ## until here one function

    if (heatmap == FALSE) {
      return(pval_df)
    } else {

      # heatmap plot
      visual <- pval_df %>%

        # calculate -log10
        dplyr::mutate(p_val_adj = -log10(p_val_adj)) %>%

        # visualization
        ggplot2::ggplot(ggplot2::aes(group, my_grouping, fill = p_val_adj)) +

        ggplot2::geom_tile() +

        ggplot2::labs(x = "Comparison", y = "Brain areas", fill = "-log10(adjusted p value)") +

        ggplot2::theme_bw()

      return(
        list(
          data = pval_df,
          visualization = visual
        )
      )


    }



}

