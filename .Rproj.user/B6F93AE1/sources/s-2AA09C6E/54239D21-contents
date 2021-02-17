test_that("pairwise_analysis() example works", {

  dat <- data.frame(
    batch = rep(c(1:3), each = 150),
    group = rep(rep(c("0", "30", "90"), each = 50), 3),
    my_grouping = rep(rep(c(1:50), 3), 3),
    cells_perthousand_box_scaled = rnorm(150*3)
   )
  res <- pairwise_activation(region_df = dat)

  res_with_visualization <- pairwise_activation(dat, heatmap = TRUE)

  expect_equal(res, res_with_visualization[["data"]])
  expect_equal(nrow(res), 100)
  expect_equal(sum(is.na(res)), 0)

})
