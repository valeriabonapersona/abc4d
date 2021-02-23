
test_that("categorize_strategy() example works", {

  x <- data.frame(
  batch = rep(c(1,1,2,2), each = 5),
  group = rep(c("control", "exp", "exp", "control"), each = 5),
  sample_id = rep(c("a", "b", "c", "d"), each = 5),
  my_grouping = rep(c("CA1", "CA2", "CA3", "DG", "BLA"), 4),
  intensity = sample(10000, 20, replace = TRUE),
  cells_perthousand = abs(rnorm(20))
  )

  y <- categorize_strategy(x, "control")

  testthat::expect_equal("strategy" %in% names(y), TRUE)

})



test_that("find_changed_strategy example works", {


  x <- data.frame(
  batch = rep(c(1,1,2,2), each = 5),
  group = rep(c("control", "exp", "exp", "control"), each = 5),
  sample_id = rep(c("a", "b", "c", "d"), each = 5),
  my_grouping = rep(c("CA1", "CA2", "CA3", "DG", "BLA"), 4),
  intensity = sample(10000, 20, replace = TRUE),
  cells_perthousand = abs(rnorm(20))
  )
  x$intensity_box_scaled <- scale(x$intensity)
  x$cells_perthousand_box_scaled <- scale(x$cells_perthousand)


  y <- find_changed_strategy(x, 2, "control", min_change_rate = 1)
  testthat::expect_equal("strategy_towards" %in% names(y), TRUE)

})
