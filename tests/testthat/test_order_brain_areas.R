test_that("find_intercept() example works", {
  x <- find_intercept(x = c(-1.1, 1.1) , y = c(-1.1, 0.9))

  expect_equal(x, -0.1)

})

test_that("order_brain_areas() example works", {
  x <- data.frame(
    batch = rep(c(1:3), each = 15),
    group = rep(rep(c("0", "30", "90"), each = 5), 3),
    my_grouping = c("CA1", "CA2", "CA3", "DG", "BLA"),
    cells_perthousand_box_scaled_ba = rnorm(45, 1.5, 2)
  )

  y <- order_brain_areas(x, order_of_groups = c("0", "30", "90"))
  expect_equal(x, -0.1)

})




y <- sim_most_active_one_batch(x)
testthat::expect_equal(nrow(y), 1)
