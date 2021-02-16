test_that("prepare_sim_weights() example works", {

  # create dataframes
    x <- data.frame(
      batch = rep(c(1,1,2,2), each = 5),
      group = rep(c("control", "exp", "exp", "control"), each = 5),
      sample_id = rep(c("a", "b", "c", "d"), each = 5),
      my_grouping = rep(c("CA1", "CA2", "CA3", "DG", "BLA"), 4),
      intensity_ave = sample(10000, 20, replace = TRUE),
      cells_perthousand = abs(rnorm(20))
      )

    y <- data.frame(
      my_grouping = c("CA1", "CA2", "CA3", "DG", "BLA"),
      parents = c(rep("hippocampus",4), "cortical subplate")
    )

    z <- data.frame(
      acronym = c("CA1", "CA2", "CA3", "DG", "BLA"),
      mean_expression = rnorm(5, 10, 1),
      sd_expression = abs(rnorm(5))
    )

    # run
  q <- prepare_sim_weights(x,y,z, "control")
  # check
  testthat::expect_equal(nrow(q), 10)
})
