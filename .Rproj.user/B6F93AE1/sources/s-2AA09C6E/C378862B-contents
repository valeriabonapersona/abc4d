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
      acronym = c("hippocampus", "cortical subplate"),
      mean_expression = rnorm(2, 10, 1),
      sd_expression = abs(rnorm(2))
    )

    # run
  q <- prepare_sim_weights(x,y,z, "control")
  # check
  testthat::expect_equal(nrow(q), 10)
})


test_that("sim_most_active_one_batch() example works", {

  # WITH EXPRESSION LEVELS
  # create dataframes
  x <- data.frame(
   group = rep(c("control", "experimental"), each = 5),
   my_grouping = rep(c("CA1", "CA2", "CA3", "DG", "BLA"), 2),
   mean_expression = c(rnorm(5, 10, 2), rnorm(5, 13, 2)),
   sd_expression = abs(rnorm(10)),
   weight = c(rep(1, 5), rnorm(5, 3, 1))
   )

  y <- sim_most_active_one_batch(x)
  testthat::expect_equal(nrow(y), 1)


  # WITHOUT EXPRESSION LEVELS
  z <- sim_most_active_one_batch(x, weight_by_expression = FALSE)

  # check
  testthat::expect_equal(nrow(z), 1)
})


test_that("sim_most_active() example works", {

  # WITH EXPRESSION LEVELS
  # create dataframes
  x <- data.frame(
    group = rep(c("control", "experimental"), each = 5),
    my_grouping = rep(c("CA1", "CA2", "CA3", "DG", "BLA"), 2),
    mean_expression = c(rnorm(5, 10, 2), rnorm(5, 13, 2)),
    sd_expression = abs(rnorm(10)),
    weight = c(rep(1, 5), rnorm(5, 2, 1))
  )

  y <- sim_most_active(x, samples_per_group = 3)
  testthat::expect_equal(nrow(y), 3)


  # WITHOUT EXPRESSION LEVELS
  x <- data.frame(
    group = rep(c("control", "experimental"), each = 5),
    my_grouping = rep(c("CA1", "CA2", "CA3", "DG", "BLA"), 2)
  )
  z <- sim_most_active(x, samples_per_group = 3, n_exp = 3, weight_by_expression = FALSE)

  # check
  testthat::expect_equal(nrow(z), 9)
})



#' x <- data.frame(
#' group = rep(c("control", "experimental"), each = 5),
#' my_grouping = rep(c("CA1", "CA2", "CA3", "DG", "BLA"), 2),
#' mean_expression = c(rnorm(5, 10, 2), rnorm(5, 13, 2)),
#' sd_expression = abs(rnorm(10)),
#' weight = c(rep(1, 5), rnorm(5, 3, 1))
#' )
#'
#' sim_most_active_one_batch(x)
