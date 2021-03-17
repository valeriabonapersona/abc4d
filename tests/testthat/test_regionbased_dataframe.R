
test_that("adapt_estimation_atlas example works", {

  # dataset
  x <- data.frame(
    Regions = c("a_1", "a_2", "a_3", "b_1", "b_2"),
    Cells = c(10000, 2100, 39847, 754, 923)
  )
  y <- data.frame(
    name = c("a_1", "a_2", "a_3", "b_1", "b_2"),
    my_grouping = c(rep("a",3), rep("b",2))
  )

  # run
  z <- adapt_estimation_atlas(x,y)

  correct_res <- dplyr::tibble(my_grouping = c("a", "b"),
                               cells = c(51947, 1677))
  expect_equal(z, correct_res)


})



test_that("summarize_per_region example works", {

  # create dataframe
  x <- data.frame(
    sample_id = c(rep("a", 100), rep("b", 75), rep("c", 50)),
    my_grouping = c(rep(c("CA1","CA2","CA3", "DG", "BLA"), each = 20),
                    rep(c("CA1", "CA2", "CA3"), each = 20),
                    rep("BLA", 15),
                    rep(c("DG", "BLA"), each = 25)),
    maxInt = abs(rnorm(225, 100))
   )

  y <- data.frame(
    my_grouping = c("CA1","CA2","CA3", "DG", "BLA"),
    cells = sample(100000, 5),
    glia = sample(10000, 5)
  )

  z <- data.frame(
    sample_id = c("a", "b", "c"),
    batch = c(1,1,2),
    group = c("A", "B", "A")
    )
  res <- summarize_per_region(x,y, cells_type = "cells")

  expect_equal(nrow(res), 11)


})

test_that("preprocess_per_region example works", {

  # create dataframe
   x <- data.frame(
       batch = rep(c(1,1,2,2), each = 5),
       group = rep(c("control", "exp", "exp", "control"), each = 5),
       sample_id = rep(c("a", "b", "c", "d"), each = 5),
       my_grouping = rep(c("CA1", "CA2", "CA3", "DG", "BLA"), 4),
       intensity = sample(10000, 20, replace = TRUE),
       cells_perthousand = abs(rnorm(20))
   )

   z <- preprocess_per_region(x)
   expect_equal(is.data.frame(z), TRUE)

})
