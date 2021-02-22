test_that("calculate_density() example works", {
   x <- data.frame(
   sample_id = rep(c(1:3), each = 100),
   xPos = c(rnorm(100), rpois(100, 2), sample(100)),
   yPos = c(rnorm(100), rpois(100, 2), sample(100)),
   zPos = c(rnorm(100), rpois(100, 2), sample(100))
   )

   y <- calculate_density(x)

   testthat::expect_equal(nrow(y), 300)

})



test_that("plotting_all_coordinates_2d() example works", {

  x <- data.frame(
  sample_id = rep(c(1:3), each = 100),
  xPos = c(rnorm(100), rpois(100, 2), sample(100)),
  yPos = c(rnorm(100), rpois(100, 2), sample(100)),
  zPos = c(rnorm(100), rpois(100, 2), sample(100)),
  group = rep(c("a", "b", "c"), each = 100)
  )

  my_cols <- c("black", "green", "red")
  g <- plotting_all_coordinates_2d(x, "group", my_cols)

  testthat::expect_equal(length(g), 3)

})



test_that("voxelize() example works", {

  x <- data.frame(
    xPos = rnorm(100, 300, 100), yPos = rnorm(100, 300, 100), zPos = rnorm(100, 300, 100),
    group = rep(c("control", "exp"), each = 50),
    sample_id = rep(c(1:10), each = 10)
  )
  y <- voxelize(x, voxel_size = 100)

  testthat::expect_equal(length(y), 2)

})
