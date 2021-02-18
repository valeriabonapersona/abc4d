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
