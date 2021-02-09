'
Test cleaning functions
'

test_that("specify_damage() provides correct output", {

  # first test
  sample <- "5"
  my_damage <- data.frame(
    sample = c(1:10),
    area = factor(sample(10, replace = TRUE)),
    hemisphere = c(rep("left", 7), rep("right", 3))
  )
  res <- specify_damage(sample, my_damage)
  expect_equal("area_hemisphere" %in% names(res), TRUE)


  # second test
   x <- data.frame(
   sample = c(1,1,2,3),
   area = c("BLA", "VTA", "CA1", "DG"),
   hemisphere = c("left", "left", "right", "right")
   )

   my_sample <- 1

   res <- specify_damage(my_sample, x)
   expect_equal(nrow(res), 2)
})

test_that("normalize() modifies string", {

  x <- rpois(100, 6)
  y <- normalize(x)

  expect_equal(identical(x,y), FALSE)

})
