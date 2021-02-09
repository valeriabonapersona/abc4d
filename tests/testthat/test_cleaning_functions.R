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

test_that("remove_blinding() returns a dataframe with the correct columns' names", {

  x <- data.frame(sample_id = c(1:20), blinded_group = c(rep(c(1:4), each = 5)))
  y <- data.frame(blinded_group = 1:4, key = c("control", "drug_a", "drug_b", "drug_c"))

  z <- remove_blinding(x,y)
  expect_named(z, c("sample_id", "group"))

})

