test_that("manage_aba_atlas functions find correct output", {

  # dataset
  x <- data.frame(
    parent_acronym = c(rep("main_parent", 3), rep("child_with_children", 2), "grandchild_with_children"),
    acronym = c("main_parent", "child_with_children", "child_without_children", "grandchild_with_children",
                "grandchild_without_children", "great_grandchild")
  )

  # expectations
  expect_equal(find_children(x, "child_with_children"), c("grandchild_with_children", "grandchild_without_children"))
  expect_equal(find_all_children(x, "child_with_children"), c("grandchild_with_children", "grandchild_without_children", "great_grandchild"))
  expect_equal(find_parents(x, "child_with_children"), "main_parent")
})

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
