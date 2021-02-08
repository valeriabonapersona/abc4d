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
