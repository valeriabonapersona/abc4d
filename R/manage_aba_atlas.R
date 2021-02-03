'
  Functions to find parents / children within the ABA tree structure
'


# find_children() ------------------------------------------------------------
#' @title Find smaller areas (children) of a larger brain area categorization
#' (parent) in the Allen Brain Atlas.
#'
#'
#' @description The Allen Brain Atlas is organized with a tree structure, where
#'  smaller brain areas (e.g. nuclei) are grouped into larger brain areas.
#'  Given a particular (parent) acronym, this function finds related sub-brain areas,
#'  which are returned in a string. The function also returns children that are
#'  parents themselves, but it does not return the smallest categorization of brain areas.
#'  To find the smallest categorizations, use the find_all_children().
#'
#'  If the parent acronym provided is not present in atlas, the function returns an empty string.
#'
#'
#' @param atlas Data frame with (at least) two character variables: one for the parents (parent_acronym),
#' and one for the children (acronym).
#' @param parent Acronym of the parent brain area. It must be a character vector of
#' length 1.
#'
#' @return
#'
#' @examples
#' x <- data.frame(
#' parent_acronym = c(rep("main_parent", 3), rep("child_with_children", 2), "grandchild_with_children"),
#' acronym = c("main_parent", "child_with_children", "child_without_children", "grandchild_with_children", "grandchild_without_children", "great_grandchild"))
#' this_parent <- "child_with_children"
#'
#' find_children(atlas = x, parent = this_parent)
#' # to return also "great_grandchild", you find_all_children
#'
find_children <- function(atlas, parent) {

  # check atlas is a dataframe
  assertthat::assert_that(is.data.frame(atlas))

  # check parent_acronym and acronym are the var names of the df
  assertthat::has_name(atlas, "parent_acronym")
  assertthat::has_name(atlas, "acronym")

  # check parent has length 1
  assertthat::are_equal(length(parent), 1)


  # get the children
  atlas %>%
    dplyr::filter(parent_acronym == parent) %>%
    dplyr::pull(acronym)

}


# find_all_children() -------------------------------------------------------

#' @title Find smaller areas (children) of a larger brain area categorization
#' (parent) in the Allen Brain Atlas.
#'
#' @description The Allen Brain Atlas is organized with a tree structure, where
#'  smaller brain areas (e.g. nuclei) are grouped into larger brain areas. Given
#'  a particular (child) acronym, this function finds all sub-brain areas, and
#'  returns them as a string. This function returns children that are parents themselves,
#'  as well their children (i.e., grandchildren (etc.) of the parent structure provided).
#'
#' @param atlas Data frame with (at least) two character variables: one for the parents (parent_acronym),
#' and one for the children (acronym).
#' @param parent Acronym of the parent brain area. It must be a character vector of
#' length 1.
#'
#' @return
#' @export
#'
#' @examples
#' x <- data.frame(
#' parent_acronym = c(rep("main_parent", 3), rep("child_with_children", 2), "grandchild_with_children"),
#' acronym = c("main_parent", "child_with_children", "child_without_children", "grandchild_with_children", "grandchild_without_children", "great_grandchild"))
#' this_parent <- "child_with_children"
#'
#' find_all_children(atlas = x, parent = this_parent)
#'

find_all_children <- function(atlas, parent) {

  # check atlas is a dataframe
  assertthat::assert_that(is.data.frame(atlas))

  # check parent_acronym and acronym are the var names of the df
  assertthat::has_name(atlas, "parent_acronym")
  assertthat::has_name(atlas, "acronym")

  # check parent has length 1
  assertthat::are_equal(length(parent), 1)

  # check parent is in parent_acronym
  if(!parent %in% atlas$parent_acronym) stop('Parent acronym provided not found in atlas')

  # find first generation
  f_1 <- find_children(atlas, parent) # find children of this parent
  save_f <- f_1

  # repeat until no more children are found
  repeat {

    f_n <- unlist(lapply(f_1, find_children, atlas = atlas)) # find grandchildren
    save_f <- c(save_f, f_n)
    f_1 <- f_n

    if (is.null(f_n)) break
  } # repeat until no more found

  return(save_f)
}




# find_parents() ----------------------------------------------------------
#' @title Find larger grouping areas (parents) of a smaller brain area
#' (child) in the Allen Brain Atlas.
#'
#' @description The Allen Brain Atlas is organized with a tree structure, where
#'  smaller brain areas (e.g. nuclei) are grouped into larger brain areas. This function
#'  allows you to find the parent of a sub-brain area of choice. Returns only one level up
#'  in the tree structure.
#'
#' @param atlas Data frame with (at least) two character variables: one for the parents (parent_acronym),
#' and one for the children (acronym).
#' @param child Acronym of the child brain area. It must be a character vector of
#' length 1.
#'
#' @export
#' @return
#'
#' @examples
#'#' x <- data.frame(
#' parent_acronym = c(rep("main_parent", 3), rep("child_with_children", 2), "grandchild_with_children"),
#' acronym = c("main_parent", "child_with_children", "child_without_children", "grandchild_with_children", "grandchild_without_children", "great_grandchild"))
#' this_child <- "child_with_children"
#'
#' find_parents(atlas = x, child = this_child)


find_parents <- function(atlas, child) {

  # check atlas is a dataframe
  assertthat::assert_that(is.data.frame(atlas))

  # check parent_acronym and acronym are the var names of the df
  assertthat::has_name(atlas, "parent_acronym")
  assertthat::has_name(atlas, "acronym")

  # check parent has length 1
  assertthat::are_equal(length(child), 1)

  # check parent is in parent_acronym
  if(!child %in% atlas$acronym) stop('Child acronym provided not found in atlas')

  atlas %>%
    dplyr::filter(acronym %in% child) %>%
    dplyr::pull(parent_acronym)

}





# # Wrappers ----------------------------------------------------------------

# find_all_children_df <- function(atlas, these_parents) {
#   res_ls <- lapply(these_parents, function(x) {find_all_children(atlas, x)})
#   names(res_ls) <- these_parents
#   res <- reshape2::melt(res_ls)
#   names(res) <- c("acronym", "parents")
#   return(res)
# }
#
#
