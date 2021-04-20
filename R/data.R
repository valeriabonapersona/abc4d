#' Allen Brain Reference Atlas tree structure of brain areas.
#'
#' A dataset containing the tree structure of the brain areas in the 25µm mouse Allen Brain Reference Atlas (ABA)
#'
#' @format A data frame with 1205 rows and 5 variables:
#' \describe{
#'   \item{id}{brain area id, consistent with ABA (25µm, mouse)}
#'   \item{name}{name of the brain areas, consistent with ABA (25µm, mouse)}
#'   \item{acronym}{acronym of the brain area, consistent with ABAs (25µm, mouse)}
#'   \item{parent_acronym}{acronym of the parent structure of the brain area, ie. upstream structure in the ABA. Consistent with acronym var}
#'   \item{category}{variable internally used}
#' }
#' @source \url{http://help.brain-map.org/display/api/Atlas+Drawings+and+Ontologies}
"atlas_tree"

