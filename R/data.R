
# Atlas tree --------------------------------------------------------------
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


# mask --------------------------------------------------------------------
#' Mask for data cleaning of whole-brain microscopy.
#'
#' A list containing two datasets to remove the edges of the brain (dataframe edge) and the space
#' around the ventricles (dataframe ventricles) for cleaning whole-brain microscopy images. The
#' datasets are on the Allen Brain Reference Atlas (25µm, mouse) brain. The abc4d package uses functions to
#' remove these coordinates from the samples.
#'
#' @format A data frame with 1205 rows and 5 variables:
#' \describe{
#'   \item{edge}{locations around edges of the brain, size of 75µm, consistent with ABA (25µm, mouse)}
#'   \item{edge}{locations around the ventricles, size of 75µm, consistent with ABA (25µm, mouse)}
#' }
#' @source \url{add link to git}
"mask"



# Atlas estimation --------------------------------------------------------------
#' Atlas with estimation of number of cells per brain area by Erö et al.
#'
#' A dataset published by Erö and colleagues (doi: 10.3389/fninf.2018.00084) containing an estimation of number
#' of cells per brain areas in Allen Brain Reference atlas (25µm, mouse) space.
#'
#' @format A data frame with 955 rows and 10 variables:
#' \describe{
#'   \item{Regions}{name of the brain areas, consistent with ABA (25µm, mouse), corresponds to name var in data("atlas_tree")}
#'   \item{Cells}{number of total cells estimated in each brain area}
#'   \item{Neurons}{number of neurons estimated in each brain area}
#'   \item{Glia}{number of glia estimated in each brain area}
#'   \item{Excitatory}{number of excitatory cells estimated in each brain area}
#'   \item{Inhibitory}{number of inhibitory cells estimated in each brain area}
#'   \item{Modulatory}{number of modulatory cells estimated in each brain area}
#'   \item{Astrocytes}{number of astrocytes estimated in each brain area}
#'   \item{Olygodendrocytes}{number of olygodendrocytes estimated in each brain area}
#'   \item{Microglia}{number of microglia estimated in each brain area}
#' }
#' @source \url{https://doi.org/10.3389/fninf.2018.00084}
"atlas_estimation"


# expression cfos from ABA API --------------------------------------------------------------
#' Expression of cfos from Allen Brain Atlas API
#'
#' A dataset with the tabulate file downloaded from the Allen Brain Atlas API of cfos expression.
#' The values refer to ISH staining. For more information, see [this link](http://help.brain-map.org/display/api/Quantified+Data+by+Structures).
#'
#' @format A data frame with 1205 rows and 5 variables:
#' \describe{
#'   \item{id}{brain area id, consistent with ABA (25µm, mouse)}
#'   \item{name}{name of the brain areas, consistent with ABA (25µm, mouse)}
#'   \item{acronym}{acronym of the brain area, consistent with ABAs (25µm, mouse)}
#'   \item{parent_acronym}{acronym of the parent structure of the brain area, ie. upstream structure in the ABA. Consistent with acronym var}
#'   \item{category}{variable internally used}
#' }
#' @source \url{http://mouse.brain-map.org/search/show?page_num=0&page_size=20&no_paging=false&exact_match=true&search_term=cFos&search_type=gene}
"cfos_expression_aba"



# example --------------------------------------------------------------
#' Miniature example of whole-brain imaging data over time from Bonapersona et al.
#'
#' A dataset containing the cleaned coordinates of c-fos positive cells in a selection of brain areas
#' (amygdala and hippocampus) after foot-shock. Experimental groups are: 30, 90, 180 minutes after
#' foot-shock and a sham control (0 min). For more information, see https://osf.io/8muvw/ .
#'
#' @format A data frame in long format with each row indicating a c-fos+ cells. Contains 177885 rows and 7 variables:
#' \describe{
#'   \item{sample_id}{unique identified for each sample}
#'   \item{xPos}{position on the x-axis in horizontal sectioning, consistent with ABA (25µm, mouse)}
#'   \item{yPos}{position on the y-axis in horizontal sectioning, consistent with ABA (25µm, mouse)}
#'   \item{zPos}{position on the z-axis in horizontal sectioning, consistent with ABA (25µm, mouse)}
#'   \item{maxInt}{maximum expression level of c-fos measured per cell, estimated with Imaris}
#'   \item{my_grouping}{grouping of the brain areas based on the researchers specified spatial resolution, consistent with acronyms from ABA (25µm, mouse)}
#'   \item{hemisphere}{hemisphere}
#' }
#' @source \url{https://osf.io/8muvw/}
"example"

