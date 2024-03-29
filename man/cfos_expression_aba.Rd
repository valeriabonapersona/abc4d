% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{cfos_expression_aba}
\alias{cfos_expression_aba}
\title{Expression of cfos from Allen Brain Atlas API}
\format{
A data frame with 1205 rows and 5 variables:
\describe{
  \item{id}{brain area id, consistent with ABA (25µm, mouse)}
  \item{name}{name of the brain areas, consistent with ABA (25µm, mouse)}
  \item{acronym}{acronym of the brain area, consistent with ABAs (25µm, mouse)}
  \item{parent_acronym}{acronym of the parent structure of the brain area, ie. upstream structure in the ABA. Consistent with acronym var}
  \item{category}{variable internally used}
}
}
\source{
\url{http://mouse.brain-map.org/search/show?page_num=0&page_size=20&no_paging=false&exact_match=true&search_term=cFos&search_type=gene}
}
\usage{
cfos_expression_aba
}
\description{
A dataset with the tabulate file downloaded from the Allen Brain Atlas API of cfos expression.
The values refer to ISH staining. For more information, see the Allen Brain Atlas API (http://help.brain-map.org/display/api/Quantified+Data+by+Structures).
}
\keyword{datasets}
