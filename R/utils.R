'

  Utilities and Dependencies

'
# I also installed: devtools, rmarkdown, usethis
usethis::use_package('magrittr') # for general data manipulation
usethis::use_package('dplyr') # for general data manipulation
usethis::use_package('tidyr') # for drop_na in find_most_active
usethis::use_package('ggplot2') # for visualizations (heatmap)
usethis::use_package('stringr') # for manipulation to strings
usethis::use_package('assertthat') # for checks in functions
usethis::use_package('rmarkdown') # is this necessary for vignettes?
usethis::use_package('recipes') # for boxcox transformation
usethis::use_package('ggrepel') # to repel labels in ggplot visual
usethis::use_package('ks') # for kernel densities
usethis::use_package("effsize") # to calculate effect sizes

