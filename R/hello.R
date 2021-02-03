# Hello, world!
#
# This is an example function named 'hello'
# which prints 'Hello, world!'.
#
# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Install Package:           'Cmd + Shift + B'
#   Check Package:             'Cmd + Shift + E'
#   Test Package:              'Cmd + Shift + T'

hello <- function() {
  print("Hello, world!")
}

# to build the package
# devtools::document()
# devtools::load_all()
# Load a function with Ctrl/Cmd + Shift + L

# if magrittr does not work, put this in the NAMESPACE file
#importFrom(magrittr,"%>%")

# test
## create a test directory
#usethis::use_testthat()

# put your tests in the tests/testthat/ folder. The files must start with test_
# then bundle with devtools::document()
# then run test with devtools::test()
