
# find_intercept ----------------------------------------------------------
#' @title Function to find where an interpolated linear curve crosses a y value
#'
#' @description The function interpolates a linear model within two points, and finds at
#' which x the curve has a specified y.
#'
#' @param x number of x variable
#' @param y number of y variable
#' @param inter number of y at which you want to calculate x
#'
#' @return
#' @export
#'
#' @examples
#'find_intercept(x = c(-1.1, 1.1) , y = c(-1.1, 0.9))


find_intercept <- function(x, y, inter = 0) {
  approximate_function <- approxfun(x, y)
  approximate_function(inter)
}

