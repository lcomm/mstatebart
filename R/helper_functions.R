#' Function to extract unique values from vector and sort in ascending order
#'
#' @param x Vector with values to sort
#' @return Vector of sorted unique values
#' @export
sort_unique <- function(x) {
  stopifnot(is.vector(x))
  return(sort(unique(x)))
}



#' Array binding
#'
#' @param ... Parameters to pass to \code{\link{abind::abind}}
#' @return 3-dimensional array
#' @export
acomb <- function(...) {
  abind::abind(..., along = 3)
}

