#' Calculate utility/favorability for ranking prognoses
#'
#' Expects output like that from \code{\link{pp_scr_xnew}}
#'
#' @param z Potential treatment level under which to calculate favorability
#' @param pp N x 10 x B posterior prediction array
#' @param disc_fact Discount factor for time spent in ill state
#' @return N x B matrix of of utilities under Z = z
#' @export
calc_util <- function(z, pp, disc_fact = 0.001) {
  if (z == 0) {
    # u <- pp[, "yr0_imp"] + (1 - disc_fact) * (pp[, "yt0_imp"] - pp[, "yr0_imp"])
    u <- pp[, 3] + (1 - disc_fact) * (pp[, 4] - pp[, 3])
  } else if (z == 1) {
    # u <- pp[, "yr1_imp"] + (1 - disc_fact) * (pp[, "yt1_imp"] - pp[, "yr1_imp"])
    u <- pp[, 7] + (1 - disc_fact) * (pp[, 8] - pp[, 7])
  }
  return(u)
}
