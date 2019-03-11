#' Simulate semicompeting risks with no recurrence and no subject-specific
#' frailty
#'
#' @param n Number of observations
#' @param scenario Parameter scenario from
#' \code{\link{rsemicompstan::return_dgp_parameters}}, except that it overrides
#' the variance of the subject-specific frailty to be (basically) zero
#' @param seed Random number seed
#' @param cens_time Time for administrative right censoring. Default is 90.
#' @export
#' @examples
#' \dontrun{
#' # Make rehosp data set for package
#' rehosp <- simulate_scenario_no_frailty(n = 6000, scenario = 1, seed = 456)
#' new_yr <- ceiling(rehosp$yr)
#' new_yt <- ceiling(rehosp$yt)
#' for (i in which((new_yt - new_yr == 0) & (rehosp$dyr == 1))) {
#'   if (new_yt[i] < 90) {
#'     new_yt[i] <- new_yt[i] + 1
#'   } else if (new_yt[i] == 90) {
#'     new_yr[i] <- new_yr[i] - 1
#'   }
#' }
#' rehosp$yr <- new_yr
#' rehosp$yt <- new_yt
#' }
simulate_scenario_no_frailty <- function(n, scenario, seed = 456, cens_time = 90) {
  params <- rsemicompstan::return_dgp_parameters(scenario = scenario)
  params$sigma <- .Machine$double.eps
  dat <- rsemicompstan::simulate_from_param(n = n, seed = seed,
                                            params = params, cens_times = rep(cens_time, n))
  dat$a <- dat$z
  return(dat)
}
