#' Run a cause-of-exit BART to model the reason for a compartment exit
#'
#' Second model of "Method 1" from competing risks BART paper
#'
#' @param input_list List of inputs containing 0/1 indicator of whether the exit
#' was due to cause 1 (\code{is_cause1}) as well as design matrix (\code{x_train})
#' @param mc.cores Number of cores to use in fitting process
#' @param ... Other parameters to pass to \code{\link{bart::mc.gbart}}
#' @return probit BART fit object
#' @export
exit_cause_bart_from_input_list <- function(input_list, mc.cores = 4, ...) {
  res <- mc.gbart(y.train  = input_list[["is_cause1"]],
                  x.train  = input_list[["x_train"]],
                  type     = "pbart",
                  mc.cores = mc.cores,
                  ...)
  return(res)
}



#' Predict whether or not the compartment exit is due to Cause 1
#' based on a probit BART fit
#'
#' Second model of "Method 1" for competing risks BART
#'
#' @param pbart_fit BART fit object with B posterior ensemble draws
#' @param x_new New data set of N observations for which to make predictions
#' @param mc.cores Number of cores to use
#' @return N x B 0/1 matrix for whether the exit was due to Cause 1
#' @export
predict_is_cause1 <- function(pbart_fit, x_new, mc.cores = 4) {
  B <- pbart_fit$ndpost
  pred <- predict(object = pbart_fit, newdata = x_new, mc.cores = mc.cores)
  coin_flips <- apply(pred$prob.test,
                      MARGIN = c(1, 2),
                      FUN = rbinom,
                      size = 1, n = 1)
  return(coin_flips)
}



#' Predict whether or not the compartment exit is due to Cause 1
#' based on a probit BART fit for a single MCMC iteration
#'
#' Second model of "Method 1" for competing risks BART
#'
#' @param pbart_fit BART fit object with B posterior ensemble draws
#' @param x_new New data set of N observations for which to make predictions
#' @param mc.cores Number of cores to use
#' @return N x B 0/1 matrix for whether the exit was due to Cause 1
#' @export
predict_is_cause1_b <- function(b, pbart_fit, x_new) {
  stopifnot(1 <= b, b <= pbart_fit$ndpost)
  pred <- predict_prob_b(b = b, object = pbart_fit, newdata = x_new)
  coin_flips <- rbinom(n = length(pred), size = 1, pred)
  return(coin_flips)
}
