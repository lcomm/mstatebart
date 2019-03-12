#' Function to turn list of tree inputs and BART fit objects into posterior
#' predictive array of exit times and dispositions for the population
#' used to fit the tree (\code{factual_imputation = TRUE}) or the design matrix
#' for the other treatment arm (\code{factual_imputation = FALSE})
#'
#' @param tree_factual_inputs List of inputs for compartment
#' @param tree_cfactual_inputs List of inputs for compartment, for group not
#' used to fit the tree
#' @param tree1_factual_fit BART fit object from time-to-exit fit
#' @param tree2_factual_fit BART fit object for cause-of-exit fit
#' @param factual_imputation Whether to impute for population used to fit the
#' trees, or for the counter-to-fact group. If \code{factual_imputation = TRUE},
#' the imputation will only be done for those whose exit time was censored.
#' @param mc.cores Number of cores to use in predictions
#' @return N x 2 x B array with (i, 1, b) element equal to the exit time (note:
#' may be \code{Inf}) for ith person in population to predict according to the
#' bth posterior draw from the ensemble. Whether the exit was due to Cause 1
#' is recorded in the element (i, 2, b)
#' @export
make_crisk_exit_array <- function(tree_factual_inputs,
                                  tree_cfactual_inputs,
                                  tree1_factual_fit,
                                  tree2_factual_fit,
                                  factual_imputation = FALSE,
                                  mc.cores = 4) {

  # Same regardless of population we are imputing for
  time_seq <- sort(unique(tree_factual_inputs[[1]][["times"]]))

  # Make design matrices based on whether we are imputing in factual
  # arm or not
  if (factual_imputation) {
    censored <- (tree_factual_inputs[[1]][["delta"]] == 0)
    stopifnot(sum(censored) > 0)
    x_new <- tree_factual_inputs[[1]][["xnew"]][censored, ]
    min_times <- tree_factual_inputs[[1]][["times"]][censored]
  } else {
    x_new <- tree_cfactual_inputs[[1]][["x_train"]]
    min_times <- rep(0, nrow(x_new))
  }

  # Predict exit times from all B ensembles
  cf_time_exit <- predict_times_from_pbart(pbart_fit = tree1_factual_fit,
                                           time_seq = time_seq,
                                           x_new = x_new,
                                           min_times = min_times,
                                           mc.cores = mc.cores)

  # Make array of predictions for each posterior ensemble
  B <- tree1_factual_fit$ndpost
  exit_array <- array(NA, c(nrow(x_new), 2, B))
  for (b in 1:B) {
    imputed_times <- cf_time_exit[, b]
    finite_time <- which(is.finite(imputed_times))
    tx_new <- cbind(imputed_times[finite_time], x_new[finite_time, ])
    colnames(tx_new) <- colnames(tree_cfactual_inputs[[2]][["x_train"]])
    cf_cause_exit <- predict_is_cause1(pbart_fit = tree2_factual_fit,
                                       x_new = tx_new,
                                       mc.cores = mc.cores)[b, ]
    exit_time <- imputed_times
    is_cause_1 <- imputed_times * NA
    is_cause_1[finite_time] <- cf_cause_exit
    if (factual_imputation) {
      exit_time[!censored] <- tree_factual_inputs[[1]][["times"]][!censored]
      is_cause_1[!censored] <- tree_factual_inputs[[1]][["delta"]][!censored]
    }
    exit_array[ , 1, b] <- exit_time
    exit_array[ , 2, b] <- is_cause_1
  }

  return(exit_array)
}


