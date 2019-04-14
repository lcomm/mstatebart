#' Make a list of time sequences for BART modesl
#'
#' @param df Data frame containing \code{z}, \code{yr}, \code{yt}, \code{dyr},
#' \code{dyt}
#' @return Length-2 list length-2 lists of unique times under each Z condition
#' for each compartment exit
#' @export
make_time_seqs <- function(df) {
  stopifnot(c("z", "yr", "yt", "dyr", "dyt") %in% colnames(df))

  # Unique time-to-exit time labels
  df0 <- df[df$z == 0, ]
  t0_k1 <- sort_unique(df0$yr)
  t0_k2 <- sort_unique((df0$yt - df0$yr)[df0$dyr == 1])
  t0_k1 <- t0_k1[t0_k1 > 0]
  t0_k2 <- t0_k2[t0_k2 > 0]

  df1 <- df[df$z == 1, ]
  t1_k1 <- sort_unique(df1$yr)
  t1_k1 <- t1_k1[t1_k1 > 0]
  t1_k2 <- sort_unique((df1$yt - df1$yr)[df1$dyr == 1])
  t1_k2 <- t1_k2[t1_k2 > 0]

  time_seqs <- list("0" = list(t0_k1, t0_k2),
                    "1" = list(t1_k1, t1_k2))
  return(time_seqs)
}




#' Make list of minimum time-to-exit for the two compartments for the
#' nonrecurrent semicompeting risks model
#'
#' If drawing from the factual, this will be the censoring times translated
#' into minimum comparment durations. If drawing from the counterfactual,
#' the minimum durations are all zero
#'
#' @param yr Time of (factual) censoring or event for illness
#' @param yt Time of (factual) censoring or event for death
#' @param dyr Indicator for whether illness was observed
#' @param dyt Indicator for whether death was observed
#' @param factual Whether or not the minimum times to event should be made to
#' agree with \code{yr} and \code{yt}
#' @return Length-2 list of minimum ttes for the 2 compartments
#' @export
make_min_ttes <- function(yr, yt, dyr, dyt, factual = TRUE) {
  N_test <- length(yr)
  if (factual) {
    # TODO(LCOMM): REVIEW THIS LOGIC
    min_k1 <- yr
    # Should not predict; we will use observed factual Compartment 1 exit time
    min_k1[pmax(dyr, dyt) == 1] <- Inf
    min_k2 <- pmax(0, yt - yr)
    # Should not predict; we will used observed factual death time
    min_k2[dyt == 1] <- Inf
    min_ttes <- list(min_k1, min_k2)
  } else {
    min_ttes <- list(rep(0, N_test), rep(0, N_test))
  }
  #
  return(min_ttes)
}



#' Draw from the posterior predictive distribution of (censored) event times
#' for the nonrecurrent semicompeting risks model
#'
#' Only performs draws for one treatment condition.
#'
#' @param b MCMC iteration to do prediction for
#' @param z Z condition for which imputations are being done (only used to label
#' variables in resulting output matrix)
#' @param time_seqs_train Length-2 list of time sequences for durations spent in
#' Compartments 1 and 2 (no recurrence)
#' @param x_train Covariate matrix used to train first model
#' @param x_test n-row covariate matrix for which to make predictions
#' @param min_ttes Length-2 list of minimum event durations for each compartment
#' (if trying to match observed censoring times; all zero otherwise)
#' @param fit_k1 Length-2 list of BART models for (1) time to Compartment 1 exit
#' and (2) whether cause of Compartment 1 exit was illness
#' @param fit_k2 Length-1 list of BART model for time to Compartment 2 exit,
#' among those who ever enter Compartment 2
#' @return n x 4 matrix of imputed yr, yt, dyr, and dyt
#' @export
predict_scr_b <- function(b,
                          z,
                          time_seqs_train,
                          x_train, # design matrix used to fit tree
                          x_test, # design matrix needed for prediction
                          min_ttes,
                          fit_k1,
                          fit_k2,
                          mc.cores = 1) {

  # Parse inputs
  N_test <- nrow(x_test)
  N_train <- nrow(x_train)
  max_t <- max(time_seqs_train[[1]])

  # Make shells
  tte_2 <- tte_1 <- cause_1 <- rep(NA, N_test)
# browser()
  # Make the counting process format to feed into Fit 1, compartment 1
  tx_k1 <- surv.pre.bart(times = time_seqs_train[[1]],
                         delta = rep(1, length(time_seqs_train[[1]])),
                         # x.train = NULL,
                         x.train = x_test[rep(1, length(time_seqs_train[[1]])), ],
                         x.test = x_test)$tx.test

  # Predict compartment 1 time-to-exits
  tte_1 <- predict_times_from_pbart_b(pbart_fit = fit_k1[[1]],
                                      b = b,
                                      time_seq = time_seqs_train[[1]],
                                      x_new = x_test,
                                      tx_new = tx_k1,
                                      min_times = min_ttes[[1]])

  # Make boolean vector of when an exit occurs before max_t
  had_event_k1 <- is.finite(tte_1)

  if (any(is.finite(tte_1))) {

    # Make newdata matrix for cause-of-exit model
    tx_new_k1_cause <- cbind(tte_1, x_test)

    # Predict whether the Compartment 1 exit was due to Cause 1 (illness)
    cause_1[had_event_k1] <- predict_is_cause1_b(pbart_fit = fit_k1[[2]],
                                                 b = b,
                                                 x_new =
                                                   tx_new_k1_cause[had_event_k1, ])

    if (any(cause_1 == 1)) {

      # ID who moved to Compartment 2
      # browser()
      in_k2 <- (cause_1 == 1)
      win_k2 <- which(in_k2)

      # Predict Compartment 2 exit times
      tte_2[win_k2] <- predict_times_from_pbart_b(pbart_fit = fit_k2,
                                                  b = b,
                                                  time_seq = time_seqs_train[[2]],
                                                  x_new = tx_new_k1_cause[win_k2, ],
                                                  min_times = min_ttes[[2]][win_k2])
    #branch for if none of Compartment 1 exits are due to illness
    } else {
      warning("No healthy-ill transitions imputed...")
    }

    tte_2[is.na(cause_1)] <- 0 # zero sojourn
    tte_2[cause_1 == 0] <- 0
    # browser()
    death_time <- tte_1 + tte_2
    yt_imp <- pmin(death_time, max_t)
    dyt_imp <- (death_time <= max_t) * 1
    yr_imp <- pmin(tte_1, max_t)
    dyr_imp <- 1 - (cause_1 %in% c(NA, 0))

  #branch for if no Compartment 1 exits are imputed at all
  } else {
    warning("No exit events imputed...")
    # browser()
    yr_imp <- yt_imp <- rep(max_t, NROW(N_test))
    dyr_imp <- dyt_imp <- rep(0, NROW(N_test))
  }

  # Overwrite times if actually observed
  # browser()
  # yr_imp[dyr == 1] <- min_ttes[[1]][dyr == 1]
  # yt_imp[dyt == 1] <- pmax(min_ttes[[1]][dyt == 1],
  #                          min_ttes[[1]][dyt == 1] +
  #                            min_ttes[[2]][dyt == 1])

  # Final prettying up
  labs <- paste0(c("yr", "yt", "dyr", "dyt"), z, "_imp")
  res <- matrix(NA, nrow = N_test, ncol = 4)
  colnames(res) <- labs
  res[ , 1] <- yr_imp
  res[ , 2] <- yt_imp
  res[ , 3] <- dyr_imp
  res[ , 4] <- dyt_imp
  # browser()
  # head(res)
  return(res)
}



#' Do all posterior predictions for a treatment subset (those observed to have
#' level \code{Z}) using the bth MCMC iteration
#'
#' @param b MCMC iteration for prediction
#' @param x Named list of design for Z=0 and Z=1 conditions
#' @param time_seqs Named list of time-to-event sequences for model evaluation
#' times
#' @param Z Scalar 0/1 subgroup for which predictions are desired
#' @param min_ttes Named list of minimum time-to-events (only used if imputing
#' data missing due to censoring in the factual arm)
#' @param fit_k1 Named list of lists containing BART fits for time and disposition
#' of exit from Compartment 1
#' @param fit_k2 Named list of BART fits for exiting Compartment 2
#' @return N_test x 10 matrix of predicted times and observation indicators, where
#' N_test is the number of subjects with (factual) treatment equal to Z
#' @export
ppredict_scr_b <- function(b, x,
                           yr, yt, dyr, dyt,
                           time_seqs, Z, min_ttes, fit_k1, fit_k2) {
  charZ <- as.character(Z)

  res <- matrix(NA, nrow = NROW(x[[charZ]]), ncol = 8)
  if (Z == 0) {
    res[, 1] <- yr
    res[, 2] <- yt
    res[, 3] <- dyr
    res[, 4] <- dyt
    res[, 5:8] <- predict_scr_b(b = b, z = 1,
                                time_seqs_train = time_seqs[["1"]],
                                x_train = x[["1"]], # design matrix used to fit tree
                                x_test = x[[charZ]], # design matrix needed for prediction
                                min_ttes = min_ttes[[charZ]][["1"]],
                                fit_k1 = fit_k1[["1"]],
                                fit_k2 = fit_k2[["1"]])
  } else if (Z == 1) {
    res[, 1:4] <- predict_scr_b(b = b, z = 0,
                                time_seqs_train = time_seqs[["0"]],
                                x_train = x[["0"]], # design matrix used to fit tree
                                x_test = x[[charZ]], # design matrix needed for prediction
                                min_ttes = min_ttes[[charZ]][["0"]],
                                fit_k1 = fit_k1[["0"]],
                                fit_k2 = fit_k2[["0"]])
    res[, 5] <- yr
    res[, 6] <- yt
    res[, 7] <- dyr
    res[, 8] <- dyt
  }

  return(res)
}



#' Do BART (non-recurrent) semicompeting risks posterior predictive draws for
#' all subjects, all z levels, and all MCMC iterations
#'
#' @param df Data frame with treatment indicator \code{z}, illness and death
#' times \code{yr} and \code{yt}, and observation indicators \code{dyr} and
#' \code{dyt}
#' @param x Named list of covariate design matrices for each treatment level
#' @param fit_k1 List of lists containing BART fits for Compartment 1 stays
#' (1 = time-to-exit model, 2 = cause-of-exit model)
#' @param fit_k2 List of BART fits for Compartment 2 stays
#' @return N x 10 x B posterior prediction array
#' @export
posterior_predict_scr_bart <- function(df, x, fit_k1, fit_k2) {

  # Need all models to get predictions
  # Therefore can only predict for the minimum number of post-warmup draws
  B <- min(fit_k1[["0"]][[1]][["ndpost"]],
           fit_k1[["1"]][[1]][["ndpost"]],
           fit_k1[["0"]][[2]][["ndpost"]],
           fit_k1[["1"]][[2]][["ndpost"]],
           fit_k2[["0"]][["ndpost"]],
           fit_k2[["1"]][["ndpost"]])

  # Unique time-to-exit time labels
  time_seqs <- make_time_seqs(df = df)

  # List of minimum prediction times
  # Highest level grouping is who-is-being-predicted (Z) and subgrouping
  # is condition-under-which-prediction-is-happening (z)
  yrs <- list("0" = df$yr[df$z == 0],
              "1" = df$yr[df$z == 1])
  yts <- list("0" = df$yt[df$z == 0],
              "1" = df$yt[df$z == 1])
  dyrs <- list("0" = df$dyr[df$z == 0],
               "1" = df$dyr[df$z == 1])
  dyts <- list("0" = df$dyt[df$z == 0],
               "1" = df$dyt[df$z == 1])
  min_ttes <- list("0" = list("0" = make_min_ttes(yr = yrs[["0"]],
                                                  yt = yts[["0"]],
                                                  dyr = dyrs[["0"]],
                                                  dyt = dyts[["0"]],
                                                  factual = TRUE),
                              "1" = make_min_ttes(yr = yrs[["0"]],
                                                  yt = yts[["0"]],
                                                  dyr = dyrs[["0"]],
                                                  dyt = dyts[["0"]],
                                                  factual = FALSE)),
                   "1" = list("0" = make_min_ttes(yr = yrs[["1"]],
                                                  yt = yts[["1"]],
                                                  dyr = dyrs[["1"]],
                                                  dyt = dyts[["1"]],
                                                  factual = FALSE),
                              "1" = make_min_ttes(yr = yrs[["1"]],
                                                  yt = yts[["1"]],
                                                  dyr = dyrs[["1"]],
                                                  dyt = dyts[["1"]],
                                                  factual = TRUE)))

  # Start cluster
  cl <- parallel::makeCluster(parallel::detectCores())
  doParallel::registerDoParallel(cl)

  # Loop over MCMC iterations
  temp <- foreach::foreach(b = 1:B,
                           .combine = "acomb",
                           .multicombine = TRUE,
                           .packages = c("BART", "mstatebart")) %dopar% {
    # Do predictions for subgroup Z separately and bind together
    cbind(df$z, 1, rbind(ppredict_scr_b(b = b,
                                        x = x,
                                        time_seqs = time_seqs,
                                        Z = 0,
                                        min_ttes = min_ttes,
                                        yr = yrs[["0"]],
                                        yt = yts[["0"]],
                                        dyr = dyrs[["0"]],
                                        dyt = dyts[["0"]],
                                        fit_k1 = fit_k1,
                                        fit_k2 = fit_k2),
                         ppredict_scr_b(b = b,
                                        x = x,
                                        time_seqs = time_seqs,
                                        Z = 1,
                                        min_ttes = min_ttes,
                                        yr = yrs[["1"]],
                                        yt = yts[["1"]],
                                        dyr = dyrs[["1"]],
                                        dyt = dyts[["1"]],
                                        fit_k1 = fit_k1,
                                        fit_k2 = fit_k2)))
  }

  parallel::stopCluster(cl)

  # Rearrange to match the original sort order
  res <- temp * NA
  which_Z0 <- 1:NROW(x[["0"]])
  which_Z1 <- (NROW(x[["0"]]) + 1):NROW(df)
  res[which(df$z == 0), , ] <- temp[which_Z0, , ]
  res[which(df$z == 1), , ] <- temp[which_Z1, , ]
  dimnames(res)[[2]] <- c("z", "frailty",
                          "yr0_imp", "yt0_imp", "dyr0_imp", "dyt0_imp",
                          "yr1_imp", "yt1_imp", "dyr1_imp", "dyt1_imp")

  return(res)
}



#' Do all posterior predictions for a new data set using the bth MCMC iteration
#'
#' @param b MCMC iteration for prediction
#' @param x_new Matrix of N_test rows for which to make predictions
#' @param time_seqs Named list of unique time-to-event sequences for model
#' evaluation times
#' @param fit_k1 Named list of lists containing BART fits for time and disposition
#' of exit from Compartment 1
#' @param fit_k2 Named list of BART fits for exiting Compartment 2
#' @return N_test x 10 matrix of predicted times and observation indicators, where
#' N_test is the number of subjects with (factual) treatment equal to Z
#' @export
ppredict_scr_b_x_new <- function(b, x_new, time_seqs, fit_k1, fit_k2) {
  res <- matrix(NA, nrow = NROW(x_new), ncol = 8)
  zeros <- list(rep(0, NROW(x_new)), rep(0, NROW(x_new)))
  res[, 1:4] <- predict_scr_b(b = b, z = 0,
                              time_seqs_train = time_seqs[["0"]],
                              x_train = x_new, # design matrix used to fit tree
                              x_test = x_new, # design matrix needed for prediction
                              min_ttes = zeros,
                              fit_k1 = fit_k1[["0"]],
                              fit_k2 = fit_k2[["0"]])
  res[, 5:8] <- predict_scr_b(b = b, z = 1,
                              time_seqs_train = time_seqs[["1"]],
                              x_train = x_new, # design matrix used to fit tree
                              x_test = x_new, # design matrix needed for prediction
                              min_ttes = zeros,
                              fit_k1 = fit_k1[["1"]],
                              fit_k2 = fit_k2[["1"]])

  return(res)
}



#' Vectorized posterior prediction for new data
#'
#' See \code{\link{pp_scr_xnew}} for version that gives array output
v_ppredict_scr_b_x_new <- Vectorize(ppredict_scr_b_x_new,
                                    vectorize.args = "b")



#' Do BART (non-recurrent) semicompeting risks posterior predictive draws for
#' new observations, all z levels, and all MCMC iterations
#'
#' @param x_new Matrix of N_test rows for which to make predictions
#' @param time_seqs Named list of unique time-to-event sequences for model
#' evaluation times from \code{\link{make_time_seqs}}
#' @param fit_k1 List of lists containing BART fits for Compartment 1 stays
#' (1 = time-to-exit model, 2 = cause-of-exit model)
#' @param fit_k2 List of BART fits for Compartment 2 stays
#' @return N x 10 x B posterior prediction array
#' @export
# posterior_predict_scr_bart_x_new <- function(x_new, time_seqs, fit_k1, fit_k2) {
#
#   # Need all models to get predictions
#   # Therefore can only predict for the minimum number of post-warmup draws
#   B <- min(fit_k1[["0"]][[1]][["ndpost"]],
#            fit_k1[["1"]][[1]][["ndpost"]],
#            fit_k1[["0"]][[2]][["ndpost"]],
#            fit_k1[["1"]][[2]][["ndpost"]],
#            fit_k2[["0"]][["ndpost"]],
#            fit_k2[["1"]][["ndpost"]])
#
#   # Start cluster
#   cl <- parallel::makeCluster(parallel::detectCores())
#   doParallel::registerDoParallel(cl)
#
#   # Loop over MCMC iterations
#   res <- array(NA, dim = c(nrow(x_new), 10, B),
#                dimnames = list(NULL,
#                                c("z", "frailty",
#                                  "yr0_imp", "yt0_imp", "dyr0_imp", "dyt0_imp",
#                                  "yr1_imp", "yt1_imp", "dyr1_imp", "dyt1_imp"),
#                                NULL))
#   res[, 3:10, ] <- foreach::foreach(b = 1:B,
#                                     # .combine = "acomb",
#                                     # .multicombine = TRUE,
#                                     combine = "c",
#                                     .packages = c("BART", "mstatebart")) %dopar% {
#
#     ppredict_scr_b_x_new(b = b, x = x_new, time_seqs, fit_k1, fit_k2)
#   }
#
#   parallel::stopCluster(cl)
#   return(res)
#
# }



#' Do BART (non-recurrent) semicompeting risks posterior predictive draws for
#' new observations, all z levels, and all MCMC iterations
#'
#' @param x_new Matrix of N_test rows for which to make predictions
#' @param time_seqs Named list of unique time-to-event sequences for model
#' evaluation times from \code{\link{make_time_seqs}}
#' @param fit_k1 List of lists containing BART fits for Compartment 1 stays
#' (1 = time-to-exit model, 2 = cause-of-exit model)
#' @param fit_k2 List of BART fits for Compartment 2 stays
#' @return N x 10 x B posterior prediction array
#' @export
pp_scr_xnew <- function(x_new, time_seqs, fit_k1, fit_k2) {

  # Need all models to get predictions
  B <- 20
  # Therefore can only predict for the minimum number of post-warmup draws
  # B <- min(fit_k1[["0"]][[1]][["ndpost"]],
  #          fit_k1[["1"]][[1]][["ndpost"]],
  #          fit_k1[["0"]][[2]][["ndpost"]],
  #          fit_k1[["1"]][[2]][["ndpost"]],
  #          fit_k2[["0"]][["ndpost"]],
  #          fit_k2[["1"]][["ndpost"]])

  nrep <- nrow(x_new)

  # Loop over MCMC iterations
  res <- array(NA, dim = c(nrep, 10, B),
               dimnames = list(NULL,
                               c("z", "frailty",
                                 "yr0_imp", "yt0_imp", "dyr0_imp", "dyt0_imp",
                                 "yr1_imp", "yt1_imp", "dyr1_imp", "dyt1_imp"),
                               NULL))

  res0 <- v_ppredict_scr_b_x_new(b = 1:B, x = x_new,
                                 time_seqs, fit_k1, fit_k2)

  for (b in 1:B) {
    res[ , 3:10, b] <- matrix(res0[, b], nrow = nrep)
  }

  return(res)

}
