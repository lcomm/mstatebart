#' Plot predicted vs. observed event times for a survival BART
#'
#' @param surv_bart_fit Probit-survival fit from
#' @param input_list List of inputs used to construct
#' @return ggplot2 plot object
#' @export
plot_surv_pred <- function(surv_bart_fit, input_list) {
  stopifnot(c("times", "x_train") %in% names(input_list))
  old_times <- input_list[["times"]]
  time_seq <- sort(unique(old_times))
  x_new <- input_list[["x_train"]]
  pred <- predict_times_from_pbart(pbart_fit = surv_bart_fit,
                                   time_seq = time_seq,
                                   x_new = x_new)
  new_times <- c(pred)

  new_times[is.infinite(new_times)] <- max(time_seq)
  to_plot <- data.frame(Type = c(rep("Predicted", length(new_times)),
                                 rep("Observed", length(old_times))),
                        Times = c(new_times,
                                  old_times))
  p <- ggplot(to_plot, aes_string(y = "Times",
                                  fill = "Type")) +
    geom_density(alpha = 0.2)
  return(p)
}

