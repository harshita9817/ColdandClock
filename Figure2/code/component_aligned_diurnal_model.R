# Fit the fixed external diurnal prediction to the DAWN-minus-DUSK contrast.
#
# The prediction is not applied to the response shared between the two onsets.
# It carries no gene-specific information there: cross-gene agreement with the
# DAWN/DUSK mean stays below r^2 = 0.07 at every duration and is 0.001 at 11 h,
# and a term fitted to that channel is indistinguishable from the gene-specific
# shared cold response added immediately afterwards, which absorbs it with no
# change in total explained variance. Scales use only the two supplied training
# trials.

fit_one_diurnal_component <- function(response, fixed_prediction, training_trials) {
  if (!is.matrix(response) || ncol(response) != 3L ||
      length(training_trials) != 2L || anyDuplicated(training_trials) ||
      !all(training_trials %in% 1:3)) {
    stop("A component fit requires a gene-by-three-trial matrix and two training trials.")
  }
  if (length(fixed_prediction) != nrow(response) ||
      any(!is.finite(response)) || any(!is.finite(fixed_prediction))) {
    stop("Non-finite or misaligned component data.")
  }
  denominator <- sum(fixed_prediction^2)
  if (!is.finite(denominator) || denominator <= 0) {
    return(list(scale = 0, prediction = numeric(length(fixed_prediction))))
  }

  scale <- sum(
    fixed_prediction * rowSums(response[, training_trials, drop = FALSE])
  ) / (length(training_trials) * denominator)
  list(scale = scale, prediction = scale * fixed_prediction)
}

fit_component_aligned_diurnal_model <- function(
    dawn_response, dusk_response, fixed_dawn, fixed_dusk, training_trials) {
  if (!identical(dim(dawn_response), dim(dusk_response))) {
    stop("DAWN and DUSK response matrices must have identical dimensions.")
  }
  contrast_response <- (dawn_response - dusk_response) / 2
  fixed_contrast <- (fixed_dawn - fixed_dusk) / 2
  contrast_fit <- fit_one_diurnal_component(
    contrast_response, fixed_contrast, training_trials
  )
  list(
    contrast = contrast_fit,
    predicted_dawn = contrast_fit$prediction,
    predicted_dusk = -contrast_fit$prediction
  )
}
