#!/usr/bin/env Rscript

# Generate the canonical gene-level table for early (3 h) and late (27 h)
# responses. The external diurnal waveform is removed first. A single joint
# limma model then uses the crossed 3/15-h cells to separate collection phase in
# the cold (ZT4/light versus ZT16/dark) from a distinct 27-h onset-history term.
# Shared, DAWN/DUSK-specific, and light/dark-associated responses are tested as
# contrasts from that fit so uncertainty in the phase estimate reaches the final
# P values; no point-estimate phase subtraction or 27-h validation rule is used
# for statistical selection.

suppressPackageStartupMessages(library(limma))
if (as.character(packageVersion("limma")) != "3.68.4") {
  stop("Pinned limma version is 3.68.4.")
}

script_argument <- grep("^--file=", commandArgs(FALSE), value = TRUE)
# Rscript writes spaces in --file as "~+~"; the gsub() undoes that so paths with spaces work.
root <- dirname(dirname(normalizePath(
  gsub("~+~", " ", sub("^--file=", "", script_argument), fixed = TRUE), mustWork = TRUE
)))
source(file.path(root, "code", "component_aligned_diurnal_model.R"))
source(file.path(root, "code", "analysis_universe.R"))
reported_times <- c(3, 27)
phase_training_times <- c(3, 15)
effect_cutoff <- 1
alpha <- 0.05
onset_categories <- c(
  "shared_induction", "dawn_specific_induction", "dusk_specific_induction"
)

normalization <- if (any(grepl("^--normalization=size-factor$", commandArgs(FALSE)))) {
  "size_factor"
} else {
  "tpm"
}
universe <- load_analysis_universe(root, normalization = normalization)
genes <- universe$genes
metadata <- universe$metadata
reference <- universe$reference
current_log <- universe$current_log

directional_p <- function(p_value, coefficient, direction) {
  if (direction == "positive") {
    ifelse(coefficient > 0, p_value / 2, 1)
  } else {
    ifelse(coefficient < 0, p_value / 2, 1)
  }
}

# Apply the same conjunction logic to coefficients and P values from either a
# simple two-condition fit or the joint duration/phase/onset model below.
classify_two_condition_categories <- function(
    coefficient, p_value, first, second,
    first_specific_name, second_specific_name, include_shared = TRUE) {
  conjunction_p <- cbind(
    shared_induction = pmax(
      directional_p(p_value[, "FIRST"], coefficient[, "FIRST"], "positive"),
      directional_p(p_value[, "SECOND"], coefficient[, "SECOND"], "positive")
    ),
    first_specific = pmax(
      directional_p(p_value[, "FIRST"], coefficient[, "FIRST"], "positive"),
      directional_p(
        p_value[, "FIRST_minus_SECOND"],
        coefficient[, "FIRST_minus_SECOND"], "positive"
      )
    ),
    second_specific = pmax(
      directional_p(p_value[, "SECOND"], coefficient[, "SECOND"], "positive"),
      directional_p(
        p_value[, "FIRST_minus_SECOND"],
        coefficient[, "FIRST_minus_SECOND"], "negative"
      )
    )
  )
  colnames(conjunction_p) <- c(
    "shared_induction", first_specific_name, second_specific_name
  )
  category_names <- colnames(conjunction_p)
  if (!include_shared) {
    category_names <- setdiff(category_names, "shared_induction")
    conjunction_p <- conjunction_p[, category_names, drop = FALSE]
  }
  conjunction_q <- apply(conjunction_p, 2, p.adjust, method = "BH")
  sign_consistent <- cbind(
    shared_induction = apply(first > 0, 1, all) & apply(second > 0, 1, all),
    first_specific = apply(first > 0, 1, all) &
      apply(second <= 0, 1, all) & apply(first - second > 0, 1, all),
    second_specific = apply(second > 0, 1, all) &
      apply(first <= 0, 1, all) & apply(first - second < 0, 1, all)
  )
  colnames(sign_consistent) <- c(
    "shared_induction", first_specific_name, second_specific_name
  )
  sign_consistent <- sign_consistent[, category_names, drop = FALSE]
  effect_pass <- cbind(
    shared_induction = coefficient[, "FIRST"] >= effect_cutoff &
      coefficient[, "SECOND"] >= effect_cutoff,
    first_specific = coefficient[, "FIRST"] >= effect_cutoff &
      coefficient[, "SECOND"] <= 0 &
      coefficient[, "FIRST_minus_SECOND"] >= effect_cutoff,
    second_specific = coefficient[, "SECOND"] >= effect_cutoff &
      coefficient[, "FIRST"] <= 0 &
      coefficient[, "FIRST_minus_SECOND"] <= -effect_cutoff
  )
  colnames(effect_pass) <- c(
    "shared_induction", first_specific_name, second_specific_name
  )
  effect_pass <- effect_pass[, category_names, drop = FALSE]
  calls <- conjunction_q <= alpha & effect_pass & sign_consistent
  if (include_shared) {
    # Statistical evidence is the BH-controlled conjunction call. Effect size
    # and trial-level sign consistency are biologically distinct from FDR
    # control. Retain the one-log2 threshold for biological relevance and report
    # trial-level sign consistency as an audit rather than another selection
    # rule.
    calls[, "shared_induction"] <-
      conjunction_q[, "shared_induction"] <= alpha &
      effect_pass[, "shared_induction"]
  }
  if (any(rowSums(calls) > 1L)) {
    stop("Response categories are not mutually exclusive.")
  }
  category <- rep("not_called", nrow(first))
  for (name in category_names) {
    category[calls[, name]] <- name
  }
  list(
    category = category, coefficient = coefficient,
    p_value = conjunction_p, q_value = conjunction_q,
    sign_consistent = sign_consistent
  )
}

call_two_condition_categories <- function(
    first, second, first_specific_name, second_specific_name,
    include_shared = TRUE) {
  response <- cbind(first, second)
  samples <- data.frame(
    condition = factor(rep(c("FIRST", "SECOND"), each = 3),
      levels = c("FIRST", "SECOND")),
    trial = factor(rep(1:3, 2))
  )
  design <- model.matrix(
    ~ 0 + condition + trial, samples,
    contrasts.arg = list(trial = contr.sum(3))
  )
  contrasts <- makeContrasts(
    FIRST = conditionFIRST,
    SECOND = conditionSECOND,
    FIRST_minus_SECOND = conditionFIRST - conditionSECOND,
    levels = design
  )
  fit <- eBayes(contrasts.fit(lmFit(response, design), contrasts))
  classify_two_condition_categories(
    fit$coefficients, fit$p.value, first, second,
    first_specific_name, second_specific_name, include_shared
  )
}

# For every duration and omitted trial, fit the fixed external prediction to the
# DAWN-minus-DUSK contrast using only the other two trials.
build_diurnal_adjusted_responses <- function(elapsed) {
  observed <- fixed <- list()
  for (onset in c("DAWN", "DUSK")) {
    time <- onset_time_h(onset, elapsed)
    observed[[onset]] <- response_matrix(universe, onset, time)
    fixed[[onset]] <- fixed_onset_hold(universe, onset, time)
  }

  dawn <- dusk <- matrix(
    NA_real_, length(genes), 3,
    dimnames = list(genes, paste0("trial", 1:3))
  )
  predicted_dawn <- predicted_dusk <- dawn
  contrast_scale <- numeric(3)
  for (held in 1:3) {
    training <- setdiff(1:3, held)
    model <- fit_component_aligned_diurnal_model(
      observed$DAWN, observed$DUSK, fixed$DAWN, fixed$DUSK, training
    )
    predicted_dawn[, held] <- model$predicted_dawn
    predicted_dusk[, held] <- model$predicted_dusk
    dawn[, held] <- observed$DAWN[, held] - model$predicted_dawn
    dusk[, held] <- observed$DUSK[, held] - model$predicted_dusk
    contrast_scale[held] <- model$contrast$scale
  }
  list(
    dawn = dawn, dusk = dusk, observed = observed, fixed = fixed,
    predicted_dawn = predicted_dawn, predicted_dusk = predicted_dusk,
    contrast_scale = contrast_scale
  )
}

all_times <- sort(unique(c(reported_times, phase_training_times)))
responses <- setNames(lapply(all_times, build_diurnal_adjusted_responses), all_times)

# The crossed 3/15-h cells identify one collection-phase coefficient, while a
# separate 27-h onset coefficient captures late onset history. This is a single
# constrained fit: the DAWN-minus-DUSK onset effect is assumed to be the same at
# 3 and 15 h, and the collection-phase coefficient is assumed to carry across
# durations. Growth trial is blocked separately within each duration.
light_by_trial <- 0.5 * (
  responses[["3"]]$dawn + responses[["15"]]$dusk
)
dark_by_trial <- 0.5 * (
  responses[["3"]]$dusk + responses[["15"]]$dawn
)
joint_response <- do.call(cbind, lapply(c(3, 15, 27), function(elapsed) {
  cbind(responses[[as.character(elapsed)]]$dawn,
        responses[[as.character(elapsed)]]$dusk)
}))
joint_samples <- do.call(rbind, lapply(c(3, 15, 27), function(elapsed) {
  data.frame(
    duration = elapsed,
    onset = rep(c("DAWN", "DUSK"), each = 3),
    trial = rep(1:3, 2), stringsAsFactors = FALSE
  )
}))
joint_samples$phase <- ifelse(
  (joint_samples$onset == "DAWN" & joint_samples$duration %in% c(3, 27)) |
    (joint_samples$onset == "DUSK" & joint_samples$duration == 15),
  "LIGHT_ZT4", "DARK_ZT16"
)
joint_samples$collection_phase <-
  ifelse(joint_samples$phase == "LIGHT_ZT4", 0.5, -0.5)
joint_samples$early_onset_difference <- ifelse(
  joint_samples$duration %in% c(3, 15),
  ifelse(joint_samples$onset == "DAWN", 0.5, -0.5), 0
)
joint_samples$late_onset_difference <- ifelse(
  joint_samples$duration == 27,
  ifelse(joint_samples$onset == "DAWN", 0.5, -0.5), 0
)
joint_samples$duration <- factor(joint_samples$duration, levels = c(3, 15, 27))
joint_samples$trial <- factor(joint_samples$trial)
joint_design <- model.matrix(
  ~ 0 + duration + collection_phase + early_onset_difference +
    late_onset_difference + duration:trial,
  joint_samples, contrasts.arg = list(trial = contr.sum(3))
)
if (qr(joint_design)$rank != ncol(joint_design)) {
  stop("The constrained duration/phase/onset design is not full rank.")
}

joint_contrasts <- matrix(
  0, nrow = ncol(joint_design), ncol = 6,
  dimnames = list(colnames(joint_design), c(
    "DAWN27", "DUSK27", "DAWN_minus_DUSK27",
    "LIGHT_early", "DARK_early", "LIGHT_minus_DARK"
  ))
)
joint_contrasts["duration27", c("DAWN27", "DUSK27")] <- 1
joint_contrasts["late_onset_difference", "DAWN27"] <- 0.5
joint_contrasts["late_onset_difference", "DUSK27"] <- -0.5
joint_contrasts["late_onset_difference", "DAWN_minus_DUSK27"] <- 1
joint_contrasts[c("duration3", "duration15"),
                c("LIGHT_early", "DARK_early")] <- 0.5
joint_contrasts["collection_phase", "LIGHT_early"] <- 0.5
joint_contrasts["collection_phase", "DARK_early"] <- -0.5
joint_contrasts["collection_phase", "LIGHT_minus_DARK"] <- 1

joint_fit <- eBayes(contrasts.fit(
  lmFit(joint_response, joint_design), joint_contrasts
))
if (!all(joint_fit$df.residual == 6)) {
  stop("Expected six residual degrees of freedom in the joint model.")
}
direct_phase_effect <- rowMeans(light_by_trial - dark_by_trial)
direct_27_dawn <- rowMeans(responses[["27"]]$dawn) - direct_phase_effect / 2
direct_27_dusk <- rowMeans(responses[["27"]]$dusk) + direct_phase_effect / 2
if (max(abs(joint_fit$coefficients[, "LIGHT_minus_DARK"] -
            direct_phase_effect)) > 1e-10 ||
    max(abs(joint_fit$coefficients[, "DAWN27"] - direct_27_dawn)) > 1e-10 ||
    max(abs(joint_fit$coefficients[, "DUSK27"] - direct_27_dusk)) > 1e-10) {
  stop("Joint-model contrasts do not reproduce the intended balanced estimates.")
}
late_coefficient <- joint_fit$coefficients[, c(
  "DAWN27", "DUSK27", "DAWN_minus_DUSK27"
), drop = FALSE]
late_p_value <- joint_fit$p.value[, c(
  "DAWN27", "DUSK27", "DAWN_minus_DUSK27"
), drop = FALSE]
colnames(late_coefficient) <- colnames(late_p_value) <-
  c("FIRST", "SECOND", "FIRST_minus_SECOND")
phase_effect <- joint_fit$coefficients[, "LIGHT_minus_DARK"]
controlled_27_dawn <- sweep(responses[["27"]]$dawn, 1, phase_effect / 2, "-")
controlled_27_dusk <- sweep(responses[["27"]]$dusk, 1, phase_effect / 2, "+")
late_calls <- classify_two_condition_categories(
  late_coefficient, late_p_value, controlled_27_dawn, controlled_27_dusk,
  "dawn_specific_induction", "dusk_specific_induction"
)

phase_coefficient <- joint_fit$coefficients[, c(
  "LIGHT_early", "DARK_early", "LIGHT_minus_DARK"
), drop = FALSE]
phase_p_value <- joint_fit$p.value[, c(
  "LIGHT_early", "DARK_early", "LIGHT_minus_DARK"
), drop = FALSE]
colnames(phase_coefficient) <- colnames(phase_p_value) <-
  c("FIRST", "SECOND", "FIRST_minus_SECOND")
phase_calls <- classify_two_condition_categories(
  phase_coefficient, phase_p_value, light_by_trial, dark_by_trial,
  "light_cold_specific_induction", "dark_cold_specific_induction",
  include_shared = FALSE
)
collection_phase_category <- phase_calls$category

output <- list()
for (elapsed in reported_times) {
  original <- responses[[as.character(elapsed)]]
  if (elapsed == 27) {
    analyzed_dawn <- controlled_27_dawn
    analyzed_dusk <- controlled_27_dusk
    response_basis <- "external_diurnal_and_joint_phase_onset_model"
    phase_category <- collection_phase_category
    calls <- late_calls
  } else {
    analyzed_dawn <- original$dawn
    analyzed_dusk <- original$dusk
    response_basis <- "external_diurnal_adjusted"
    phase_category <- rep("not_evaluated", length(genes))
    calls <- call_two_condition_categories(
      analyzed_dawn, analyzed_dusk,
      "dawn_specific_induction", "dusk_specific_induction"
    )
  }
  called <- calls$category
  category_q <- rep(NA_real_, length(genes))
  for (name in onset_categories) {
    selected <- called == name
    category_q[selected] <- calls$q_value[selected, name]
  }
  reported <- ifelse(called == "not_called", "", called)
  if (elapsed == 27) {
    reported <- ifelse(
      phase_category == "not_called", reported,
      ifelse(reported == "", phase_category, paste(reported, phase_category, sep = ";"))
    )
  }
  reported[reported == ""] <- "not_called"

  one <- data.frame(
    gene_id = genes,
    hours_after_onset = elapsed,
    response_basis = response_basis,
    mean_external_tpm = reference$mean_external_tpm,
    jtk_adjusted_p = reference$jtk_adjusted_p,
    rhythmic_adjusted_p_le_0_01 = reference$rhythmic_adjusted_p_le_0_01,
    fixed_diurnal_prediction_dawn = original$fixed$DAWN,
    fixed_diurnal_prediction_dusk = original$fixed$DUSK,
    mean_crossfit_diurnal_contrast_scale = mean(original$contrast_scale),
    mean_adjusted_dawn_response_log2 = rowMeans(analyzed_dawn),
    mean_adjusted_dusk_response_log2 = rowMeans(analyzed_dusk),
    dawn_minus_dusk_log2 = rowMeans(analyzed_dawn - analyzed_dusk),
    called_category = called,
    collection_phase_specific_category = phase_category,
    reported_categories = reported,
    category_q_value = category_q,
    limma_version = as.character(packageVersion("limma")),
    stringsAsFactors = FALSE
  )
  for (name in onset_categories) {
    one[[paste0(name, "_p_value")]] <- calls$p_value[, name]
    one[[paste0(name, "_q_value")]] <- calls$q_value[, name]
    one[[paste0(name, "_sign_consistent")]] <- calls$sign_consistent[, name]
  }
  for (trial in 1:3) {
    one[[paste0("trial_", trial, "_adjusted_dawn_response_log2")]] <-
      analyzed_dawn[, trial]
    one[[paste0("trial_", trial, "_adjusted_dusk_response_log2")]] <-
      analyzed_dusk[, trial]
    one[[paste0("trial_", trial, "_diurnal_only_adjusted_dawn_response_log2")]] <-
      original$dawn[, trial]
    one[[paste0("trial_", trial, "_diurnal_only_adjusted_dusk_response_log2")]] <-
      original$dusk[, trial]
  }
  output[[as.character(elapsed)]] <- one
}

result <- do.call(rbind, output)
rownames(result) <- NULL
late_rows <- result$hours_after_onset == 27
late_index <- match(result$gene_id[late_rows], genes)
result$collection_phase_effect_log2 <- NA_real_
result$mean_early_light_zt4_response_log2 <- NA_real_
result$mean_early_dark_zt16_response_log2 <- NA_real_
result$collection_phase_effect_log2[late_rows] <- phase_effect[late_index]
result$mean_early_light_zt4_response_log2[late_rows] <-
  phase_coefficient[late_index, "FIRST"]
result$mean_early_dark_zt16_response_log2[late_rows] <-
  phase_coefficient[late_index, "SECOND"]
phase_category_q <- rep(NA_real_, length(genes))
for (name in c("light_cold_specific_induction", "dark_cold_specific_induction")) {
  selected <- collection_phase_category == name
  phase_category_q[selected] <- phase_calls$q_value[selected, name]
  result[[paste0(name, "_p_value")]] <- NA_real_
  result[[paste0(name, "_q_value")]] <- NA_real_
  result[[paste0(name, "_sign_consistent")]] <- NA
  result[[paste0(name, "_p_value")]][late_rows] <-
    phase_calls$p_value[late_index, name]
  result[[paste0(name, "_q_value")]][late_rows] <-
    phase_calls$q_value[late_index, name]
  result[[paste0(name, "_sign_consistent")]][late_rows] <-
    phase_calls$sign_consistent[late_index, name]
}
result$collection_phase_category_q_value <- NA_real_
result$collection_phase_category_q_value[late_rows] <- phase_category_q[late_index]

path <- file.path(root, "data/generatable", sprintf(
  "diurnal_adjusted_response_categories%s.csv.gz",
  if (normalization == "size_factor") "_size_factor" else ""
))
connection <- gzfile(path, "wt", compression = 9)
write.csv(result, connection, row.names = FALSE, quote = TRUE)
close(connection)

cat("Onset/shared categories after all applicable controls:\n")
print(with(result, table(hours_after_onset, called_category)))
cat("\nJoint-model collection-phase-specific induction:\n")
print(table(result$collection_phase_specific_category[late_rows]))
cat("Wrote", nrow(result), "rows to", path, "\n")
