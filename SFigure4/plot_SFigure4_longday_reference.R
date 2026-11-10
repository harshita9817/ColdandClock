#!/usr/bin/env Rscript

# Fig. S4: the Figure 3 analyses repeated with the long-day (14L:10D) JGI BTx623
# diel series as the external diurnal reference. This is figures/Figure3/
# plot_figure3.R with three changes: Dawn/Dusk legend labels, a panel E of
# response-category counts, and SFigure4 output names. It must be run by
# run_SFigure4_longday_reference.sh, which places it inside a temporary copy of
# the repository whose external reference has been rebuilt from the long-day
# series.
#
# Original Figure 3 description: combines the fixed external diurnal-arrest
# test, biological examples, and held-out comparison of three response models. The analysis
# universe is the intersection of genes averaging at least 3 TPM in the
# independent series and genes tested by DESeq2; JTK_CYCLE significance is not
# an inclusion criterion.

script_argument <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_argument) != 1L) {
  stop("Run this figure script with Rscript so repository-relative paths can be resolved.")
}
script_path <- normalizePath(sub("^--file=", "", script_argument), mustWork = TRUE)
figure_dir <- dirname(script_path)
analysis_dir <- dirname(dirname(figure_dir))
source(file.path(analysis_dir, "scripts", "component_aligned_diurnal_model.R"))
source(file.path(analysis_dir, "scripts", "analysis_universe.R"))
durations <- c(3, 7, 11, 15, 19, 23, 27)

universe <- load_analysis_universe(analysis_dir)
genes <- universe$genes
metadata <- universe$metadata
current_log <- universe$current_log

observed <- prediction <- list()
gene_rows <- contrast_rows <- list()
for (elapsed in durations) {
  for (treatment in c("DAWN", "DUSK")) {
    time_h <- onset_time_h(treatment, elapsed)
    key <- paste(elapsed, treatment, sep = "_")
    observed[[key]] <- response_matrix(universe, treatment, time_h)
    prediction[[key]] <- fixed_onset_hold(universe, treatment, time_h)
    mean_observed <- rowMeans(observed[[key]])
    gene_rows[[length(gene_rows) + 1L]] <- data.frame(
      gene_id = genes, treatment = treatment, hours_after_onset = elapsed,
      observed_cold_minus_control_log2_tpm = mean_observed,
      fixed_onset_hold_prediction = prediction[[key]], stringsAsFactors = FALSE
    )
    correlation <- cor(mean_observed, prediction[[key]])
    contrast_rows[[length(contrast_rows) + 1L]] <- data.frame(
      treatment = treatment, hours_after_onset = elapsed, genes = length(genes),
      pearson_r = correlation,
      variance_aligned_r_squared = correlation^2, stringsAsFactors = FALSE
    )
  }
}
gene_data <- do.call(rbind, gene_rows)
contrast_data <- do.call(rbind, contrast_rows)

# The fixed prediction is applied only to the DAWN-minus-DUSK contrast. Report
# agreement in both channels so that choice is visible: per-onset correlations
# blend them and hide how little the prediction says about the shared response.
channel_rows <- list()
for (elapsed in durations) {
  dawn <- rowMeans(observed[[paste(elapsed, "DAWN", sep = "_")]])
  dusk <- rowMeans(observed[[paste(elapsed, "DUSK", sep = "_")]])
  fixed_dawn <- prediction[[paste(elapsed, "DAWN", sep = "_")]]
  fixed_dusk <- prediction[[paste(elapsed, "DUSK", sep = "_")]]
  channel_rows[[length(channel_rows) + 1L]] <- data.frame(
    hours_after_onset = elapsed,
    mean_channel_r_squared = cor((dawn + dusk) / 2, (fixed_dawn + fixed_dusk) / 2)^2,
    contrast_channel_r_squared = cor((dawn - dusk) / 2, (fixed_dawn - fixed_dusk) / 2)^2,
    stringsAsFactors = FALSE
  )
}
channel_data <- do.call(rbind, channel_rows)
dawn_11 <- gene_data[
  gene_data$treatment == "DAWN" & gene_data$hours_after_onset == 11,
]
dusk_11 <- gene_data[
  gene_data$treatment == "DUSK" & gene_data$hours_after_onset == 11,
]
dusk_11 <- dusk_11[match(dawn_11$gene_id, dusk_11$gene_id), ]
matched_11h_correlation <- cor(
  dawn_11$observed_cold_minus_control_log2_tpm,
  dusk_11$observed_cold_minus_control_log2_tpm
)

# Sequential prediction is trained on two growth trials and scored on the third.
fold_rows <- list()
for (elapsed in durations) {
  dawn <- observed[[paste(elapsed, "DAWN", sep = "_")]]
  dusk <- observed[[paste(elapsed, "DUSK", sep = "_")]]
  dawn_prediction <- prediction[[paste(elapsed, "DAWN", sep = "_")]]
  dusk_prediction <- prediction[[paste(elapsed, "DUSK", sep = "_")]]
  for (held_trial in 1:3) {
    training <- setdiff(1:3, held_trial)
    diurnal_model <- fit_component_aligned_diurnal_model(
      dawn, dusk, dawn_prediction, dusk_prediction, training
    )
    dawn_residual <- sweep(
      dawn[, training, drop = FALSE], 1, diurnal_model$predicted_dawn, "-"
    )
    dusk_residual <- sweep(
      dusk[, training, drop = FALSE], 1, diurnal_model$predicted_dusk, "-"
    )
    shared <- (rowSums(dawn_residual) + rowSums(dusk_residual)) / 4
    specific <- (rowMeans(dawn_residual) - rowMeans(dusk_residual)) / 2
    held_observed <- c(dawn[, held_trial], dusk[, held_trial])
    pred1 <- c(diurnal_model$predicted_dawn, diurnal_model$predicted_dusk)
    pred2 <- pred1 + rep(shared, 2)
    pred3 <- pred2 + c(specific, -specific)
    fold_rows[[length(fold_rows) + 1L]] <- data.frame(
      hours_after_onset = elapsed, held_out_trial = held_trial,
      zero_sse = sum(held_observed^2),
      diurnal_sse = sum((held_observed - pred1)^2),
      shared_sse = sum((held_observed - pred2)^2),
      full_sse = sum((held_observed - pred3)^2),
      diurnal_contrast_scale = diurnal_model$contrast$scale,
      stringsAsFactors = FALSE
    )
  }
}
fold_data <- do.call(rbind, fold_rows)
held_out <- do.call(rbind, lapply(split(fold_data, fold_data$hours_after_onset), function(one) {
  sse <- colSums(one[, c("zero_sse", "diurnal_sse", "shared_sse", "full_sse")])
  data.frame(
    hours_after_onset = unique(one$hours_after_onset),
    diurnal_increment = (sse[1] - sse[2]) / sse[1],
    shared_increment = (sse[2] - sse[3]) / sse[1],
    onset_increment = (sse[3] - sse[4]) / sse[1],
    held_out_residual = sse[4] / sse[1], stringsAsFactors = FALSE
  )
}))
held_out <- held_out[order(held_out$hours_after_onset), ]

colors <- c(DAWN = "#D55E00", DUSK = "#0072B2")
component_colors <- c(
  "Shared response" = "#1B7837",
  "Specific response" = "#D95F0E",
  "Diurnal" = "#756BB1",
  "Residual" = "#BDBDBD"
)
selected_go <- c(
  "GO:0009765" = "photosynthesis, light harvesting",
  "GO:0042254" = "ribosome biogenesis",
  "GO:0070417" = "cellular response to cold"
)
go_sets <- read.delim(
  gzfile(file.path(analysis_dir, "data", "generatable", "go_tested_sets.tsv.gz")),
  check.names = FALSE, stringsAsFactors = FALSE
)
go_rows <- list()
for (go_id in names(selected_go)) {
  members <- strsplit(go_sets$gene_ids[go_sets$go_id == go_id], ";", fixed = TRUE)[[1]]
  members <- intersect(members, genes)
  for (treatment in names(colors)) {
    for (elapsed in durations) {
      one <- gene_data[
        gene_data$gene_id %in% members & gene_data$treatment == treatment &
          gene_data$hours_after_onset == elapsed,
      ]
      go_rows[[length(go_rows) + 1L]] <- data.frame(
        go_id = go_id, go_name = unname(selected_go[go_id]), treatment = treatment,
        hours_after_onset = elapsed, genes = nrow(one),
        mean_observed_change = mean(one$observed_cold_minus_control_log2_tpm),
        mean_fixed_prediction = mean(one$fixed_onset_hold_prediction),
        stringsAsFactors = FALSE
      )
    }
  }
}
go_data <- do.call(rbind, go_rows)
go_limits <- range(c(
  go_data$mean_observed_change[go_data$go_id %in% names(selected_go)],
  go_data$mean_fixed_prediction[go_data$go_id %in% names(selected_go)], 0
), finite = TRUE)
go_padding <- diff(go_limits) * 0.05
go_limits <- go_limits + c(-go_padding, go_padding)

if (!all(durations %in% contrast_data$hours_after_onset) ||
    !all(durations %in% held_out$hours_after_onset) ||
    !all(names(selected_go) %in% go_data$go_id)) {
  stop("A required duration or GO category is missing from the saved results.")
}

panel_label <- function(label) {
  # The label is placed in the y-axis-title column at the exact height of the
  # top of the y-axis line. mtext positions text by its midpoint, so subtract
  # half its user-coordinate height to align the top edge instead.
  top_aligned_position <- par("usr")[4] - strheight(label, cex = 1) / 2
  mtext(
    label, side = 2, line = 2.75, at = top_aligned_position, las = 1,
    font = 2, cex = 1, adj = 0.5
  )
}

draw_scatter <- function(treatment, panel) {
  one <- gene_data[
    gene_data$treatment == treatment & gene_data$hours_after_onset == 11,
  ]
  x <- one$fixed_onset_hold_prediction
  y <- one$observed_cold_minus_control_log2_tpm
  correlation <- cor(x, y)
  limits_x <- range(x, finite = TRUE)
  limits_y <- range(y, finite = TRUE)
  plot(
    x, y, pch = 16, cex = 0.32,
    col = adjustcolor(colors[[treatment]], alpha.f = 0.12),
    xlab = "",
    ylab = "",
    xlim = limits_x, ylim = limits_y, main = ""
  )
  y_range <- par("usr")[4] - par("usr")[3]
  mtext(
    expression("Observed response (" * log[2] * ")"), side = 2,
    line = 2.05, at = mean(par("usr")[3:4]) - 0.07 * y_range
  )
  abline(lm(y ~ x), col = "#222222", lwd = 1.4)
  mtext("Cold response predicted\nfrom diurnal cycling", side = 1, line = 2.65)
  legend(
    "topleft", legend = sprintf("r = %.2f", correlation),
    text.col = "#222222", bty = "n", cex = 1
  )
  if (nzchar(panel)) panel_label(panel)
}

draw_agreement <- function() {
  plot(
    NA, xlim = range(durations), ylim = c(0, 0.65), xaxt = "n",
    xlab = "Hours after cold onset",
    ylab = expression("Variance explained (" * R^2 * ")"), main = ""
  )
  axis(1, at = durations)
  abline(h = 0, col = "#999999", lty = 2)
  for (treatment in names(colors)) {
    one <- contrast_data[contrast_data$treatment == treatment, ]
    one <- one[order(one$hours_after_onset), ]
    lines(
      one$hours_after_onset, one$variance_aligned_r_squared,
      type = "o", pch = 16, lwd = 1.7, col = colors[[treatment]]
    )
  }
  legend(
    "topright", legend = c("Dawn", "Dusk"), col = colors, pch = 16,
    lwd = 1.7, bty = "n", horiz = TRUE, cex = 1
  )
  panel_label("B")
}

draw_go <- function(go_id, label, panel, show_y_label = FALSE, show_legend = FALSE) {
  one_term <- go_data[go_data$go_id == go_id, ]
  plot(
    NA, xlim = range(durations), ylim = go_limits,
    axes = FALSE, xlab = "Hours after cold onset",
    ylab = if (show_y_label) expression("Mean change (" * log[2] * ")") else "",
    main = ""
  )
  axis(1, at = durations, labels = durations)
  if (show_y_label) axis(2)
  box()
  abline(h = 0, col = "#BBBBBB", lwd = 0.8)
  for (treatment in names(colors)) {
    one <- one_term[one_term$treatment == treatment, ]
    one <- one[order(one$hours_after_onset), ]
    lines(
      one$hours_after_onset, one$mean_observed_change,
      type = "o", pch = 16, lwd = 1.25, col = colors[[treatment]]
    )
    lines(
      one$hours_after_onset, one$mean_fixed_prediction,
      lty = 2, lwd = 1.25, col = colors[[treatment]]
    )
  }
  mtext(label, side = 3, line = 0.15, adj = 0.5, font = 2, cex = 1)
  if (nzchar(panel)) panel_label(panel)
  if (show_legend) {
    legend(
      "top", inset = c(0, 0.01), legend = c("Observed", "Predicted"),
      col = "#222222", lty = c(1, 2), pch = c(16, NA),
      bty = "n", cex = 1, horiz = TRUE, x.intersp = 0.7
    )
  }
}

draw_partition <- function() {
  one <- held_out[match(durations, held_out$hours_after_onset), ]
  components <- t(as.matrix(one[, c(
    "shared_increment", "onset_increment", "diurnal_increment", "held_out_residual"
  )]))
  rownames(components) <- names(component_colors)
  if (any(components < 0) || max(abs(colSums(components) - 1)) > 1e-8) {
    stop("Held-out components are not a non-negative unit-sum partition.")
  }
  barplot(
    components, col = component_colors, border = "white", space = 0.38,
    ylim = c(0, 1), names.arg = paste0(durations, " h"),
    ylab = "Variance Explained",
    xlab = "Hours after cold onset"
  )
  legend(
    "top", inset = c(0, -0.19), legend = names(component_colors),
    fill = component_colors, border = NA, bty = "n", cex = 0.95,
    horiz = TRUE, xpd = NA, x.intersp = 0.4, text.width = NA
  )
  panel_label("D")
}

categories <- read.csv(
  gzfile(file.path(analysis_dir, "data", "generatable", "diurnal_adjusted_response_categories.csv.gz")),
  stringsAsFactors = FALSE
)
late <- categories[categories$hours_after_onset == 27, ]
phase <- categories[categories$collection_phase_specific_category != "not_evaluated", ]
count_data <- data.frame(
  group = c(rep("27 h onset", 3), rep("3 + 15 h light vs. dark", 2)),
  category = c("shared_induction", "dawn_specific_induction", "dusk_specific_induction",
               "light_cold_specific_induction", "dark_cold_specific_induction"),
  genes = c(sum(late$called_category == "shared_induction"),
            sum(late$called_category == "dawn_specific_induction"),
            sum(late$called_category == "dusk_specific_induction"),
            sum(phase$collection_phase_specific_category == "light_cold_specific_induction"),
            sum(phase$collection_phase_specific_category == "dark_cold_specific_induction")),
  stringsAsFactors = FALSE
)
count_colors <- c(Shared = "#E69F00", "Dawn-specific" = "#CC79A7", "Dusk-specific" = "#56B4E9",
                  "Light-cold" = "#F0E442", "Dark-cold" = "#5E3C99")
count_text <- c("#222222", "#222222", "#222222", "#222222", "white")

draw_counts <- function() {
  m <- matrix(0, nrow = 5, ncol = 2)
  m[1:3, 1] <- count_data$genes[1:3]; m[4:5, 2] <- count_data$genes[4:5]
  top <- max(colSums(m))
  ymax <- ceiling(top * 1.8 / 200) * 200
  mids <- barplot(m, col = count_colors, border = "white", space = 0.3, ylim = c(0, ymax),
                  xlim = c(0, 3.0), names.arg = rep("", 2), ylab = "Induced genes")
  mtext(c("27 h\nOnset", "3 + 15 h\nLight vs. dark"), side = 1, at = mids, line = 1.6, cex = 0.9)
  for (j in 1:2) {
    y <- cumsum(m[, j]); y0 <- c(0, head(y, -1))
    for (i in which(m[, j] > 0)) text(mids[j], (y0[i] + y[i]) / 2, m[i, j], cex = 0.9, font = 2, col = count_text[i])
  }
  legend("topleft", legend = names(count_colors), fill = count_colors, border = NA,
         bty = "n", cex = 0.9, x.intersp = 0.4, y.intersp = 0.9, inset = c(-0.03, -0.02))
  panel_label("E")
}

make_figure <- function(device_function, output_file, ...) {
  device_function(output_file, ...)
  layout(
    matrix(c(
      1, 2, 2, 2,
      0, 0, 0, 0,
      3, 3, 4, 5,
      0, 0, 0, 0,
      6, 6, 6, 7
    ), nrow = 5, byrow = TRUE),
    widths = c(0.333, 0.049, 0.309, 0.309),
    heights = c(1.02, 0.09, 0.92, 0.09, 1.02)
  )
  par(
    family = "sans", col = "#222222", fg = "#222222", cex = 1,
    mgp = c(2.05, 0.55, 0), tcl = -0.25
  )

  par(mar = c(4.1, 3.8, 0.55, 0.35))
  draw_scatter("DUSK", "A")
  par(mar = c(3.65, 3.8, 0.55, 0.35))
  draw_agreement()

  par(mar = c(3.15, 3.8, 1.35, 0.25))
  draw_go(names(selected_go)[1], selected_go[1], "C", TRUE, FALSE)
  par(mar = c(3.15, 0.65, 1.35, 0.25))
  draw_go(names(selected_go)[2], selected_go[2], "")
  par(mar = c(3.15, 0.65, 1.35, 0.25))
  draw_go(names(selected_go)[3], selected_go[3], "", FALSE, TRUE)

  par(mar = c(3.9, 4.15, 2.05, 0.45))
  draw_partition()
  par(mar = c(3.9, 4.15, 2.05, 0.3))
  draw_counts()
  dev.off()
}

png_file <- file.path(figure_dir, "SFigure4_longday_reference.png")
svg_file <- file.path(figure_dir, "SFigure4_longday_reference.svg")

make_figure(
  function(file, ...) png(filename = file, ...),
  png_file, width = 6.5, height = 6.9, units = "in", res = 300, pointsize = 9
)
make_figure(
  function(file, ...) svg(filename = file, ...),
  svg_file, width = 6.5, height = 6.9, pointsize = 9
)

write.csv(
  contrast_data, file.path(figure_dir, "SFigure4_fixed_prediction_summary.csv"),
  row.names = FALSE, quote = FALSE
)
write.csv(
  channel_data, file.path(figure_dir, "SFigure4_channel_agreement.csv"),
  row.names = FALSE, quote = FALSE
)
write.csv(
  go_data, file.path(figure_dir, "SFigure4_go_summary.csv"),
  row.names = FALSE, quote = TRUE
)
write.csv(
  held_out, file.path(figure_dir, "SFigure4_held_out_partition.csv"),
  row.names = FALSE, quote = FALSE
)
write.csv(
  fold_data, file.path(figure_dir, "SFigure4_held_out_fold_sse.csv"),
  row.names = FALSE, quote = FALSE
)

write.csv(count_data, file.path(figure_dir, "SFigure4_category_counts.csv"), row.names = FALSE, quote = FALSE)
cat("Fig. S4 analysis universe:", length(genes), "genes\n")
cat("DAWN-DUSK observed-response correlation at 11 h:", matched_11h_correlation, "\n")
cat("Fig. S4 outputs written to:\n", figure_dir, "\n")
