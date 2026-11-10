#!/usr/bin/env Rscript

# Figure 2: the frozen-clock model of the sorghum cold response.
#   A  fixed diurnal-arrest prediction vs observed response, 11 h after dusk onset
#   B  variance explained by the prediction across the time course (Dawn, Dusk)
#   C  observed and predicted mean change of three GO terms
#   D  held-out variance partition: shared, specific, diurnal, residual
#   E  induced genes in the shared / Dawn-specific / Dusk-specific categories at 27 h
#
# Writes, inside the Figure2 folder:
#   Figure2.png, Figure2.svg            assembled figure
#   panels/Figure2<A-E>.png             each panel on its own
#   results/panel_<A-E>/*.csv           the numbers plotted in each panel
#
# Run from any directory:  Rscript code/plot_figure2.R
#
# Manuscript figure combining the fixed external diurnal-arrest test, biological
# examples, and held-out comparison of three response models. The analysis
# universe is the intersection of genes averaging at least 3 TPM in the
# independent series and genes tested by DESeq2; JTK_CYCLE significance is not
# an inclusion criterion.

script_argument <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_argument) != 1L) {
  stop("Run this figure script with Rscript so repository-relative paths can be resolved.")
}
# Rscript writes spaces in --file as "~+~"; the gsub() undoes that so paths with spaces work.
script_path <- normalizePath(
  gsub("~+~", " ", sub("^--file=", "", script_argument), fixed = TRUE), mustWork = TRUE
)
analysis_dir <- dirname(dirname(script_path))   # the Figure2 folder
figure_dir <- analysis_dir
results_dir <- file.path(analysis_dir, "results")
panels_dir <- file.path(analysis_dir, "panels")
source(file.path(analysis_dir, "code", "component_aligned_diurnal_model.R"))
source(file.path(analysis_dir, "code", "analysis_universe.R"))
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
onset_labels <- c(DAWN = "Dawn", DUSK = "Dusk")
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

# Residual response-category sizes (shared and onset-specific induction only;
# collection-phase categories are omitted) from the downstream classification,
# which uses the same gene universe and external waveforms.
count_hours <- c(3, 27)
count_categories <- c(
  "Shared" = "shared_induction",
  "Dawn-specific" = "dawn_specific_induction",
  "Dusk-specific" = "dusk_specific_induction"
)
# Colors not used elsewhere in the figure, so these categories are not read as
# the onsets (B, C) or the model components (D).
count_colors <- c(
  "Shared" = "#E69F00", "Dawn-specific" = "#CC79A7", "Dusk-specific" = "#56B4E9"
)
categories <- read.csv(
  gzfile(file.path(
    analysis_dir, "data", "generatable", "diurnal_adjusted_response_categories.csv.gz"
  )),
  stringsAsFactors = FALSE
)
categories <- categories[
  categories$hours_after_onset %in% count_hours & categories$gene_id %in% genes,
]
category_counts <- sapply(count_hours, function(elapsed) {
  called <- categories$called_category[categories$hours_after_onset == elapsed]
  vapply(count_categories, function(category) sum(called == category), numeric(1))
})
colnames(category_counts) <- paste0(count_hours, " h")
# Only the 27-h bar is plotted; both hours are still written to the counts CSV.
plotted_counts <- category_counts[, "27 h", drop = FALSE]

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
    "topright", legend = onset_labels[names(colors)], col = colors, pch = 16,
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
  # Per-entry text widths keep the horizontal legend within the narrower panel.
  legend_cex <- 0.85
  legend(
    "top", inset = c(0, -0.17), legend = names(component_colors),
    fill = component_colors, border = NA, bty = "n", cex = legend_cex,
    horiz = TRUE, xpd = NA, x.intersp = 0.5,
    text.width = strwidth(names(component_colors), cex = legend_cex)
  )
  panel_label("D")
}

draw_category_counts <- function() {
  # End the y axis at its last tick so the axis line spans the full plot
  # height, matching panel D.
  y_ticks <- pretty(c(0, max(colSums(plotted_counts))))
  y_ticks <- y_ticks[seq_len(which(y_ticks >= max(colSums(plotted_counts)))[1])]
  midpoints <- barplot(
    plotted_counts, col = count_colors, border = "white", space = 0.12, width = 1,
    xlim = c(0, 2.35), ylim = range(y_ticks), names.arg = colnames(plotted_counts),
    ylab = "Induced genes", xlab = "", las = 1, yaxt = "n"
  )
  axis(2, at = y_ticks, las = 1)
  # Centre the axis title on the whole panel; the plot region alone is too
  # narrow to hold it.
  mtext(
    "Hours after cold onset", side = 1, line = par("mgp")[1],
    at = max(midpoints)
  )
  tops <- apply(plotted_counts, 2, cumsum)
  bottoms <- tops - plotted_counts
  for (column in seq_len(ncol(plotted_counts))) {
    for (row in seq_len(nrow(plotted_counts))) {
      n <- plotted_counts[row, column]
      if (n >= 30) {
        text(
          midpoints[column], (tops[row, column] + bottoms[row, column]) / 2,
          labels = n, col = "#222222", font = 2, cex = 0.75
        )
      }
    }
    small <- plotted_counts[, column] < 30
    if (any(small)) {
      text(
        midpoints[column], tops[nrow(plotted_counts), column],
        labels = paste(plotted_counts[small, column], tolower(names(which(small))),
                       collapse = "\n"),
        pos = 3, cex = 0.75, col = "#222222"
      )
    }
  }
  # The legend sits to the right of the single 27-h bar.
  legend(
    x = max(midpoints) + 0.55, y = par("usr")[4],
    legend = names(count_colors), fill = count_colors, border = NA,
    bty = "n", cex = 0.85, x.intersp = 0.5, y.intersp = 0.9, xpd = NA
  )
  panel_label("E")
}

make_figure <- function(device_function, output_file, ...) {
  device_function(output_file, ...)
  layout(
    matrix(c(
      1, 2, 2, 2, 2, 2,
      0, 0, 0, 0, 0, 0,
      3, 3, 4, 4, 4, 5,
      0, 0, 0, 0, 0, 0,
      6, 6, 6, 6, 7, 7
    ), nrow = 5, byrow = TRUE),
    # Panel E holds a single 27-h bar plus its key; the fifth column (part of the
    # second GO panel) lets E be a little wider than the third GO panel.
    widths = c(0.333, 0.049, 0.218, 0.021, 0.07, 0.309),
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
  par(mar = c(3.9, 3.6, 2.05, 0.05))
  draw_category_counts()
  dev.off()
}

png_file <- file.path(figure_dir, "Figure2.png")
svg_file <- file.path(figure_dir, "Figure2.svg")

make_figure(
  function(file, ...) png(filename = file, ...),
  png_file, width = 6.5, height = 6.9, units = "in", res = 300, pointsize = 9
)
make_figure(
  # Prefer svglite: the Cairo svg() device fails where XQuartz is missing.
  function(file, ...) {
    if (requireNamespace("svglite", quietly = TRUE)) {
      svglite::svglite(file, ...)
    } else {
      svg(filename = file, ...)
    }
  },
  svg_file, width = 6.5, height = 6.9, pointsize = 9
)

# ---- Each panel on its own -------------------------------------------------
# Same drawing functions and margins as the assembled figure.
panel_png <- function(name, width, height, draw) {
  png(file.path(panels_dir, paste0("Figure2", name, ".png")),
      width = width, height = height, units = "in", res = 300, pointsize = 9)
  par(family = "sans", col = "#222222", fg = "#222222", cex = 1,
      mgp = c(2.05, 0.55, 0), tcl = -0.25)
  draw()
  dev.off()
}
dir.create(panels_dir, showWarnings = FALSE, recursive = TRUE)
panel_png("A", 2.4, 2.35, function() {
  par(mar = c(4.1, 3.8, 0.55, 0.35)); draw_scatter("DUSK", "A")
})
panel_png("B", 4.3, 2.35, function() {
  par(mar = c(3.65, 3.8, 0.55, 0.35)); draw_agreement()
})
panel_png("C", 6.5, 2.2, function() {
  layout(matrix(1:3, nrow = 1), widths = c(0.382, 0.309, 0.309))
  par(cex = 1)  # layout() shrinks text for three columns; keep it at full size
  par(mar = c(3.15, 3.8, 1.35, 0.25))
  draw_go(names(selected_go)[1], selected_go[1], "C", TRUE, FALSE)
  par(mar = c(3.15, 0.65, 1.35, 0.25))
  draw_go(names(selected_go)[2], selected_go[2], "")
  draw_go(names(selected_go)[3], selected_go[3], "", FALSE, TRUE)
})
panel_png("D", 4.1, 2.45, function() {
  par(mar = c(3.9, 4.15, 2.05, 0.45)); draw_partition()
})
panel_png("E", 2.4, 2.45, function() {
  par(mar = c(3.9, 3.6, 2.05, 0.05)); draw_category_counts()
})

# ---- The numbers plotted in each panel ---------------------------------------
save_table <- function(data, panel, name) {
  folder <- file.path(results_dir, paste0("panel_", panel))
  dir.create(folder, showWarnings = FALSE, recursive = TRUE)
  write.csv(data, file.path(folder, name), row.names = FALSE)
}
# A: every gene, 11 h after dusk onset (wall-clock T23).
panel_a <- gene_data[gene_data$treatment == "DUSK" & gene_data$hours_after_onset == 11,
                     c("gene_id", "fixed_onset_hold_prediction",
                       "observed_cold_minus_control_log2_tpm")]
names(panel_a) <- c("gene_id", "predicted_response_log2", "observed_response_log2")
panel_a_fit <- lm(observed_response_log2 ~ predicted_response_log2, data = panel_a)
save_table(panel_a, "A", "panelA_gene_points.csv")
save_table(data.frame(
  treatment = "DUSK", hours_after_onset = 11, genes = nrow(panel_a),
  pearson_r = cor(panel_a$predicted_response_log2, panel_a$observed_response_log2),
  line_intercept = unname(coef(panel_a_fit)[1]), line_slope = unname(coef(panel_a_fit)[2])
), "A", "panelA_correlation_and_line.csv")
# B: one point per treatment and time.
save_table(contrast_data, "B", "panelB_variance_explained.csv")
save_table(channel_data, "B", "panelB_channel_agreement_not_plotted.csv")
# C: one row per GO term, treatment and time.
save_table(go_data, "C", "panelC_go_term_mean_change.csv")
# D: one bar per time; fold-level sums of squares behind it.
save_table(held_out, "D", "panelD_variance_partition.csv")
save_table(fold_data, "D", "panelD_held_out_fold_sums_of_squares.csv")
# E: only 27 h is plotted; 3 h is kept for reference.
save_table(data.frame(
  hours_after_onset = rep(count_hours, each = length(count_categories)),
  category = rep(unname(count_categories), length(count_hours)),
  label = rep(names(count_categories), length(count_hours)),
  genes = as.vector(category_counts),
  plotted = rep(count_hours, each = length(count_categories)) == 27
), "E", "panelE_induced_gene_counts.csv")
writeLines(c(
  sprintf("Analysis universe: %d genes", length(genes)),
  sprintf("Dawn-Dusk observed-response correlation at 11 h: %.4f", matched_11h_correlation),
  sprintf("R version: %s", R.version.string)
), file.path(results_dir, "run_summary.txt"))

cat("Figure 2 analysis universe:", length(genes), "genes\n")
cat("DAWN-DUSK observed-response correlation at 11 h:", matched_11h_correlation, "\n")
cat("Figure 2 outputs written to:\n", figure_dir, "\n")
