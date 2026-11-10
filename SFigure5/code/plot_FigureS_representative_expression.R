#!/usr/bin/env Rscript

# Supplementary figure (Fig. S5): expression of representative genes.
#
# Stand-alone version of the supplementary panel in plot_Figure4_2.R, so the
# supplement can be redrawn without re-running (and re-writing) Figure 4.
# The gene choice and the plotted values follow plot_Figure4_2.R exactly; the
# differences are:
#   - only the 27-h examples (shared, Dawn-specific, Dusk-specific) and SbPPCK1
#     are shown; the two 3-h panels are dropped;
#   - on-figure labels use "Dawn"/"Dusk" instead of "DAWN"/"DUSK".
#
# Run with Rscript from any directory. Requires ggplot2 and systemfonts.
# Outputs, inside this figure's folder:
#   figure/FigureS_representative_expression.png / .pdf / .svg
#   data/plotted/FigureS_representative_expression_summary.csv
#   data/plotted/FigureS_representative_expression_genes.csv

suppressPackageStartupMessages(library(ggplot2))
arg <- grep("^--file=", commandArgs(FALSE), value=TRUE)
if (length(arg) != 1L) stop("Run this script with Rscript.")
# Rscript writes spaces in --file as "~+~"; undo that so paths with spaces work.
script <- normalizePath(gsub("~+~", " ", sub("^--file=", "", arg), fixed=TRUE), mustWork=TRUE)
# Layout: code/ (this script), data/ (inputs; plotted numbers go to data/plotted/),
# figure/ (the images). Paths are relative to this figure's folder.
source_dir <- dirname(dirname(script))   # this figure's folder
table_dir <- file.path(source_dir, "data", "plotted")
out_dir <- file.path(source_dir, "figure")
dir.create(table_dir, recursive=TRUE, showWarnings=FALSE)
dir.create(out_dir, showWarnings=FALSE)

paths <- c(
  responses=file.path(source_dir, "data/generatable/diurnal_adjusted_response_categories.csv.gz"),
  tpm=file.path(source_dir, "data/externalsourcerequired/kallisto_tpms.csv.gz"),
  metadata=file.path(source_dir, "data/provided/sorghum_cold_timecourse_sample_metadata.csv")
)
stopifnot(all(file.exists(paths)))
responses <- read.csv(gzfile(paths["responses"]), check.names=FALSE, stringsAsFactors=FALSE)
metadata <- read.csv(paths["metadata"], check.names=FALSE, stringsAsFactors=FALSE)
stopifnot(nrow(responses) == 30260L, length(unique(responses$gene_id)) == 15130L)

category_labels <- c(shared_induction="Shared", dawn_specific_induction="Dawn-specific",
                     dusk_specific_induction="Dusk-specific")

# ---- Genes: strongest 27-h example per category, plus SbPPCK1 ---------------
# Same rule as plot_Figure4_2.R: lowest category BH q, ties broken by the
# largest absolute DAWN-DUSK difference, then gene id.
spec <- data.frame(hours=27, category=names(category_labels), stringsAsFactors=FALSE)
reps <- do.call(rbind, lapply(seq_len(nrow(spec)), function(i) {
  d <- responses[responses$hours_after_onset == spec$hours[i] &
                   responses$called_category == spec$category[i], ]
  d <- d[order(d$category_q_value, -abs(d$dawn_minus_dusk_log2), d$gene_id), ]
  d[1, c("gene_id", "hours_after_onset", "called_category", "category_q_value")]
}))
reps$selection_reason <- "lowest category BH q; tie by largest absolute DAWN-DUSK difference"
ppck_id <- "Sobic.004G338000"
ppck_rep <- responses[responses$gene_id == ppck_id & responses$hours_after_onset == 27,
                      c("gene_id", "hours_after_onset", "called_category", "category_q_value")]
stopifnot(nrow(ppck_rep) == 1L)
ppck_rep$selection_reason <- "pre-specified SbPPCK1 biological example"
reps <- rbind(reps, ppck_rep)
reps$facet <- paste0(reps$hours_after_onset, " h\n",
                     category_labels[reps$called_category], "\n", reps$gene_id)
reps$facet[reps$gene_id == ppck_id] <- paste0("SbPPCK1\n", ppck_id)
rownames(reps) <- NULL
write.csv(reps, file.path(table_dir, "FigureS_representative_expression_genes.csv"), row.names=FALSE)

# ---- Raw TPM for the displayed genes -----------------------------------------
tpm <- read.csv(gzfile(paths["tpm"]), check.names=FALSE, stringsAsFactors=FALSE)
tpm$gene_id <- sub("\\.[0-9]+$", "", tpm$GeneID)
tpm <- tpm[match(reps$gene_id, tpm$gene_id), ]
stopifnot(!anyNA(tpm$gene_id))
find_sample <- function(treatment, time, trial) {
  z <- metadata$Treatment == treatment & metadata$Time == time & metadata$Rep. == trial
  stopifnot(sum(z) == 1L)
  metadata$`Sample name`[z]
}

# One experimental-time axis starting at dawn cold onset. Dusk cold began 12 h
# later, so its 3-h sample sits at experimental hour 15.
elapsed <- c(3, 7, 11, 15, 19, 23, 27)
rows <- list()
for (i in seq_len(nrow(reps))) {
  for (onset in c("DAWN", "DUSK")) for (e in elapsed) for (trial in 1:3) {
    wall <- e + if (onset == "DUSK") 12 else 0
    cold <- find_sample(onset, wall, trial)
    control <- find_sample("CONT", wall, trial)
    rows[[length(rows) + 1L]] <- data.frame(
      gene_id=reps$gene_id[i], trial=trial, experiment_h=wall,
      series=if (onset == "DAWN") "Cold: Dawn" else "Cold: Dusk",
      sample_name=cold, tpm=as.numeric(tpm[i, cold]), stringsAsFactors=FALSE)
    rows[[length(rows) + 1L]] <- data.frame(
      gene_id=reps$gene_id[i], trial=trial, experiment_h=wall, series="Control",
      sample_name=control, tpm=as.numeric(tpm[i, control]), stringsAsFactors=FALSE)
  }
  # Sequenced T-1 pretreatment controls, as descriptive reference points.
  for (trial in 1:3) {
    control <- find_sample("CONT", -1, trial)
    rows[[length(rows) + 1L]] <- data.frame(
      gene_id=reps$gene_id[i], trial=trial, experiment_h=-1, series="Control",
      sample_name=control, tpm=as.numeric(tpm[i, control]), stringsAsFactors=FALSE)
  }
}
raw <- do.call(rbind, rows)
# Controls at wall-clock times shared by the two onset series appear twice;
# keep one observation per gene and sample.
raw <- raw[!duplicated(raw[, c("gene_id", "sample_name", "series")]), ]

# ---- Mean +/- SE across the three growth trials ------------------------------
expr_summary <- aggregate(tpm ~ gene_id + series + experiment_h, raw,
                          function(x) c(mean=mean(x), se=sd(x) / sqrt(length(x)), n=length(x)))
expr_summary <- data.frame(
  gene_id=expr_summary$gene_id,
  facet=reps$facet[match(expr_summary$gene_id, reps$gene_id)],
  series=expr_summary$series, experiment_h=expr_summary$experiment_h,
  mean_tpm=expr_summary$tpm[, "mean"], se=expr_summary$tpm[, "se"],
  n=expr_summary$tpm[, "n"], stringsAsFactors=FALSE)
stopifnot(all(expr_summary$n == 3))
series_levels <- c("Control", "Cold: Dawn", "Cold: Dusk")
expr_summary$series <- factor(expr_summary$series, levels=series_levels)
expr_summary$facet <- factor(expr_summary$facet, levels=reps$facet)
expr_summary <- expr_summary[order(expr_summary$facet, expr_summary$series, expr_summary$experiment_h), ]
write.csv(expr_summary, file.path(table_dir, "FigureS_representative_expression_summary.csv"),
          row.names=FALSE)

# Each cold line starts from a control point so the treatment visibly branches
# from the control trajectory: dawn from the T-1 control, dusk from the T11
# control (one hour before dusk cold began). Symbols and error bars are drawn
# only for actual cold samples.
anchor <- function(hour, series) {
  a <- expr_summary[expr_summary$series == "Control" & expr_summary$experiment_h == hour, ]
  a$series <- factor(series, levels=series_levels)
  a
}
line_summary <- rbind(expr_summary, anchor(-1, "Cold: Dawn"), anchor(11, "Cold: Dusk"))

# ---- Plot --------------------------------------------------------------------
# Arial where installed, otherwise Liberation Sans (same letter widths).
font_family <- if (any(grepl("^Arial$", systemfonts::system_fonts()$family))) "Arial" else "Liberation Sans"
supp_cols <- c("Control"="#6E6E6E", "Cold: Dawn"="#C51B7D", "Cold: Dusk"="#008B72")
supp_types <- c("Control"="solid", "Cold: Dawn"="solid", "Cold: Dusk"="22")
supp_shapes <- c("Control"=1, "Cold: Dawn"=16, "Cold: Dusk"=15)

p <- ggplot(expr_summary, aes(experiment_h, mean_tpm, color=series,
                              linetype=series, shape=series, group=series)) +
  geom_errorbar(aes(ymin=pmax(0, mean_tpm - se), ymax=mean_tpm + se), width=.4, linewidth=.3) +
  geom_line(data=line_summary, linewidth=.65) +
  geom_point(size=1.7) +
  facet_wrap(~facet, ncol=1, scales="free_y") +
  scale_color_manual(values=supp_cols) +
  scale_linetype_manual(values=supp_types) +
  scale_shape_manual(values=supp_shapes) +
  scale_x_continuous(breaks=c(-1, 3, 7, 11, 15, 19, 23, 27, 31, 35, 39), limits=c(-1, 39)) +
  labs(x="Experimental time (hours after dawn cold onset)", y="Expression (TPM)",
       color=NULL, linetype=NULL, shape=NULL) +
  theme_classic(base_size=9, base_family=font_family) +
  theme(axis.title=element_text(size=9), axis.text=element_text(size=9),
        legend.text=element_text(size=9), legend.position="top",
        strip.background=element_rect(fill="#F2F2F2", color="#AAAAAA"),
        strip.text=element_text(face="bold"),
        panel.spacing=unit(2.2, "mm"), plot.margin=margin(5, 5, 5, 5)) +
  guides(color=guide_legend(nrow=1, byrow=TRUE),
         linetype=guide_legend(nrow=1, byrow=TRUE),
         shape=guide_legend(nrow=1, byrow=TRUE))

# Four panels instead of six: keep the per-panel height of the earlier 9-in version.
width <- 6.5
height <- 6.4
base <- file.path(out_dir, "FigureS_representative_expression")
ggsave(paste0(base, ".png"), p, width=width, height=height, dpi=300, bg="white")
# cairo_pdf embeds system fonts; the base pdf() device does not know Liberation Sans.
ggsave(paste0(base, ".pdf"), p, width=width, height=height, bg="white", device=cairo_pdf)
ggsave(paste0(base, ".svg"), p, width=width, height=height, bg="white")
cat("Fig. S5 written to", out_dir, "\n")
print(reps[, c("gene_id", "called_category", "category_q_value")], row.names=FALSE)
