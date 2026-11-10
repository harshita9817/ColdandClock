# Differential expression across the sorghum dawn/dusk cold time course.
#
# One DESeq2 model is fit to the 72 libraries used in the 14 reported
# post-treatment contrasts. Growth trial is included as a blocking factor, and
# each observed timepoint-treatment combination is represented by a group
# factor. No additional per-contrast feature filtering is applied.

suppressPackageStartupMessages({
  library(DESeq2)
})

alpha <- 0.05
lfc_cutoff <- 1
dawn_times <- c(3, 7, 11, 15, 19, 23, 27)
dusk_times <- c(15, 19, 23, 27, 31, 35, 39)
excluded_nonconvergent_genes <- "Sobic.004G220500.1"

script_argument <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_argument) != 1L) {
  stop("Run this analysis with Rscript so the repository-relative paths can be resolved.")
}
# Rscript writes spaces in --file as "~+~"; the gsub() undoes that so paths with spaces work.
script_path <- normalizePath(
  gsub("~+~", " ", sub("^--file=", "", script_argument), fixed = TRUE), mustWork = TRUE
)
analysis_dir <- dirname(dirname(script_path))

metadata_file <- file.path(
  analysis_dir, "data", "provided", "sorghum_cold_timecourse_sample_metadata.csv"
)
count_file <- file.path(
  analysis_dir, "data", "externalsourcerequired", "kallisto_counts_filtered.csv.gz"
)
output_dir <- file.path(analysis_dir, "data", "generatable")
output_file <- file.path(output_dir, "deseq2_timepoint_contrasts.csv.gz")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

metadata <- read.csv(metadata_file, check.names = FALSE, stringsAsFactors = FALSE)
required_metadata_columns <- c(
  "Sample name", "Species", "Time", "Treatment", "Rep.", "Collection time"
)
missing_metadata_columns <- setdiff(required_metadata_columns, names(metadata))
if (length(missing_metadata_columns) > 0L) {
  stop("Metadata is missing required columns: ", paste(missing_metadata_columns, collapse = ", "))
}
if (anyDuplicated(metadata[["Sample name"]])) {
  stop("Metadata contains duplicate sample names.")
}

counts_data <- read.csv(gzfile(count_file), check.names = FALSE, stringsAsFactors = FALSE)
if (names(counts_data)[1] != "GeneID") {
  stop("The first count-table column must be named GeneID.")
}
if (anyNA(counts_data[["GeneID"]]) || any(counts_data[["GeneID"]] == "")) {
  stop("The count table contains missing or empty gene IDs.")
}
if (anyDuplicated(counts_data[["GeneID"]])) {
  stop("The count table contains duplicate gene IDs.")
}

metadata_samples <- metadata[["Sample name"]]
count_samples <- names(counts_data)[-1]
if (!identical(count_samples, metadata_samples)) {
  stop("Count columns must exactly match the metadata sample IDs and order.")
}

counts <- as.matrix(counts_data[, -1, drop = FALSE])
storage.mode(counts) <- "numeric"
if (anyNA(counts) || any(!is.finite(counts))) {
  stop("The count matrix contains missing or non-finite values.")
}
if (any(counts < 0)) {
  stop("The count matrix contains negative values.")
}
if (any(counts != round(counts))) {
  stop("The count matrix contains non-integer values.")
}
if (any(counts > .Machine$integer.max)) {
  stop("The count matrix contains values larger than R's integer limit.")
}
storage.mode(counts) <- "integer"
rownames(counts) <- counts_data[["GeneID"]]

all_contrast_times <- sort(union(dawn_times, dusk_times))
use_sample <-
  (metadata[["Treatment"]] == "DAWN" & metadata[["Time"]] %in% dawn_times) |
  (metadata[["Treatment"]] == "DUSK" & metadata[["Time"]] %in% dusk_times) |
  (metadata[["Treatment"]] == "CONT" & metadata[["Time"]] %in% all_contrast_times)
analysis_metadata <- metadata[use_sample, , drop = FALSE]

if (nrow(analysis_metadata) != 72L) {
  stop("Expected 72 libraries in the joint analysis, found ", nrow(analysis_metadata), ".")
}

growth_trial_date <- sub(
  "^S-([0-9]{6})-.*$", "\\1", analysis_metadata[["Sample name"]]
)
if (any(nchar(growth_trial_date) != 6L)) {
  stop("Could not derive a six-digit growth-trial date from every sample ID.")
}
analysis_metadata$growth_trial_date <- factor(growth_trial_date)
analysis_metadata$group <- factor(sprintf(
  "T%02d_%s", analysis_metadata[["Time"]], analysis_metadata[["Treatment"]]
))
analysis_metadata$group <- relevel(analysis_metadata$group, ref = "T03_CONT")
rownames(analysis_metadata) <- analysis_metadata[["Sample name"]]

expected_groups <- c(
  sprintf("T%02d_CONT", all_contrast_times),
  sprintf("T%02d_DAWN", dawn_times),
  sprintf("T%02d_DUSK", dusk_times)
)
observed_group_counts <- table(analysis_metadata$group)
if (!setequal(names(observed_group_counts), expected_groups)) {
  stop("The observed timepoint-treatment groups do not match the expected design.")
}
if (any(observed_group_counts != 3L)) {
  stop("Every timepoint-treatment group must contain exactly three libraries.")
}
group_by_trial <- table(analysis_metadata$group, analysis_metadata$growth_trial_date)
if (any(group_by_trial != 1L)) {
  stop("Every timepoint-treatment group must contain one library from each growth trial.")
}

analysis_samples <- analysis_metadata[["Sample name"]]
analysis_counts <- counts[, analysis_samples, drop = FALSE]
all_zero <- rowSums(analysis_counts) == 0L
if (any(all_zero)) {
  stop(
    sum(all_zero),
    " genes have zero counts across all 72 analysis libraries; decide how to handle them before fitting."
  )
}
if (!all(excluded_nonconvergent_genes %in% rownames(analysis_counts))) {
  stop("An explicitly excluded nonconvergent gene is absent from the count matrix.")
}
analysis_counts <- analysis_counts[
  !rownames(analysis_counts) %in% excluded_nonconvergent_genes,
  ,
  drop = FALSE
]

dds <- DESeqDataSetFromMatrix(
  countData = analysis_counts,
  colData = analysis_metadata,
  design = ~ growth_trial_date + group
)
dds <- DESeq(dds)
nonconvergent_genes <- rownames(dds)[!mcols(dds)$betaConv]
if (length(nonconvergent_genes) > 0L) {
  stop(
    "DESeq2 coefficient fitting did not converge for: ",
    paste(nonconvergent_genes, collapse = ", ")
  )
}

contrasts <- rbind(
  data.frame(treatment = "DAWN", time_h = dawn_times),
  data.frame(treatment = "DUSK", time_h = dusk_times)
)

contrast_results <- vector("list", nrow(contrasts))
for (index in seq_len(nrow(contrasts))) {
  treatment <- contrasts$treatment[index]
  time_h <- contrasts$time_h[index]
  treatment_group <- sprintf("T%02d_%s", time_h, treatment)
  control_group <- sprintf("T%02d_CONT", time_h)

  result <- results(
    dds,
    contrast = c("group", treatment_group, control_group),
    alpha = alpha,
    independentFiltering = FALSE
  )
  result_data <- as.data.frame(result)
  result_data$gene_id <- rownames(result_data)
  result_data$treatment <- treatment
  result_data$time_h <- time_h
  result_data$contrast <- paste0(treatment, "_vs_CONT_T", time_h)
  result_data$fold_change <- 2^result_data$log2FoldChange
  result_data$de_status <- "Not DEG"

  is_deg <-
    !is.na(result_data$padj) &
    result_data$padj < alpha &
    abs(result_data$log2FoldChange) >= lfc_cutoff
  result_data$de_status[is_deg & result_data$log2FoldChange > 0] <- "Up"
  result_data$de_status[is_deg & result_data$log2FoldChange < 0] <- "Down"

  contrast_results[[index]] <- result_data[, c(
    "gene_id", "treatment", "time_h", "contrast", "baseMean",
    "log2FoldChange", "fold_change", "lfcSE", "stat", "pvalue", "padj",
    "de_status"
  )]
}

all_results <- do.call(rbind, contrast_results)
names(all_results) <- c(
  "gene_id", "treatment", "time_h", "contrast", "base_mean",
  "log2_fold_change", "fold_change", "lfc_se", "wald_statistic", "p_value",
  "adjusted_p_value", "de_status"
)
rownames(all_results) <- NULL

output_connection <- gzfile(output_file, open = "wt", compression = 9)
write.csv(all_results, output_connection, row.names = FALSE, quote = FALSE)
close(output_connection)

summary_counts <- as.data.frame(table(
  all_results$treatment, all_results$time_h, all_results$de_status
))
names(summary_counts) <- c("treatment", "time_h", "de_status", "genes")
summary_counts <- summary_counts[summary_counts$genes > 0, , drop = FALSE]

cat("DESeq2 timepoint contrasts written to:\n", output_file, "\n\n")
cat("Model: ~ growth_trial_date + group\n")
cat("Libraries:", ncol(analysis_counts), "\n")
cat("Genes:", nrow(analysis_counts), "\n")
cat("Excluded nonconvergent genes:", paste(excluded_nonconvergent_genes, collapse = ", "), "\n")
cat("Contrasts:", nrow(contrasts), "\n\n")
print(summary_counts, row.names = FALSE)
cat("\nR version:", R.version.string, "\n")
cat("DESeq2 version:", as.character(packageVersion("DESeq2")), "\n")
