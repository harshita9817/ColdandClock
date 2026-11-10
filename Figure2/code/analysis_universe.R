# Shared definition of the cold-response analysis universe and the fixed
# external diurnal prediction.
#
# Figure 3 and the response-category table are two views of one analysis, so
# both take the gene set, the log2 expression matrix, the external waveform, the
# sample lookup, and the onset-hold prediction from here instead of rebuilding
# them independently. The universe size is derived once, in
# load_analysis_universe(), and is not asserted as a literal anywhere: it is a
# consequence of the expression threshold and the DESeq2 tested set, so a change
# in either should propagate rather than halt the pipeline. Every run reports
# the size it derived.

PSEUDOCOUNT <- 0.1
MINIMUM_EXTERNAL_MEAN_TPM <- 3
COLD_ONSET_ZT <- c(DAWN = 1, DUSK = 13)
DUSK_ONSET_OFFSET_H <- 12

# Evaluate a periodic cubic spline through the eight measured external phases.
# Accepts one gene's waveform as a vector, or a gene-by-phase matrix.
periodic_waveform <- function(values, phases, phase) {
  if (is.matrix(values)) {
    return(vapply(
      seq_len(nrow(values)),
      function(index) periodic_waveform(values[index, ], phases, phase),
      numeric(1)
    ))
  }
  splinefun(
    c(phases, 24), c(values, values[1]), method = "periodic"
  )(phase %% 24)
}

find_sample <- function(metadata, treatment, time_h, trial) {
  selected <- metadata$Treatment == treatment & metadata$Time == time_h &
    metadata$Rep. == trial
  if (sum(selected) != 1L) {
    stop("Sample lookup failed for ", treatment, ", T", time_h, ", trial ", trial)
  }
  metadata$`Sample name`[selected]
}

# Wall-clock hour of a sample collected `elapsed` hours after the named onset.
onset_time_h <- function(treatment, elapsed) {
  elapsed + ifelse(treatment == "DUSK", DUSK_ONSET_OFFSET_H, 0)
}

# `normalization` selects what is held constant across libraries. "tpm" keeps the
# total transcript pool fixed, which is what kallisto TPM already does.
# "size_factor" instead holds the typical gene fixed, applying median-of-ratios
# size factors to the TPM matrix. Gene length is handled by TPM either way, so
# this isolates the per-library scaling convention. Neither quantity is
# measurable without external spike-ins; the choice is a stated convention.
load_analysis_universe <- function(root, quiet = FALSE,
                                   normalization = c("tpm", "size_factor")) {
  normalization <- match.arg(normalization)
  reference <- read.delim(
    gzfile(file.path(root, "data/generatable/external_diurnal_reference.tsv.gz")),
    check.names = FALSE, stringsAsFactors = FALSE
  )
  metadata <- read.csv(
    file.path(root, "data/provided/sorghum_cold_timecourse_sample_metadata.csv"),
    check.names = FALSE, stringsAsFactors = FALSE
  )
  current <- read.csv(
    gzfile(file.path(root, "data/externalsourcerequired/kallisto_tpms.csv.gz")),
    check.names = FALSE, stringsAsFactors = FALSE
  )
  deg_results <- read.csv(
    gzfile(file.path(root, "data/generatable/deseq2_timepoint_contrasts.csv.gz")),
    check.names = FALSE, stringsAsFactors = FALSE
  )
  if (!identical(names(current)[-1], metadata$`Sample name`)) {
    stop("Current TPM columns and sample metadata are not aligned.")
  }
  if (!all(c("gene_id", "treatment", "time_h") %in% names(deg_results))) {
    stop("DESeq2 result table is missing required tested-universe columns.")
  }

  current_id <- sub("\\.[0-9]+$", "", current$GeneID)
  deg_id <- sub("\\.[0-9]+$", "", deg_results$gene_id)
  genes <- Reduce(intersect, list(
    reference$gene_id[reference$mean_external_tpm >= MINIMUM_EXTERNAL_MEAN_TPM],
    current_id,
    unique(deg_id)
  ))
  if (length(genes) == 0L) {
    stop("The expression/DESeq2 analysis universe is empty.")
  }

  reference <- reference[match(genes, reference$gene_id), ]
  current_log <- log2(
    as.matrix(current[match(genes, current_id), -1, drop = FALSE]) + PSEUDOCOUNT
  )
  storage.mode(current_log) <- "numeric"
  rownames(current_log) <- genes
  colnames(current_log) <- names(current)[-1]

  waveform_columns <- grep("^waveform_log2_tpm_zt", names(reference), value = TRUE)
  phases <- as.numeric(sub("waveform_log2_tpm_zt", "", waveform_columns))
  waveform <- as.matrix(reference[, waveform_columns, drop = FALSE])
  storage.mode(waveform) <- "numeric"
  rownames(waveform) <- genes

  if (normalization == "size_factor") {
    log_geometric_mean <- rowMeans(current_log)
    size_factor <- apply(current_log, 2, function(library) {
      median(library - log_geometric_mean)
    })
    current_log <- sweep(current_log, 2, size_factor, "-")
  }

  if (!quiet) {
    cat(sprintf("Normalization: %s\n", normalization))
    cat(sprintf(
      "Analysis universe: %d genes (mean external TPM >= %g, DESeq2 tested)\n",
      length(genes), MINIMUM_EXTERNAL_MEAN_TPM
    ))
  }
  list(
    genes = genes, current_log = current_log, metadata = metadata,
    reference = reference, waveform = waveform, phases = phases
  )
}

# Cold-minus-control log2 response for one treatment and wall-clock hour,
# as a gene-by-growth-trial matrix.
response_matrix <- function(universe, treatment, time_h) {
  matrix <- vapply(1:3, function(trial) {
    universe$current_log[, find_sample(universe$metadata, treatment, time_h, trial)] -
      universe$current_log[, find_sample(universe$metadata, "CONT", time_h, trial)]
  }, numeric(length(universe$genes)))
  rownames(matrix) <- universe$genes
  matrix
}

# Expression change predicted if a cold-treated leaf holds its onset state while
# its time-matched control continues through the external diurnal waveform.
fixed_onset_hold <- function(universe, treatment, time_h) {
  prediction <-
    periodic_waveform(universe$waveform, universe$phases, COLD_ONSET_ZT[[treatment]]) -
    periodic_waveform(universe$waveform, universe$phases, (time_h + 1) %% 24)
  names(prediction) <- universe$genes
  prediction
}
