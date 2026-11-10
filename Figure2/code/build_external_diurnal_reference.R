#!/usr/bin/env Rscript

# Build the canonical gene-level reference from the independent 24-sample
# BTx623 diurnal series. Every detected gene receives JTK_CYCLE results and an
# empirical eight-phase waveform; rhythmic significance is descriptive and is
# never an inclusion requirement for downstream cold-response modeling.

suppressPackageStartupMessages(library(MetaCycle))
if (as.character(packageVersion("MetaCycle")) != "1.2.1") {
  stop("This analysis is pinned to MetaCycle 1.2.1.")
}

arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
# Rscript writes spaces in --file as "~+~"; the gsub() undoes that so paths with spaces work.
root <- dirname(dirname(normalizePath(gsub("~+~", " ", sub("^--file=", "", arg), fixed = TRUE), mustWork = TRUE)))
input <- file.path(root, "data/externalsourcerequired/sorghum_2020_diurnal_tpms.csv.gz")
output <- file.path(root, "data/generatable/external_diurnal_reference.tsv.gz")
dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)

x <- read.csv(gzfile(input), check.names = FALSE, stringsAsFactors = FALSE)
if (ncol(x) != 25L) stop("Expected one gene-ID column and 24 samples.")
gene_id <- as.character(x[[1]])
normalized_gene_id <- sub("\\.$", "", gene_id)
if (anyDuplicated(gene_id) || anyDuplicated(normalized_gene_id)) {
  stop("External source or terminal-period-normalized IDs are not unique.")
}
sample_number <- as.integer(sub("^Sorghum", "", names(x)[-1]))
o <- order(sample_number); sample_number <- sample_number[o]
expression <- as.matrix(x[, 1L + o, drop = FALSE]); storage.mode(expression) <- "numeric"
rownames(expression) <- gene_id
collection_hour <- sample_number * 3
zt <- collection_hour %% 24
cycle <- ceiling(collection_hour / 24)
detected <- rowSums(expression > 0) > 0

jtk_input <- data.frame(gene_id=gene_id[detected], expression[detected,,drop=FALSE], check.names=FALSE)
jtk_run <- meta2d(
  infile="in_memory_diurnal_tpms.csv", filestyle="csv",
  timepoints=collection_hour, minper=20, maxper=28, cycMethod="JTK",
  analysisStrategy="selfUSE", outputFile=FALSE, outIntegration="both",
  adjustPhase="predictedPer", releaseNote=FALSE, parallelize=FALSE, inDF=jtk_input
)
jtk <- as.data.frame(jtk_run$JTK, stringsAsFactors=FALSE)
if (!identical(as.character(jtk$CycID), jtk_input$gene_id)) stop("JTK output alignment failed.")

# Remove complete-day offsets while retaining each gene's grand mean, then
# average the three observations at each measured phase. This is done for all
# detected genes, not only genes called rhythmic.
logx <- log2(expression[detected,,drop=FALSE] + 0.1)
grand <- rowMeans(logx); adjusted <- logx
for (d in sort(unique(cycle))) {
  use <- cycle == d
  adjusted[,use] <- sweep(logx[,use,drop=FALSE], 1, rowMeans(logx[,use,drop=FALSE]), "-")
  adjusted[,use] <- sweep(adjusted[,use,drop=FALSE], 1, grand, "+")
}
phases <- sort(unique(zt))
wave <- vapply(phases, function(p) rowMeans(adjusted[,zt==p,drop=FALSE]), numeric(nrow(adjusted)))
colnames(wave) <- paste0("waveform_log2_tpm_zt", phases)

result <- data.frame(
  gene_id=jtk$CycID,
  normalized_gene_id=sub("\\.$", "", jtk$CycID),
  mean_external_tpm=rowMeans(expression[detected,,drop=FALSE]),
  detected_samples=rowSums(expression[detected,,drop=FALSE] > 0),
  jtk_bh_q=jtk$BH.Q, jtk_adjusted_p=jtk$ADJ.P,
  jtk_period_h=jtk$PER, jtk_lag_h=jtk$LAG, jtk_amplitude_tpm=jtk$AMP,
  rhythmic_adjusted_p_le_0_01=jtk$ADJ.P <= 0.01 & jtk$PER > 0,
  wave, metacycle_version=as.character(packageVersion("MetaCycle")),
  stringsAsFactors=FALSE, check.names=FALSE
)
con <- gzfile(output, "wt", compression=9)
write.table(result, con, sep="\t", row.names=FALSE, quote=FALSE); close(con)
cat("Wrote", nrow(result), "genes to", output, "\n")
