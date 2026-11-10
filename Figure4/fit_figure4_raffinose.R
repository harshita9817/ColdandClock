#!/usr/bin/env Rscript

# Figure 4B: raffinose change after 27 hours of dawn- and dusk-onset chilling.
#
# Input: Figure4_raffinose_samples.csv, one row per GC-MS injection with
# log2 raffinose peak height normalized to sample weight (NA for the three
# failed injections). Times (T) are hours from dawn onset, so 27 h after dusk
# onset is T39.
#
# Model: robust (bisquare M-estimator) regression of log2 raffinose on
# treatment-by-time cell plus growth trial (date), GC-MS run, injection order,
# and the ribitol internal standard. Contrasts compare chilled samples with
# time-matched controls 27 h after each onset, and the dusk-minus-dawn
# difference of those changes.
#
# Output: Figure4_raffinose_model.csv

suppressMessages(library(MASS))
arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
figdir <- dirname(normalizePath(sub("^--file=", "", arg), mustWork = TRUE))

d <- read.csv(file.path(figdir, "Figure4_raffinose_samples.csv"), colClasses = c(date = "character"))
d <- d[!is.na(d$log2_raffinose), ]
X <- model.matrix(~ 0 + cell + date + run + order, data = d)
colnames(X) <- sub("^cell", "", colnames(X))
X <- cbind(X, rib = d$rib)
q <- qr(X)
X <- X[, q$pivot[seq_len(q$rank)]]

fit <- rlm(X, d$log2_raffinose, method = "M", psi = psi.bisquare, init = "ls", maxit = 200)
b <- coef(fit); names(b) <- colnames(X)
V <- vcov(fit)
df <- nrow(d) - ncol(X)

contrast <- function(w_list) {
  w <- setNames(numeric(ncol(X)), colnames(X))
  for (n in names(w_list)) w[n] <- w_list[[n]]
  est <- sum(w * b)
  se <- sqrt(drop(t(w) %*% V %*% w))
  tval <- est / se
  c(log2_change = est, ci_low = est - qt(.975, df) * se, ci_high = est + qt(.975, df) * se,
    t = tval, df = df, p_two_sided = 2 * pt(-abs(tval), df), p_one_sided_greater = pt(-tval, df))
}
res <- rbind(
  dawn_onset = contrast(list(DAWN_T27 = 1, CONT_T27 = -1)),
  dusk_onset = contrast(list(DUSK_T39 = 1, CONT_T39 = -1)),
  dusk_minus_dawn = contrast(list(DUSK_T39 = 1, CONT_T39 = -1, DAWN_T27 = -1, CONT_T27 = 1))
)
out <- data.frame(contrast = rownames(res), res, percent_change = (2^res[, "log2_change"] - 1) * 100,
                  row.names = NULL, check.names = FALSE)
write.csv(out, file.path(figdir, "Figure4_raffinose_model.csv"), row.names = FALSE)
print(format(out, digits = 3))
