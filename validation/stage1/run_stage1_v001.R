# Stage-1 estimator runner, R arm (v001). Hub-neutral: reads manifest + series CSVs, writes long-format estimates.
# Usage: Rscript run_stage1_v001.R <manifest.csv> <series_dir> <out.csv> <f_lo_cyc> <f_hi_cyc>
# Likelihood estimators run on the k-times differenced series (ORACLE k from the manifest; see plan open decision 6):
# alpha_hat = alpha_hat(diff) + 2k. Frequencies for multitaper in cycles/sample.
a <- commandArgs(TRUE); man <- read.csv(a[1], stringsAsFactors = FALSE); flo <- as.numeric(a[4]); fhi <- as.numeric(a[5])
`%||%` <- function(a, b) if (is.null(a)) b else a
slope <- function(f, p) { m <- f >= flo & f <= fhi & p > 0; -unname(coef(lm(log10(p[m]) ~ log10(f[m])))[2]) }
EST <- list(
  longmemo_whittle_oraclek = function(x, k) { d <- if (k > 0) diff(x, differences = k) else x
                               2 * longmemo::WhittleEst(d, model = "fGn")$coefficients[1, 1] - 1 + 2 * k },
  fracdiff_mle_oraclek = function(x, k) { d <- if (k > 0) diff(x, differences = k) else x
                               2 * fracdiff::fracdiff(d, nar = 0, nma = 0)$d + 2 * k },
  arfima_mle_oraclek = function(x, k) { d <- if (k > 0) diff(x, differences = k) else x
                               2 * arfima::arfima(d, order = c(0, 0, 0), quiet = TRUE)$modes[[1]]$dfrac + 2 * k },
  multitaper_slope = function(x, k) { s <- multitaper::spec.mtm(ts(x), plot = FALSE); slope(s$freq, s$spec) },
  nlts_dfa1 = function(x, k) { d <- nonlinearTseries::dfa(x, window.size.range = c(10, 400), npoints = 15, do.plot = FALSE)
                               2 * nonlinearTseries::estimate(d, do.plot = FALSE) - 1 })
rows <- list()
for (i in seq_len(nrow(man))) {
  X <- as.matrix(read.csv(file.path(a[2], man$series_file[i]), header = FALSE)); k <- man$k[i]
  for (r in seq_len(ncol(X))) for (nm in names(EST)) {
    t0 <- proc.time()[["elapsed"]]
    v <- tryCatch(suppressMessages(suppressWarnings(as.numeric(EST[[nm]](X[, r], k)))), error = function(e) structure(NaN, msg = conditionMessage(e)))
    st <- if (is.finite(v)) "ok" else paste("FAIL", substr(attr(v, "msg") %||% "non-finite", 1, 120))
    rows[[length(rows) + 1]] <- data.frame(cell_id = man$cell_id[i], rep = r, estimator = nm, language = "R",
      alpha_hat = signif(v, 6), seconds = signif(proc.time()[["elapsed"]] - t0, 4), status = st)
  }
}
write.csv(do.call(rbind, rows), a[3], row.names = FALSE)
