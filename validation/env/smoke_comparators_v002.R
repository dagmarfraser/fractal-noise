# Multi-seed smoke test (v002): every R comparator on K shared exact-fGn series (alpha = 0.5).
# Usage: Rscript smoke_comparators_v002.R <series.csv (N x K)> <seeds.csv (K)>. Prints RES|tool|k|value|status.
# NB tseries::surrogate(fft = TRUE) is NOT variance-preserving (sd ~2900x input at N = 4096): slopes only.
args <- commandArgs(TRUE); X <- as.matrix(read.csv(args[1], header = FALSE)); seeds <- scan(args[2], quiet = TRUE)
`%||%` <- function(a, b) if (is.null(a)) b else a
run <- function(name, expr, k) {
  v <- tryCatch(suppressMessages(suppressWarnings(as.numeric(expr))), error = function(e) structure(NA, msg = conditionMessage(e)))
  st <- if (is.na(v)) paste("FAIL", substr(attr(v, "msg") %||% "NA result", 1, 150)) else "ok"
  cat(sprintf("RES|%s|%d|%.6g|%s\n", name, k, v, st))
}
slope <- function(f, p) { m <- f > 0.005 & f < 0.2; -coef(lm(log10(p[m]) ~ log10(f[m])))[2] }
for (k in seq_len(ncol(X))) {
  x <- X[, k]; set.seed(seeds[k] %% .Machine$integer.max)
  run("longmemo_whittle_alpha", 2 * longmemo::WhittleEst(x, model = "fGn")$coefficients[1, 1] - 1, k)
  run("fracdiff_mle_alpha", 2 * fracdiff::fracdiff(x, nar = 0, nma = 0)$d, k)
  run("arfima_mle_alpha", 2 * arfima::arfima(x, order = c(0, 0, 0), quiet = TRUE)$modes[[1]]$dfrac, k)
  run("tseries_surrogate_slope", { s <- spec.pgram(as.numeric(tseries::surrogate(x, ns = 1, fft = TRUE)), spans = c(9, 9), plot = FALSE, taper = 0.1); slope(s$freq, s$spec) }, k)
  run("multitaper_slope", { s <- multitaper::spec.mtm(ts(x), plot = FALSE); slope(s$freq, s$spec) }, k)
  run("nonlinearTseries_dfa_alpha", { d <- nonlinearTseries::dfa(x, window.size.range = c(10, 400), npoints = 15, do.plot = FALSE)
                                      2 * nonlinearTseries::estimate(d, do.plot = FALSE) - 1 }, k)
}
