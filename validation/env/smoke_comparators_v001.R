# One real call per R comparator on the shared exact-fGn series (alpha = 0.5, H = 0.75).
# Usage: Rscript smoke_comparators_v001.R <series.csv>. Prints RES|tool|value|status lines.
args <- commandArgs(TRUE); x <- scan(args[1], quiet = TRUE)
run <- function(name, expr) {
  v <- tryCatch(suppressMessages(suppressWarnings(as.numeric(expr))), error = function(e) structure(NA, msg = conditionMessage(e)))
  st <- if (is.na(v)) paste("FAIL", substr(attr(v, "msg") %||% "NA result", 1, 150)) else "ok"
  cat(sprintf("RES|%s|%.6g|%s\n", name, v, st))
}
`%||%` <- function(a, b) if (is.null(a)) b else a
run("longmemo_whittle_alpha", { w <- longmemo::WhittleEst(x, model = "fGn"); 2 * w$coefficients[1, 1] - 1 })
run("fracdiff_mle_alpha",     2 * fracdiff::fracdiff(x, nar = 0, nma = 0)$d)
run("arfima_mle_alpha",       { a <- arfima::arfima(x, order = c(0, 0, 0), quiet = TRUE); 2 * a$modes[[1]]$dfrac })
# NB tseries::surrogate(fft = TRUE) output is NOT variance-preserving (sd ~ 2900x input at N = 4096,
# unnormalised inverse FFT): score the scale-free spectral slope instead, and record the quirk.
run("tseries_surrogate_slope", { y <- as.numeric(tseries::surrogate(x, ns = 1, fft = TRUE)); s <- spec.pgram(y, spans = c(9, 9), plot = FALSE, taper = 0.1)
                                 m <- s$freq > 0.005 & s$freq < 0.2; -coef(lm(log10(s$spec[m]) ~ log10(s$freq[m])))[2] })
run("multitaper_slope", { s <- multitaper::spec.mtm(ts(x), plot = FALSE); m <- s$freq > 0.005 & s$freq < 0.2
                          -coef(lm(log10(s$spec[m]) ~ log10(s$freq[m])))[2] })
run("nonlinearTseries_dfa_alpha", { d <- nonlinearTseries::dfa(x, window.size.range = c(10, 400), npoints = 15, do.plot = FALSE)
                                    2 * nonlinearTseries::estimate(d, do.plot = FALSE) - 1 })
