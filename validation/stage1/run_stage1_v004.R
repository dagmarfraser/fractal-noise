# Stage-1 estimator runner, R arm (v004). v004 vs v003: base_dfa4/base_dfa5 (Kiyono 2015 ceiling alpha < 2m+1) and subset
# 'dfa45' (only those two; supplementary pass on an existing run). Earlier header follows.
# Stage-1 estimator runner, R arm (v003). Hub-neutral: reads manifest + series CSVs, writes long-format estimates.
# Usage: Rscript run_stage1_v004.R <manifest.csv> <series_dir> <out.csv> <fs> <band_label> <cores> <all|spectral|dfa45>
# Band per cell from manifest columns f_lo_<band_label>, f_hi_<band_label> (Hz; converted to cycles/sample here).
# 'spectral' runs only multitaper_slope (the one band-dependent R estimator). Output adds a band column.
# v003 vs v002: per-cell band, band column, spectral-only mode.
# Likelihood estimators run on the k-times differenced series, alpha_hat = alpha_hat(diff) + 2k, in two arms:
#   *_oraclek : k from the manifest (generator truth)
#   *_estk    : k chosen from the data: pre-estimate alpha with DFA-3 (valid to alpha ~7), k = floor((alpha~ + 1)/2),
#               clamped to 0..3. Chosen k is in k_used. (A Whittle-H "closest to 0.5" rule was tried and rejected:
#               longmemo bounds H to [0.1, 0.99] and once-differenced k=2 series with small H return interior H ~0.89,
#               so it over/under-picked k at alpha 2.75/3.25; see kselectDiag_v001.R.)
# v002 vs v001: estimated-k arm, DFA orders 1-3 (base R, non-overlapping boxes), mclapply over series, k_used column.
a <- commandArgs(TRUE); man <- read.csv(a[1], stringsAsFactors = FALSE); fs <- as.numeric(a[4]); blab <- a[5]
cores <- as.integer(a[6]); subset <- a[7]
if (!subset %in% c("all", "spectral", "dfa45")) stop("subset must be all|spectral|dfa45, got ", subset)
for (cn in paste0(c("f_lo_", "f_hi_"), blab)) if (!cn %in% names(man)) stop("manifest lacks column ", cn)
`%||%` <- function(a, b) if (is.null(a)) b else a
dk <- function(x, k) if (k > 0) diff(x, differences = k) else x
slope <- function(f, p, b) { m <- f >= b[1] & f <= b[2] & p > 0; -unname(coef(lm(log10(p[m]) ~ log10(f[m])))[2]) }
whH <- function(d) longmemo::WhittleEst(d, model = "fGn")$coefficients[1, 1]
dfaN <- function(x, order) { y <- cumsum(x - mean(x)); n <- length(y)
  s <- unique(round(exp(seq(log(10), log(n %/% 10), length.out = 15))))
  Fs <- sapply(s, function(w) { nb <- n %/% w; Y <- matrix(y[seq_len(nb * w)], nrow = w)   # one box per column
    q <- qr(cbind(1, poly(seq_len(w), order))); sqrt(mean(qr.resid(q, Y)^2)) })
  2 * unname(coef(lm(log(Fs) ~ log(s)))[2]) - 1 }
selectK <- function(x) min(max(floor((dfaN(x, 3) + 1) / 2), 0), 3)
# Each returns alpha of the differenced series; a parameter on (within 1e-3 of) its optimiser bound is an error,
# not an estimate (longmemo H in [0.1, 0.99]; fracdiff d in [0, 0.5]; arfima d in (-0.5, 0.5)), so it is counted as FAIL.
onBound <- function(p, lo, hi, nm) { if (p - lo < 1e-3 || hi - p < 1e-3) stop(sprintf("%s on or outside support bound, censored (%.4f vs [%g, %g])", nm, p, lo, hi)); p }
LIK <- list(longmemo_whittle = function(d) 2 * onBound(whH(d), 0.1, 0.99, "H") - 1,
            fracdiff_mle     = function(d) 2 * onBound(fracdiff::fracdiff(d, nar = 0, nma = 0)$d, 0, 0.5, "d"),
            arfima_mle       = function(d) 2 * onBound(arfima::arfima(d, order = c(0, 0, 0), quiet = TRUE)$modes[[1]]$dfrac, -0.5, 0.5, "d"))
SPC <- list(multitaper_slope = function(x, b) { s <- multitaper::spec.mtm(ts(x), plot = FALSE); slope(s$freq, s$spec, b) })
OTH <- list(nlts_dfa1 = function(x) { d <- nonlinearTseries::dfa(x, window.size.range = c(10, length(x) %/% 10), npoints = 15, do.plot = FALSE)
                                      2 * nonlinearTseries::estimate(d, do.plot = FALSE) - 1 },
            base_dfa1 = function(x) dfaN(x, 1), base_dfa2 = function(x) dfaN(x, 2), base_dfa3 = function(x) dfaN(x, 3),
            base_dfa4 = function(x) dfaN(x, 4), base_dfa5 = function(x) dfaN(x, 5))
run1 <- function(f) { t0 <- proc.time()[["elapsed"]]
  v <- tryCatch(suppressMessages(suppressWarnings(as.numeric(f()))), error = function(e) structure(NaN, msg = conditionMessage(e)))
  st <- if (length(v) == 1 && is.finite(v)) "ok" else paste("FAIL", substr(attr(v, "msg") %||% "non-finite", 1, 120))
  list(v = signif(v[1], 6), s = signif(proc.time()[["elapsed"]] - t0, 4), st = st) }
one <- function(cid, r, x, k, b) {
  if (subset == "dfa45") { out <- list()
    for (nm in c("base_dfa4", "base_dfa5")) { z <- run1(function() OTH[[nm]](x))
      out[[length(out) + 1]] <- data.frame(cell_id = cid, rep = r, estimator = nm, language = "R", alpha_hat = z$v,
        seconds = z$s, status = z$st, k_used = NA, band = "none") }
    return(do.call(rbind, out)) }
  out <- list(); add <- function(nm, z, ku, bl = "none") out[[length(out) + 1]] <<- data.frame(cell_id = cid, rep = r, estimator = nm,
    language = "R", alpha_hat = z$v, seconds = z$s, status = z$st, k_used = ku, band = bl)
  for (nm in names(SPC)) add(nm, run1(function() SPC[[nm]](x, b / fs)), NA, blab)
  if (subset == "spectral") return(do.call(rbind, out))
  ks <- run1(function() selectK(x)); add("kselect_dfa3", ks, if (ks$st == "ok") ks$v else NA)
  for (nm in names(LIK)) {
    add(paste0(nm, "_oraclek"), run1(function() LIK[[nm]](dk(x, k)) + 2 * k), k)
    if (ks$st == "ok") add(paste0(nm, "_estk"), run1(function() LIK[[nm]](dk(x, ks$v)) + 2 * ks$v), ks$v)
    else add(paste0(nm, "_estk"), list(v = NaN, s = NA, st = paste("FAIL k selection:", ks$st)), NA)
  }
  for (nm in names(OTH)) add(nm, run1(function() OTH[[nm]](x)), NA)
  do.call(rbind, out) }
jobs <- list()
for (i in seq_len(nrow(man))) { X <- as.matrix(read.csv(file.path(a[2], man$series_file[i]), header = FALSE))
  b <- c(man[[paste0("f_lo_", blab)]][i], man[[paste0("f_hi_", blab)]][i])
  for (r in seq_len(ncol(X))) jobs[[length(jobs) + 1]] <- list(cid = man$cell_id[i], r = r, x = X[, r], k = man$k[i], b = b) }
res <- parallel::mclapply(jobs, function(j) one(j$cid, j$r, j$x, j$k, j$b), mc.cores = cores)
bad <- vapply(res, function(z) inherits(z, "try-error") || !is.data.frame(z), logical(1))
if (any(bad)) stop(sprintf("%d worker job(s) crashed, e.g.: %s", sum(bad), as.character(res[[which(bad)[1]]])))
write.csv(do.call(rbind, res), a[3], row.names = FALSE)
