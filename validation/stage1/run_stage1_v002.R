# Stage-1 estimator runner, R arm (v002). Hub-neutral: reads manifest + series CSVs, writes long-format estimates.
# Usage: Rscript run_stage1_v002.R <manifest.csv> <series_dir> <out.csv> <f_lo_cyc> <f_hi_cyc> [cores]
# Likelihood estimators run on the k-times differenced series, alpha_hat = alpha_hat(diff) + 2k, in two arms:
#   *_oraclek : k from the manifest (generator truth)
#   *_estk    : k chosen from the data: Whittle fGn fit on diff^j(x), j = 0..4; pick the j whose H-hat is closest
#               to 0.5 (under-differencing pushes H-hat to 1, over-differencing to 0). Chosen k is in k_used.
# v002 vs v001: estimated-k arm, DFA orders 1-3 (base R, non-overlapping boxes), mclapply over series, k_used column.
a <- commandArgs(TRUE); man <- read.csv(a[1], stringsAsFactors = FALSE); flo <- as.numeric(a[4]); fhi <- as.numeric(a[5])
cores <- if (length(a) >= 6) as.integer(a[6]) else 1L
`%||%` <- function(a, b) if (is.null(a)) b else a
dk <- function(x, k) if (k > 0) diff(x, differences = k) else x
slope <- function(f, p) { m <- f >= flo & f <= fhi & p > 0; -unname(coef(lm(log10(p[m]) ~ log10(f[m])))[2]) }
whH <- function(d) longmemo::WhittleEst(d, model = "fGn")$coefficients[1, 1]
selectK <- function(x) { h <- sapply(0:4, function(j) tryCatch(whH(dk(x, j)), error = function(e) NA_real_))
  if (all(is.na(h))) stop("selectK: Whittle failed at every k"); which.min(abs(h - 0.5)) - 1 }
dfaN <- function(x, order) { y <- cumsum(x - mean(x)); n <- length(y)
  s <- unique(round(exp(seq(log(10), log(n %/% 10), length.out = 15))))
  Fs <- sapply(s, function(w) { nb <- n %/% w; Y <- matrix(y[seq_len(nb * w)], nrow = w)   # one box per column
    q <- qr(cbind(1, poly(seq_len(w), order))); sqrt(mean(qr.resid(q, Y)^2)) })
  2 * unname(coef(lm(log(Fs) ~ log(s)))[2]) - 1 }
LIK <- list(longmemo_whittle = function(d) 2 * whH(d) - 1,
            fracdiff_mle     = function(d) 2 * fracdiff::fracdiff(d, nar = 0, nma = 0)$d,
            arfima_mle       = function(d) 2 * arfima::arfima(d, order = c(0, 0, 0), quiet = TRUE)$modes[[1]]$dfrac)
OTH <- list(multitaper_slope = function(x) { s <- multitaper::spec.mtm(ts(x), plot = FALSE); slope(s$freq, s$spec) },
            nlts_dfa1 = function(x) { d <- nonlinearTseries::dfa(x, window.size.range = c(10, length(x) %/% 10), npoints = 15, do.plot = FALSE)
                                      2 * nonlinearTseries::estimate(d, do.plot = FALSE) - 1 },
            base_dfa1 = function(x) dfaN(x, 1), base_dfa2 = function(x) dfaN(x, 2), base_dfa3 = function(x) dfaN(x, 3))
run1 <- function(f) { t0 <- proc.time()[["elapsed"]]
  v <- tryCatch(suppressMessages(suppressWarnings(as.numeric(f()))), error = function(e) structure(NaN, msg = conditionMessage(e)))
  st <- if (length(v) == 1 && is.finite(v)) "ok" else paste("FAIL", substr(attr(v, "msg") %||% "non-finite", 1, 120))
  list(v = signif(v[1], 6), s = signif(proc.time()[["elapsed"]] - t0, 4), st = st) }
one <- function(cid, r, x, k) {
  out <- list(); add <- function(nm, z, ku) out[[length(out) + 1]] <<- data.frame(cell_id = cid, rep = r, estimator = nm,
    language = "R", alpha_hat = z$v, seconds = z$s, status = z$st, k_used = ku)
  ks <- run1(function() selectK(x)); add("kselect_whittle", ks, if (ks$st == "ok") ks$v else NA)
  for (nm in names(LIK)) {
    add(paste0(nm, "_oraclek"), run1(function() LIK[[nm]](dk(x, k)) + 2 * k), k)
    if (ks$st == "ok") add(paste0(nm, "_estk"), run1(function() LIK[[nm]](dk(x, ks$v)) + 2 * ks$v), ks$v)
    else add(paste0(nm, "_estk"), list(v = NaN, s = NA, st = paste("FAIL k selection:", ks$st)), NA)
  }
  for (nm in names(OTH)) add(nm, run1(function() OTH[[nm]](x)), NA)
  do.call(rbind, out) }
jobs <- list()
for (i in seq_len(nrow(man))) { X <- as.matrix(read.csv(file.path(a[2], man$series_file[i]), header = FALSE))
  for (r in seq_len(ncol(X))) jobs[[length(jobs) + 1]] <- list(cid = man$cell_id[i], r = r, x = X[, r], k = man$k[i]) }
res <- parallel::mclapply(jobs, function(j) one(j$cid, j$r, j$x, j$k), mc.cores = cores)
bad <- vapply(res, function(z) inherits(z, "try-error") || !is.data.frame(z), logical(1))
if (any(bad)) stop(sprintf("%d worker job(s) crashed, e.g.: %s", sum(bad), as.character(res[[which(bad)[1]]])))
write.csv(do.call(rbind, res), a[3], row.names = FALSE)
