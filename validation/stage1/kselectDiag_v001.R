# Diagnostic: Whittle fGn H-hat on diff^j(x), j = 0..4, for each series of chosen pilot cells (median, 10%, 90%).
a <- commandArgs(TRUE); sdir <- a[1]
for (cf in a[-1]) { X <- as.matrix(read.csv(file.path(sdir, cf), header = FALSE))
  H <- sapply(seq_len(ncol(X)), function(r) sapply(0:4, function(j) { d <- if (j > 0) diff(X[, r], differences = j) else X[, r]
    tryCatch(longmemo::WhittleEst(d, model = "fGn")$coefficients[1, 1], error = function(e) NA) }))
  cat(cf, "\n"); print(round(t(apply(H, 1, quantile, c(.1, .5, .9), na.rm = TRUE)), 3)) }
