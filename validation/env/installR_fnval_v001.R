# fractal-noise validation, R arm. Installs missing packages, then writes the lock (installed versions).
# Run once: Rscript installR_fnval_v001.R
pk  <- c("fracdiff", "longmemo", "arfima", "tseries", "nonlinearTseries", "multitaper", "jsonlite")
need <- pk[!vapply(pk, requireNamespace, logical(1), quietly = TRUE)]
lib <- Sys.getenv("R_LIBS_USER"); dir.create(lib, recursive = TRUE, showWarnings = FALSE); .libPaths(c(lib, .libPaths()))
if (length(need)) install.packages(need, lib = lib, repos = "https://cloud.r-project.org")
bad <- pk[!vapply(pk, requireNamespace, logical(1), quietly = TRUE)]
if (length(bad)) stop("Install FAILED for: ", paste(bad, collapse = ", "))
v <- vapply(pk, function(p) as.character(packageVersion(p)), character(1))
write.csv(data.frame(package = pk, version = v), "r_lock_v001.csv", row.names = FALSE)
cat(sprintf("R %s\n", getRversion())); print(v)
