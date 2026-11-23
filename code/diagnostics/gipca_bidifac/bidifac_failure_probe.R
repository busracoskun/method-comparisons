# Trace the SVD calls made by BIDIFAC on TCGA bootstrap replicate 4 to locate the failing call
# (internal softSVD versus the final score matrix). See README.md for inputs and outputs

diag_out <- Sys.getenv("DIAGNOSTIC_OUTPUT_DIR", "private_outputs/diagnostics/gipca_bidifac")
for (part in c("private", "public")) {
  dir.create(file.path(diag_out, part), recursive = TRUE, showWarnings = FALSE)
}

root <- Sys.getenv("DIAGNOSTIC_OUTPUT_DIR", "private_outputs/diagnostics/gipca_bidifac")
p <- Sys.getenv("BENCH_ARCHIVE_DIR")
# Evaluate only the vendor function definitions, skipping its library() calls
e <- new.env(parent = globalenv())
for (expr in parse(Sys.getenv("BIDIFAC_FILE", file.path(p, "vendor/bidifac/BIDIFAC.R")))) {
  if (!is.call(expr) || !identical(expr[[1]], as.name("library"))) eval(expr, e)
}

v <- readRDS(file.path(p, "variables/tcga_bootstrap_analysis_variables.rds"))

# Record properties of each SVD input and where it was called from
trace(svd, tracer = quote({
  if (exists("bidifac_diagnostic", .GlobalEnv)) {
    bidifac_diagnostic$calls <- bidifac_diagnostic$calls + 1L
    bidifac_diagnostic$last <- list(
      dim = dim(x),
      nonfinite = sum(!is.finite(x)),
      max_abs = max(abs(x)),
      context = if (any(vapply(
        sys.calls(), function(c) grepl("softSVD", paste(deparse(c), collapse = "")), logical(1)
      ))) "internal_softSVD" else "final_Smat"
    )
  }
}), print = FALSE)

# Refit BIDIFAC on the failing replicate
for (b in c(4L)) {
  set.seed(100 + b)
  ix <- sample(rownames(v$X_stand), nrow(v$X_stand), TRUE)
  X <- v$X_stand[ix, , drop = FALSE]
  bl <- unique(v$layer_vec[colnames(X)])
  dat <- matrix(vector("list", length(bl)), 1, length(bl))
  for (j in seq_along(bl)) dat[[1, j]] <- X[, v$layer_vec[colnames(X)] == bl[j], drop = FALSE]

  bidifac_diagnostic <- new.env()
  bidifac_diagnostic$calls <- 0L
  cat("BOOT", b, "\n")
  tryCatch(withCallingHandlers(
    {
      fit <- e$impute.BIDIFAC(dat, rmt = FALSE, sigma = matrix(1, 1, length(bl)), pbar = FALSE)
      Smat <- e$data.rearrange(fit$S, rmt = FALSE)$out
      sv <- svd(Smat)
      cat("RETURNED\n")
    },
    error = function(err) {
      print(bidifac_diagnostic$last)
      cat("SVD calls", bidifac_diagnostic$calls, "\n")
      print(sys.calls())
    }
  ), error = function(err) cat("ERROR", conditionMessage(err), "\n"))
}

untrace(svd)
