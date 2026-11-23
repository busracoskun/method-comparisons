# Reproduce the GIPCA failures on ADNI bootstrap replicates 3, 8, and 865 and report the
# non-finite SVD input that stops the fit. See README.md

diag_out <- Sys.getenv("DIAGNOSTIC_OUTPUT_DIR", "private_outputs/diagnostics/gipca_bidifac")
for (part in c("private", "public")) {
  dir.create(file.path(diag_out, part), recursive = TRUE, showWarnings = FALSE)
}

root <- Sys.getenv("DIAGNOSTIC_OUTPUT_DIR", "private_outputs/diagnostics/gipca_bidifac")
for (f in list.files(
  Sys.getenv("GIPCA_SOURCE_DIR", "external_sources/gipca/Functions"),
  full.names = TRUE, pattern = "\\.R$"
)) {
  source(f)
}

v <- readRDS(file.path(Sys.getenv("BENCH_ARCHIVE_DIR"), "variables/analysis_variables.rds"))
options(warn = 1)
# Print aggregate matrix properties and the call stack at the failure, never values or IDs
trace(svd, tracer = quote({
  if (any(!is.finite(x))) {
    cat("NONFINITE SVD INPUT", dim(x), "NA", sum(is.na(x)), "Inf", sum(is.infinite(x)), "\n")
    print(sys.calls())
    for (e in sys.frames()) {
      if (exists("l", e, inherits = FALSE) && exists("U.joint", e, inherits = FALSE)) {
        cat("GIPCA iteration", get("l", e), "\n")
        for (nm in c("U.joint", "V_new.joint", "Theta.joint")) {
          if (exists(nm, e, inherits = FALSE)) {
            a <- get(nm, e)
            cat(nm, dim(a), "NA", sum(is.na(a)), "\n")
          }
        }
      }
    }
  }
}), print = FALSE)

# Refit GIPCA on each failing replicate
for (b in c(3L, 8L, 865L)) {
  set.seed(100 + b)
  ix <- sample(rownames(v$X_stand), nrow(v$X_stand), replace = TRUE)
  X <- v$X_stand[ix, , drop = FALSE]
  blocks <- split(colnames(X), v$layer_vec[colnames(X)])
  dat <- lapply(blocks, function(cc) as.matrix(X[, cc, drop = FALSE]))
  cat("\nBOOT", b, "\n")
  # In the benchmark, earlier NA methods consume random numbers first, so this seed reproduces
  # the bootstrap sample but not the exact original RNG state at the GIPCA fit
  tryCatch(withCallingHandlers(
    {
      fit <- EPCAJIVEMissbio(
        dat,
        rankj = 3, ranka = rep(0L, length(dat)), D = vapply(dat, ncol, integer(1)),
        family = rep("gaussian", length(dat)), tol = .1, max.iter = 500, n.size = NULL, lambda = 0
      )
      cat("RETURNED\n")
    },
    error = function(e) {
      cat("ERROR STACK\n")
      print(sys.calls())
    }
  ), error = function(e) cat("ERROR", conditionMessage(e), "\n"))
}

untrace(svd)
