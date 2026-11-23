# TCGA modality availability by disease group (from input-file membership) and the number of
# features per modality retained in saved CC_SF/MI_SF fits and screening reconstructions
# Reads saved objects only; SIS is never rerun and no model is fit

# Paths and prepared variables
prepared <- Sys.getenv("BENCH_PREPARED_DIR", "private_outputs")
out <- Sys.getenv("DIAGNOSTIC_OUTPUT_DIR", "private_outputs/diagnostics/tcga_screening")
fits <- Sys.getenv("TCGA_R_FITS_DIR", file.path(prepared, "tcga", "fits"))
replay <- Sys.getenv("SCREENING_REPLAY_DIR", "private_outputs/diagnostics/gpcca/screening_replay")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
v <- readRDS(file.path(prepared, "tcga/variables/tcga_bootstrap_analysis_variables.rds"))
# Availability of each modality for the 841 saved analysis IDs, from input-file membership
base <- Sys.getenv("TCGA_INPUT_DIR", "data/tcga")

ids <- rownames(v$X_na)
clin <- read.csv(file.path(base, "clin_common.csv"), colClasses = "character", check.names = FALSE)
rownames(clin) <- clin[[1]]

files <- c(mrna = "irsem_fs.csv", methylation = "imeth_fs.csv", rppa = "irppa_fs.csv")
avail <- lapply(files, function(f) {
  p <- file.path(base, f)
  hdr <- names(read.csv(p, nrows = 0))
  as.character(read.csv(p, colClasses = c("character", rep("NULL", length(hdr) - 1)))[[1]])
})

# Check the saved cohort against the input files
stopifnot(setequal(ids, intersect(Reduce(union, avail), rownames(clin))), length(ids) == 841)
y <- toupper(as.character(v$diag_labels[ids]))
stopifnot(identical(y, toupper(clin[ids, "disease_code"])))
D <- data.frame(
  diagnosis = y, mrna = ids %in% avail$mrna, methylation = ids %in% avail$methylation,
  rppa = ids %in% avail$rppa
)
for (mod in names(avail)) {
  stopifnot(identical(
    unname(rowSums(!is.na(v$X_na[, v$layer_vec == mod, drop = FALSE])) > 0), D[[mod]]
  ))
}

# Patient counts for every availability pattern
grid <- expand.grid(mrna = c(TRUE, FALSE), methylation = c(TRUE, FALSE), rppa = c(TRUE, FALSE))
rows <- list()
for (i in seq_len(nrow(grid))) {
  g <- grid[i, ]
  keep <- D$mrna == g$mrna & D$methylation == g$methylation & D$rppa == g$rppa
  r <- data.frame(
    mrna = g$mrna, methylation = g$methylation, rppa = g$rppa, GBM = sum(keep & y == "GBM"),
    LGG = sum(keep & y == "LGG")
  )
  r$available <- paste(names(files)[as.logical(g)], collapse = " + ")
  r$missing <- paste(names(files)[!as.logical(g)], collapse = " + ")
  if (r$missing == "") r$missing <- "None"
  if (r$available == "") r$available <- "None"
  rows[[i]] <- r
}

t <- do.call(rbind, rows)
t$total <- t$GBM + t$LGG
t$GBM_percent <- 100 * t$GBM / sum(y == "GBM")
t$LGG_percent <- 100 * t$LGG / sum(y == "LGG")
t$cohort_percent <- 100 * t$total / nrow(D)
write.csv(t, file.path(out, "tcga_missingness_by_diagnosis.csv"), row.names = FALSE)

# Features retained in the saved full-cohort fits
ff <- list.files(fits, pattern = "_(CC_SF|MI_SF)_0.*\\.rds$", full.names = TRUE, recursive = TRUE)
fr <- list()
for (p in ff) {
  f <- readRDS(p)
  cols <- f$cols
  if (is.null(cols)) next
  lv <- v$layer_vec[cols]
  stopifnot(!anyNA(lv))
  z <- data.frame(
    engine = "R-compatible saved objects", method = f$model, style = f$style,
    imputation = if (is.null(f$imp)) 0 else f$imp, features = length(cols),
    mrna = sum(lv == "mrna"), methylation = sum(lv == "methylation"), rppa = sum(lv == "rppa"),
    file = basename(p)
  )
  z$feature_set_signature <- paste(sort(cols), collapse = "|")
  fr[[length(fr) + 1]] <- z
}

f <- do.call(rbind, fr)
write.csv(f, file.path(out, "saved_full_fit_feature_counts_internal.csv"), row.names = FALSE)
write.csv(
  data.frame(feature = names(v$layer_vec), modality = as.character(v$layer_vec)),
  file.path(out, "feature_modality_map.csv"), row.names = FALSE
)

# Features retained in previously saved CC_SF screening reconstructions
rr <- list()
for (p in list.files(replay, pattern = "^boot_[0-9]+\\.rds$", full.names = TRUE)) {
  r <- readRDS(p)
  cols <- r$selected
  lv <- v$layer_vec[cols]
  rr[[length(rr) + 1]] <- data.frame(
    boot = r$boot, mrna = sum(lv == "mrna"), methylation = sum(lv == "methylation"),
    rppa = sum(lv == "rppa"), total = length(cols),
    provenance = "previously saved seeded CC_SF reconstruction, not original per-fit metadata"
  )
}
write.csv(
  do.call(rbind, rr), file.path(out, "cc_sf_saved_reconstruction_feature_counts.csv"),
  row.names = FALSE
)
