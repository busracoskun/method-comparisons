# Cross-tabulate diagnosis and modality-availability categories in the prepared cohorts, and
# compute the Rand and adjusted Rand indices between the two partitions (the reference lines
# in Figures 2 and 3). Also verifies that the saved labels match the raw input files

# Paths
out <- Sys.getenv("DIAGNOSTIC_OUTPUT_DIR", "private_outputs/diagnostics/diagnosis_missingness")
root <- Sys.getenv("BENCH_PREPARED_DIR", "private_outputs")
adni <- Sys.getenv("ADNI_INPUT_DIR", "data/adni")
tcga <- Sys.getenv("TCGA_INPUT_DIR", "data/tcga")
dir.create(out, recursive = TRUE, showWarnings = FALSE)

# Path of an input file, or its .gz version when only that exists
input <- function(directory, name) {
  p <- file.path(directory, name)
  if (!file.exists(p) && file.exists(paste0(p, ".gz"))) paste0(p, ".gz") else p
}

# Patient IDs (first column) of an input file, reading no other columns
first_ids <- function(p) {
  hdr <- names(read.csv(p, nrows = 0, check.names = FALSE))
  as.character(
    read.csv(p, colClasses = c("character", rep("NULL", length(hdr) - 1)), check.names = FALSE)[[1]]
  )
}

choose2 <- function(x) x * (x - 1) / 2

# Rand index and adjusted Rand index between partitions y and z, from pair counts
calc <- function(y, z) {
  t <- table(y, z)
  n <- sum(t)
  both <- sum(choose2(t))
  yd <- sum(choose2(rowSums(t)))
  zd <- sum(choose2(colSums(t)))
  pairs <- choose2(n)
  ri <- (pairs - yd - zd + 2 * both) / pairs
  ex <- yd * zd / pairs
  den <- (yd + zd) / 2 - ex
  c(
    n = n, diagnosis_categories = nrow(t), missingness_categories = ncol(t), pairs = pairs,
    same_both = both, different_both = pairs - yd - zd + both, RI = ri,
    ARI = if (den == 0) NA_real_ else (both - ex) / den
  )
}

rows <- list()
counts <- list()
defs <- list()
checks <- list()
for (ds in c("adni", "tcga")) {
  p <- if (ds == "adni") {
    file.path(root, ds, "variables/analysis_variables.rds")
  } else {
    file.path(root, ds, "variables/tcga_bootstrap_analysis_variables.rds")
  }

  d <- readRDS(p)

  # Saved labels for the analysis cohort
  ids <- rownames(d$X_na)
  y <- as.character(d$diag_labels[ids])
  z <- as.character(d$missingness_label[ids])
  stopifnot(
    identical(names(d$diag_labels), ids), identical(names(d$missingness_label), ids),
    !anyDuplicated(ids)
  )

  # Patients expected from the raw input files
  if (ds == "adni") {
    available <- list(
      gene = first_ids(input(adni, "rnaseq.csv")), meta = first_ids(input(adni, "meta.csv"))
    )
    clin <- read.csv(input(adni, "clin.csv"), check.names = FALSE, colClasses = "character")
    rownames(clin) <- clin$RID
    expected <- intersect(Reduce(union, available), rownames(clin))
    rawy <- clin[ids, "diagnosis"]
  } else {
    available <- list(
      mrna = first_ids(input(tcga, "irsem_fs.csv")),
      methylation = first_ids(input(tcga, "imeth_fs.csv")),
      rppa = first_ids(input(tcga, "irppa_fs.csv"))
    )
    clin <- read.csv(input(tcga, "clin_common.csv"), check.names = FALSE, colClasses = "character")
    rownames(clin) <- clin[[1]]
    expected <- intersect(Reduce(union, available), rownames(clin))
    expected <- expected[
      !is.na(clin[expected, "disease_code"]) & nzchar(clin[expected, "disease_code"])
    ]
    rawy <- clin[ids, "disease_code"]
  }

  stopifnot(setequal(ids, expected), identical(y, rawy))

  # Availability pattern of each patient
  pattern <- vapply(ids, function(id) {
    missing <- names(available)[!vapply(available, function(v) id %in% v, logical(1))]
    if (length(missing) == 0) "CC" else paste0("missing_", paste(missing, collapse = "+"))
  }, character(1))
  # A modality counts as present when the patient is in its input file; check that this matches
  # the observed blocks of the aligned matrix
  for (mod in names(available)) {
    cols <- which(d$layer_vec == mod)
    observed <- rowSums(!is.na(d$X_na[, cols, drop = FALSE])) > 0
    stopifnot(identical(unname(observed), unname(ids %in% available[[mod]])))
  }

  # Each saved missingness label must map to one availability pattern
  for (code in unique(z)) {
    mapped <- unique(pattern[z == code])
    stopifnot(length(mapped) == 1)
    defs[[length(defs) + 1]] <- data.frame(
      dataset = ds, saved_label = code, missingness_definition = mapped,
      n_saved_cohort = sum(z == code), n_diagnosis_observed = sum(z == code & !is.na(y)),
      stringsAsFactors = FALSE
    )
  }

  # Rand indices and counts for each cohort definition
  keep <- !is.na(y) & nzchar(y)
  for (scope in c(
    "clinically_labelled_full", "clinically_labelled_complete_cases",
    "historical_saved_all_unknown_separate"
  )) {
    use <- if (scope == "clinically_labelled_full") {
      keep
    } else if (scope == "clinically_labelled_complete_cases") {
      keep & ids %in% rownames(d$X_cc)
    } else {
      rep(TRUE, length(ids))
    }
    yy <- y[use]
    yy[is.na(yy) | !nzchar(yy)] <- "UNKNOWN_DIAGNOSIS"
    zz <- z[use]
    rows[[length(rows) + 1]] <- data.frame(
      dataset = ds, scope = scope, as.list(calc(yy, zz)), stringsAsFactors = FALSE
    )
    t <- as.data.frame(table(diagnosis = toupper(yy), missingness = zz))
    t$dataset <- ds
    t$scope <- scope
    counts[[length(counts) + 1]] <- t
  }

  checks[[length(checks) + 1]] <- data.frame(
    dataset = ds, n_saved = length(ids), n_clinically_labelled = sum(keep), n_unknown = sum(!keep),
    n_saved_complete_cases = nrow(d$X_cc), input_ID_sets_match = TRUE,
    diagnoses_match_clinical_input = TRUE, modality_labels_one_to_one_with_availability = TRUE,
    aligned_block_presence_verified = TRUE
  )
}

write.csv(
  do.call(rbind, rows), file.path(out, "diagnosis_missingness_indices.csv"), row.names = FALSE
)
write.csv(
  do.call(rbind, counts), file.path(out, "diagnosis_by_missingness_counts.csv"), row.names = FALSE
)
write.csv(
  do.call(rbind, defs), file.path(out, "missingness_category_definitions.csv"), row.names = FALSE
)
write.csv(do.call(rbind, checks), file.path(out, "cohort_validation.csv"), row.names = FALSE)
