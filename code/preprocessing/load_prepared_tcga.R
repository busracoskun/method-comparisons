# Load and align the prepared TCGA assay and clinical inputs
# Sourced by load_dataset() into the dataset environment; defines X_all, X_stand, layer_vec,
# feature_name_map, diag_labels, mmse_vec (Karnofsky score), missingness_label, and common_ids

# Input files
TCGA_MRNA_FILE <- file.path(TCGA_INPUT_DIR, "irsem_fs.csv")
TCGA_METHYLATION_FILE <- file.path(TCGA_INPUT_DIR, "imeth_fs.csv")
TCGA_RPPA_FILE <- file.path(TCGA_INPUT_DIR, "irppa_fs.csv")
TCGA_MIRNA_FILE <- file.path(TCGA_INPUT_DIR, "imirna_fs.csv")
TCGA_CLINICAL_FILE <- file.path(TCGA_INPUT_DIR, "clin_common.csv")

required_tcga_files <- c(TCGA_MRNA_FILE, TCGA_METHYLATION_FILE, TCGA_RPPA_FILE, TCGA_CLINICAL_FILE)
missing_tcga_files <- required_tcga_files[!file.exists(required_tcga_files)]
if (length(missing_tcga_files)) {
  stop("Required prepared TCGA file(s) were not found:\n", paste(missing_tcga_files, collapse = "\n"))
}

# HELPERS

# Read one modality and check its IDs and values
read_prepared_omics <- function(file_path, modality) {
  cat("Reading", modality, "from:", file_path, "\n")
  dat <- read.csv(
    file_path,
    header = TRUE, check.names = FALSE, stringsAsFactors = FALSE,
    na.strings = c("", "NA", "N/A", "NaN", "NULL"), quote = "\"", comment.char = ""
  )
  if (ncol(dat) < 2L) {
    stop(modality, " must contain one ID column and at least one feature.")
  }

  sample_ids <- as.character(dat[[1]])
  if (anyNA(sample_ids) || any(!nzchar(sample_ids))) {
    stop(modality, " contains a missing or blank sample ID.")
  }
  if (anyDuplicated(sample_ids)) {
    stop(modality, " contains duplicated sample IDs.")
  }

  feature_data <- dat[, -1, drop = FALSE]
  if (anyDuplicated(colnames(feature_data))) {
    stop(modality, " contains duplicated feature names within its own file.")
  }
  nonnumeric <- names(feature_data)[!vapply(feature_data, is.numeric, logical(1))]
  if (length(nonnumeric)) {
    stop(modality, " contains nonnumeric feature columns: ", paste(nonnumeric, collapse = ", "))
  }

  M <- as.matrix(feature_data)
  storage.mode(M) <- "double"
  rownames(M) <- sample_ids
  if (any(is.infinite(M))) {
    stop(modality, " contains positive or negative infinity.")
  }
  cat("  Samples:", nrow(M), "| Features:", ncol(M), "| Missing cells:", sum(is.na(M)), "\n")
  M
}

# Read the clinical table and check its IDs
read_prepared_clinical <- function(file_path) {
  cat("Reading clinical data from:", file_path, "\n")
  clin <- read.csv(
    file_path,
    header = TRUE, check.names = FALSE, stringsAsFactors = FALSE,
    na.strings = c("", "NA", "N/A", "NaN", "NULL"), quote = "\"", comment.char = ""
  )
  if (ncol(clin) < 2L) {
    stop("clin_common.csv must contain an ID column and clinical variables.")
  }

  sample_ids <- as.character(clin[[1]])
  if (anyNA(sample_ids) || any(!nzchar(sample_ids))) {
    stop("clin_common.csv contains a missing or blank sample ID.")
  }
  if (anyDuplicated(sample_ids)) {
    stop("clin_common.csv contains duplicated sample IDs.")
  }
  rownames(clin) <- sample_ids
  clin
}

# Place a modality on the shared patient index; absent patients become NA rows
align_prepared_matrix <- function(M, ids) {
  out <- matrix(NA_real_, nrow = length(ids), ncol = ncol(M), dimnames = list(ids, colnames(M)))
  shared <- intersect(ids, rownames(M))
  out[shared, ] <- M[shared, , drop = FALSE]
  out
}

# READ INPUTS

# Read the modalities and the clinical table
omics_list <- list(
  mrna = read_prepared_omics(TCGA_MRNA_FILE, "mRNA"),
  methylation = read_prepared_omics(TCGA_METHYLATION_FILE, "methylation"),
  rppa = read_prepared_omics(TCGA_RPPA_FILE, "RPPA")
)
if (file.exists(TCGA_MIRNA_FILE)) {
  omics_list$mirna <- read_prepared_omics(TCGA_MIRNA_FILE, "miRNA")
} else {
  cat("Optional prepared miRNA file not found; continuing without miRNA:\n", TCGA_MIRNA_FILE, "\n")
}

clin <- read_prepared_clinical(TCGA_CLINICAL_FILE)
required_clinical_columns <- c("disease_code", "karnofsky_performance_score", "missingness")
missing_clinical_columns <- setdiff(required_clinical_columns, colnames(clin))
if (length(missing_clinical_columns)) {
  missing_list <- paste(missing_clinical_columns, collapse = ", ")
  stop("clin_common.csv is missing required columns: ", missing_list)
}

# Patients: union of assayed patients with a clinical record and a disease code
ids_by_modality <- lapply(omics_list, rownames)
omics_union_ids <- Reduce(union, ids_by_modality)
ids_score <- intersect(omics_union_ids, rownames(clin))
diagnosis_raw <- as.character(clin[ids_score, "disease_code"])
keep_diagnosis <- !is.na(diagnosis_raw) & nzchar(diagnosis_raw)
if (any(!keep_diagnosis)) {
  cat("Dropping", sum(!keep_diagnosis), "patients without disease_code.\n")
}
ids_score <- ids_score[keep_diagnosis]
diagnosis_raw <- diagnosis_raw[keep_diagnosis]
if (length(ids_score) < 2L) {
  stop("Fewer than two TCGA patients remain after clinical matching.")
}

# COMBINE MODALITIES

# Feature names can repeat across modalities, so they are made unique in memory
# and the original names are kept in feature_name_map
omics_aligned <- lapply(omics_list, align_prepared_matrix, ids = ids_score)
original_feature_names <- unlist(lapply(omics_aligned, colnames), use.names = FALSE)
internal_feature_names <- make.unique(original_feature_names, sep = ".")
feature_modalities <- rep(names(omics_aligned), times = vapply(omics_aligned, ncol, integer(1)))
feature_name_map <- data.frame(
  internal_name = internal_feature_names, original_name = original_feature_names,
  modality = feature_modalities, stringsAsFactors = FALSE
)

cursor <- 1L
for (modality in names(omics_aligned)) {
  width <- ncol(omics_aligned[[modality]])
  idx <- cursor:(cursor + width - 1L)
  colnames(omics_aligned[[modality]]) <- internal_feature_names[idx]
  cursor <- cursor + width
}

X_all <- do.call(cbind, omics_aligned)
storage.mode(X_all) <- "double"

# Standardize each feature using its observed values only
feature_means <- colMeans(X_all, na.rm = TRUE)
feature_sds <- vapply(
  seq_len(ncol(X_all)), function(j) stats::sd(X_all[, j], na.rm = TRUE), numeric(1)
)
names(feature_means) <- names(feature_sds) <- colnames(X_all)
if (any(!is.finite(feature_sds) | feature_sds <= 0)) {
  stop("TCGA contains zero-variance or invalid features before standardization.")
}
X_stand <- sweep(X_all, 2, feature_means, FUN = "-")
X_stand <- sweep(X_stand, 2, feature_sds, FUN = "/")
storage.mode(X_stand) <- "double"

layer_vec <- setNames(feature_modalities, internal_feature_names)

# CLINICAL LABELS

diag_labels <- factor(diagnosis_raw)
names(diag_labels) <- ids_score
missingness_raw <- as.character(clin[ids_score, "missingness"])
if (anyNA(missingness_raw) || any(!nzchar(missingness_raw))) {
  stop("The prepared clinical missingness column contains missing values.")
}
missingness_levels <- c("CC", sort(setdiff(unique(missingness_raw), "CC")))
missingness_label <- factor(missingness_raw, levels = missingness_levels)
names(missingness_label) <- ids_score

# The shared method code names the regression outcome mmse_vec; for TCGA it holds the
# Karnofsky performance score
mmse_vec <- suppressWarnings(as.numeric(clin[ids_score, "karnofsky_performance_score"]))
names(mmse_vec) <- ids_score
regression_name <- "karnofsky_performance_score"

# Patients observed in every modality
common_ids <- Reduce(intersect, ids_by_modality)
common_ids <- intersect(common_ids, ids_score)

# CONSISTENCY CHECKS

if (!identical(rownames(X_stand), names(diag_labels))) {
  stop("TCGA diagnosis labels are not aligned to X_stand rows.")
}
if (!identical(rownames(X_stand), names(missingness_label))) {
  stop("TCGA missingness labels are not aligned to X_stand rows.")
}
if (!identical(colnames(X_stand), names(layer_vec))) {
  stop("TCGA layer_vec is not aligned to X_stand columns.")
}
if (anyDuplicated(rownames(X_stand))) {
  stop("TCGA X_stand contains duplicated sample IDs.")
}
if (anyDuplicated(colnames(X_stand))) {
  stop("TCGA X_stand contains duplicated internal feature names.")
}
if (nlevels(droplevels(diag_labels)) != 2L) {
  stop(
    "This TCGA version uses binary SIS feature selection and requires exactly ",
    "two disease_code groups. Found: ",
    paste(levels(droplevels(diag_labels)), collapse = ", ")
  )
}
if (length(common_ids) < 2L) {
  stop("Fewer than two TCGA patients are present in every included modality.")
}

# Summary of the prepared input
cat("============================================================\n")
cat("TCGA PREPARED INPUT SUMMARY\n")
cat("============================================================\n")
cat("Input directory:", TCGA_INPUT_DIR, "\n")
cat("Modalities:", paste(names(omics_list), collapse = ", "), "\n")
cat("Patients used:", nrow(X_stand), "\n")
cat("Features used:", ncol(X_stand), "\n")
cat("Patients present in every modality:", length(common_ids), "\n")
cat(
  "Duplicate cross-modality feature names disambiguated in memory:",
  sum(duplicated(original_feature_names)), "\n"
)
cat("Diagnosis counts:\n")
print(table(diag_labels, useNA = "ifany"))
cat("Missingness counts:\n")
print(table(missingness_label, useNA = "ifany"))
cat("============================================================\n")
