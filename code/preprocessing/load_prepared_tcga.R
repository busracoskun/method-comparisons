# Preserved source expressions; see docs/reproducibility.md.
# Extracted from final_tcga.R:91
TCGA_MRNA_FILE <- file.path(TCGA_INPUT_DIR, "irsem_fs.csv")

# Extracted from final_tcga.R:92
TCGA_METHYLATION_FILE <- file.path(TCGA_INPUT_DIR, "imeth_fs.csv")

# Extracted from final_tcga.R:93
TCGA_RPPA_FILE <- file.path(TCGA_INPUT_DIR, "irppa_fs.csv")

# Extracted from final_tcga.R:94
TCGA_MIRNA_FILE <- file.path(TCGA_INPUT_DIR, "imirna_fs.csv")

# Extracted from final_tcga.R:95
TCGA_CLINICAL_FILE <- file.path(TCGA_INPUT_DIR, "clin_common.csv")

# Extracted from final_tcga.R:1403
required_tcga_files <- c(TCGA_MRNA_FILE, TCGA_METHYLATION_FILE, TCGA_RPPA_FILE, TCGA_CLINICAL_FILE)

# Extracted from final_tcga.R:1410
missing_tcga_files <- required_tcga_files[!file.exists(required_tcga_files)]

# Extracted from final_tcga.R:1411
if (length(missing_tcga_files)) {
  stop(
    "Required prepared TCGA file(s) were not found:\n",
    paste(missing_tcga_files, collapse = "\n")
  )
}

# Extracted from final_tcga.R:1418
read_prepared_omics <- function(file_path, modality) {
  cat("Reading", modality, "from:", file_path, "\n")
  
  dat <- read.csv(
    file_path,
    header = TRUE,
    check.names = FALSE,
    stringsAsFactors = FALSE,
    na.strings = c("", "NA", "N/A", "NaN", "NULL"),
    quote = "\"",
    comment.char = ""
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
  
  nonnumeric <- names(feature_data)[
    !vapply(feature_data, is.numeric, logical(1))
  ]
  if (length(nonnumeric)) {
    stop(
      modality,
      " contains nonnumeric feature columns: ",
      paste(nonnumeric, collapse = ", ")
    )
  }
  
  M <- as.matrix(feature_data)
  storage.mode(M) <- "double"
  rownames(M) <- sample_ids
  
  if (any(is.infinite(M))) {
    stop(modality, " contains positive or negative infinity.")
  }
  
  cat(
    "  Samples:", nrow(M),
    "| Features:", ncol(M),
    "| Missing cells:", sum(is.na(M)),
    "\n"
  )
  
  M
}

# Extracted from final_tcga.R:1477
read_prepared_clinical <- function(file_path) {
  cat("Reading clinical data from:", file_path, "\n")
  
  clin <- read.csv(
    file_path,
    header = TRUE,
    check.names = FALSE,
    stringsAsFactors = FALSE,
    na.strings = c("", "NA", "N/A", "NaN", "NULL"),
    quote = "\"",
    comment.char = ""
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

# Extracted from final_tcga.R:1506
align_prepared_matrix <- function(M, ids) {
  out <- matrix(NA_real_, nrow = length(ids), ncol = ncol(M), dimnames = list(ids, colnames(M)))
  
  shared <- intersect(ids, rownames(M))
  out[shared, ] <- M[shared, , drop = FALSE]
  out
}

# Extracted from final_tcga.R:1519
omics_list <- list(
  mrna = read_prepared_omics(TCGA_MRNA_FILE, "mRNA"),
  methylation = read_prepared_omics(
    TCGA_METHYLATION_FILE,
    "methylation"
  ),
  rppa = read_prepared_omics(TCGA_RPPA_FILE, "RPPA")
)

# Extracted from final_tcga.R:1530
if (file.exists(TCGA_MIRNA_FILE)) {
  omics_list$mirna <- read_prepared_omics(TCGA_MIRNA_FILE, "miRNA")
} else {
  cat(
    "Optional prepared miRNA file not found; continuing without miRNA:\n",
    TCGA_MIRNA_FILE,
    "\n"
  )
}

# Extracted from final_tcga.R:1540
clin <- read_prepared_clinical(TCGA_CLINICAL_FILE)

# Extracted from final_tcga.R:1542
required_clinical_columns <- c(
  "disease_code",
  "karnofsky_performance_score",
  "missingness"
)

# Extracted from final_tcga.R:1548
missing_clinical_columns <- setdiff(required_clinical_columns, colnames(clin))

# Extracted from final_tcga.R:1553
if (length(missing_clinical_columns)) {
  stop(
    "clin_common.csv is missing required columns: ",
    paste(missing_clinical_columns, collapse = ", ")
  )
}

# Extracted from final_tcga.R:1560
ids_by_modality <- lapply(omics_list, rownames)

# Extracted from final_tcga.R:1561
omics_union_ids <- Reduce(union, ids_by_modality)

# Extracted from final_tcga.R:1565
ids_score <- intersect(omics_union_ids, rownames(clin))

# Extracted from final_tcga.R:1567
diagnosis_raw <- as.character(clin[ids_score, "disease_code"])

# Extracted from final_tcga.R:1568
keep_diagnosis <- !is.na(diagnosis_raw) & nzchar(diagnosis_raw)

# Extracted from final_tcga.R:1570
if (any(!keep_diagnosis)) {
  cat(
    "Dropping",
    sum(!keep_diagnosis),
    "patients without disease_code.\n"
  )
}

# Extracted from final_tcga.R:1578
ids_score <- ids_score[keep_diagnosis]

# Extracted from final_tcga.R:1579
diagnosis_raw <- diagnosis_raw[keep_diagnosis]

# Extracted from final_tcga.R:1581
if (length(ids_score) < 2L) {
  stop("Fewer than two TCGA patients remain after clinical matching.")
}

# Extracted from final_tcga.R:1585
omics_aligned <- lapply(omics_list, align_prepared_matrix, ids = ids_score)

# Extracted from final_tcga.R:1595
original_feature_names <- unlist(lapply(omics_aligned, colnames), use.names = FALSE)

# Extracted from final_tcga.R:1600
internal_feature_names <- make.unique(
  original_feature_names,
  sep = "."
)

# Extracted from final_tcga.R:1605
feature_modalities <- rep(names(omics_aligned), times = vapply(omics_aligned, ncol, integer(1)))

# Extracted from final_tcga.R:1610
feature_name_map <- data.frame(
  internal_name = internal_feature_names,
  original_name = original_feature_names,
  modality = feature_modalities,
  stringsAsFactors = FALSE
)

# Extracted from final_tcga.R:1617
cursor <- 1L

# Extracted from final_tcga.R:1618
for (modality in names(omics_aligned)) {
  width <- ncol(omics_aligned[[modality]])
  idx <- cursor:(cursor + width - 1L)
  colnames(omics_aligned[[modality]]) <- internal_feature_names[idx]
  cursor <- cursor + width
}

# Extracted from final_tcga.R:1625
X_all <- do.call(cbind, omics_aligned)

# Extracted from final_tcga.R:1626
storage.mode(X_all) <- "double"

# Extracted from final_tcga.R:1631
feature_means <- colMeans(X_all, na.rm = TRUE)

# Extracted from final_tcga.R:1632
feature_sds <- vapply(seq_len(ncol(X_all)), function(j) stats::sd(X_all[, j], na.rm = TRUE), numeric(1))

# Extracted from final_tcga.R:1633
names(feature_means) <- names(feature_sds) <- colnames(X_all)

# Extracted from final_tcga.R:1635
if (any(!is.finite(feature_sds) | feature_sds <= 0)) {
  stop("TCGA contains zero-variance or invalid features before standardization.")
}

# Extracted from final_tcga.R:1639
X_stand <- sweep(X_all, 2, feature_means, FUN = "-")

# Extracted from final_tcga.R:1640
X_stand <- sweep(X_stand, 2, feature_sds, FUN = "/")

# Extracted from final_tcga.R:1641
storage.mode(X_stand) <- "double"

# Extracted from final_tcga.R:1643
layer_vec <- setNames(feature_modalities, internal_feature_names)

# Extracted from final_tcga.R:1648
diag_labels <- factor(diagnosis_raw)

# Extracted from final_tcga.R:1649
names(diag_labels) <- ids_score

# Extracted from final_tcga.R:1651
missingness_raw <- as.character(
  clin[ids_score, "missingness"]
)

# Extracted from final_tcga.R:1655
if (anyNA(missingness_raw) || any(!nzchar(missingness_raw))) {
  stop("The prepared clinical missingness column contains missing values.")
}

# Extracted from final_tcga.R:1659
missingness_levels <- c(
  "CC",
  sort(setdiff(unique(missingness_raw), "CC"))
)

# Extracted from final_tcga.R:1664
missingness_label <- factor(missingness_raw, levels = missingness_levels)

# Extracted from final_tcga.R:1668
names(missingness_label) <- ids_score

# Extracted from final_tcga.R:1670
mmse_vec <- suppressWarnings(
  as.numeric(
    clin[ids_score, "karnofsky_performance_score"]
  )
)

# Extracted from final_tcga.R:1675
names(mmse_vec) <- ids_score

# Extracted from final_tcga.R:1676
regression_name <- "karnofsky_performance_score"

# Extracted from final_tcga.R:1681
common_ids <- Reduce(intersect, ids_by_modality)

# Extracted from final_tcga.R:1682
common_ids <- intersect(common_ids, ids_score)

# Extracted from final_tcga.R:1684
if (!identical(rownames(X_stand), names(diag_labels))) {
  stop("TCGA diagnosis labels are not aligned to X_stand rows.")
}

# Extracted from final_tcga.R:1687
if (!identical(rownames(X_stand), names(missingness_label))) {
  stop("TCGA missingness labels are not aligned to X_stand rows.")
}

# Extracted from final_tcga.R:1690
if (!identical(colnames(X_stand), names(layer_vec))) {
  stop("TCGA layer_vec is not aligned to X_stand columns.")
}

# Extracted from final_tcga.R:1693
if (anyDuplicated(rownames(X_stand))) {
  stop("TCGA X_stand contains duplicated sample IDs.")
}

# Extracted from final_tcga.R:1696
if (anyDuplicated(colnames(X_stand))) {
  stop("TCGA X_stand contains duplicated internal feature names.")
}

# Extracted from final_tcga.R:1699
if (nlevels(droplevels(diag_labels)) != 2L) {
  stop(
    "This TCGA version uses binary SIS feature selection and requires exactly ",
    "two disease_code groups. Found: ",
    paste(levels(droplevels(diag_labels)), collapse = ", ")
  )
}

# Extracted from final_tcga.R:1706
if (length(common_ids) < 2L) {
  stop("Fewer than two TCGA patients are present in every included modality.")
}

# Extracted from final_tcga.R:1710
cat("============================================================\n")

# Extracted from final_tcga.R:1711
cat("TCGA PREPARED INPUT SUMMARY\n")

# Extracted from final_tcga.R:1712
cat("============================================================\n")

# Extracted from final_tcga.R:1713
cat("Input directory:", TCGA_INPUT_DIR, "\n")

# Extracted from final_tcga.R:1714
cat("Modalities:", paste(names(omics_list), collapse = ", "), "\n")

# Extracted from final_tcga.R:1715
cat("Patients used:", nrow(X_stand), "\n")

# Extracted from final_tcga.R:1716
cat("Features used:", ncol(X_stand), "\n")

# Extracted from final_tcga.R:1717
cat("Patients present in every modality:", length(common_ids), "\n")

# Extracted from final_tcga.R:1718
cat(
  "Duplicate cross-modality feature names disambiguated in memory:",
  sum(duplicated(original_feature_names)),
  "\n"
)

# Extracted from final_tcga.R:1723
cat("Diagnosis counts:\n")

# Extracted from final_tcga.R:1724
print(table(diag_labels, useNA = "ifany"))

# Extracted from final_tcga.R:1725
cat("Missingness counts:\n")

# Extracted from final_tcga.R:1726
print(table(missingness_label, useNA = "ifany"))

# Extracted from final_tcga.R:1727
cat("============================================================\n")

