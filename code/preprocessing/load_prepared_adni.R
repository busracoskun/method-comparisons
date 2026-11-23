# Load and align the prepared ADNI gene-expression, metabolomics, and clinical inputs
# Sourced by load_dataset(); defines X_stand, layer_vec, diag_labels, mmse_vec,
# missingness_label, and common_ids

# Load ADNI data
path <- ADNI_INPUT_DIR

imeta_raw <- read.csv(input_file(path, "meta.csv"))
igene_raw <- read.csv(input_file(path, "rnaseq.csv"))
clin <- read.csv(input_file(path, "clin.csv"))
clin <- as.data.frame(clin)

rownames(clin) <- clin$RID
rownames(imeta_raw) <- imeta_raw$RID
rownames(igene_raw) <- igene_raw$RID

imeta_raw <- imeta_raw[, -1, drop = FALSE]
igene_raw <- igene_raw[, -1, drop = FALSE]

# Retain the 500 gene-expression features with the highest variance
gene_vars <- matrixStats::colVars(as.matrix(igene_raw))
igene_raw <- igene_raw[, names(sort(gene_vars, decreasing = TRUE))[1:500], drop = FALSE]

# Align patients across modalities, retaining incomplete observations
ids_gene <- rownames(igene_raw)
ids_meta <- rownames(imeta_raw)
all_ids <- Reduce(union, list(ids_gene, ids_meta))

igene <- igene_raw[match(all_ids, ids_gene), , drop = FALSE]
imeta <- imeta_raw[match(all_ids, ids_meta), , drop = FALSE]

rownames(igene) <- all_ids
rownames(imeta) <- all_ids

# Combine and standardize molecular features
X_all <- cbind(igene, imeta)
X_stand <- scale(X_all)

layer_vec <- setNames(c(rep("gene", ncol(igene)), rep("meta", ncol(imeta))), colnames(X_stand))

# Extract clinical outcomes for patients with available metadata
ids_with_clin <- intersect(all_ids, rownames(clin))

diag_labels <- factor(clin[ids_with_clin, "diagnosis"], exclude = NULL)
names(diag_labels) <- ids_with_clin
mmse_vec <- setNames(clin[ids_with_clin, "MMSE"], ids_with_clin)

# Label patients according to observed modality availability
missingness_label <- rep(NA_character_, length(all_ids))
names(missingness_label) <- all_ids

missingness_label[intersect(ids_gene, ids_meta)] <- "CC"
missingness_label[setdiff(ids_gene, ids_meta)] <- "missing_meta"
missingness_label[setdiff(ids_meta, ids_gene)] <- "missing_gene"
missingness_label <- factor(missingness_label, levels = c("CC", "missing_meta", "missing_gene"))

# Restrict analysis matrices and outcomes to patients with clinical data
ids_score <- intersect(names(diag_labels), names(missingness_label))

X_stand <- X_stand[ids_score, , drop = FALSE]
missingness_label <- missingness_label[ids_score]
diag_labels <- diag_labels[ids_score]
mmse_vec <- mmse_vec[ids_score]

common_ids <- intersect(intersect(ids_gene, ids_meta), ids_score)
