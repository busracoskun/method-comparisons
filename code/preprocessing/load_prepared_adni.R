# Preserved source expressions; see docs/reproducibility.md.
# Extracted from final.R:1441
path <- ADNI_INPUT_DIR

# Extracted from final.R:1442
imeta_raw <- read.csv(input_file(path, "meta.csv"))

# Extracted from final.R:1443
igene_raw <- read.csv(input_file(path, "rnaseq.csv"))

# Extracted from final.R:1444
clin <- read.csv(input_file(path, "clin.csv"))

# Extracted from final.R:1445
clin <- as.data.frame(clin)

# Extracted from final.R:1446
rownames(clin) <- clin$RID

# Extracted from final.R:1451
rownames(imeta_raw) <- imeta_raw$RID

# Extracted from final.R:1452
imeta_raw <- imeta_raw[, -1, drop = FALSE]

# Extracted from final.R:1453
rownames(igene_raw) <- igene_raw$RID

# Extracted from final.R:1454
igene_raw <- igene_raw[, -1, drop = FALSE]

# Extracted from final.R:1457
igene_raw <- igene_raw[, names(sort(matrixStats::colVars(as.matrix(igene_raw)), decreasing = TRUE))[1:500], drop = FALSE]

# Extracted from final.R:1459
ids_gene <- rownames(igene_raw)

# Extracted from final.R:1460
ids_meta <- rownames(imeta_raw)

# Extracted from final.R:1461
all_ids <- Reduce(union, list(ids_gene, ids_meta))

# Extracted from final.R:1464
igene <- igene_raw[match(all_ids, ids_gene), , drop = FALSE]

# Extracted from final.R:1465
imeta <- imeta_raw[match(all_ids, ids_meta), , drop = FALSE]

# Extracted from final.R:1466
rownames(igene) <- all_ids

# Extracted from final.R:1467
rownames(imeta) <- all_ids

# Extracted from final.R:1470
X_all <- cbind(igene, imeta)

# Extracted from final.R:1471
X_stand <- scale(X_all)

# Extracted from final.R:1474
layer_vec <- setNames(c(rep("gene", ncol(igene)), rep("meta", ncol(imeta))), colnames(X_stand))

# Extracted from final.R:1479
ids_with_clin <- intersect(all_ids, rownames(clin))

# Extracted from final.R:1480
diag_labels <- factor(clin[ids_with_clin, "diagnosis"], exclude = NULL)

# Extracted from final.R:1481
names(diag_labels) <- ids_with_clin

# Extracted from final.R:1482
mmse_vec <- setNames(clin[ids_with_clin, "MMSE"], ids_with_clin)

# Extracted from final.R:1484
missingness_label <- rep(NA_character_, length(all_ids))

# Extracted from final.R:1485
names(missingness_label) <- all_ids

# Extracted from final.R:1486
missingness_label[intersect(ids_gene, ids_meta)] <- "CC"

# Extracted from final.R:1487
missingness_label[setdiff(ids_gene, ids_meta)] <- "missing_meta"

# Extracted from final.R:1488
missingness_label[setdiff(ids_meta, ids_gene)] <- "missing_gene"

# Extracted from final.R:1489
missingness_label <- factor(missingness_label, levels = c("CC", "missing_meta", "missing_gene"))

# Extracted from final.R:1491
ids_score <- intersect(names(diag_labels), names(missingness_label))

# Extracted from final.R:1492
X_stand <- X_stand[ids_score, , drop = FALSE]

# Extracted from final.R:1493
missingness_label <- missingness_label[ids_score]

# Extracted from final.R:1494
diag_labels <- diag_labels[ids_score]

# Extracted from final.R:1495
mmse_vec <- mmse_vec[ids_score]

# Extracted from final.R:1496
common_ids <- intersect(intersect(ids_gene, ids_meta), ids_score)

