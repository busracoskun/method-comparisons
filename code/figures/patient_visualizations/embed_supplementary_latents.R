# Fit (or load) a reference UMAP on latent scores and project further score sets onto it
# (usage: <reference_scores.csv> <model.uwot> <output.csv> [<query.csv> <output.csv>]...)

# Arguments and the reference scores
args <- commandArgs(trailingOnly = TRUE)
.libPaths(c(file.path(Sys.getenv("HOME"), "R", "4.5"), .libPaths()))
stopifnot(requireNamespace("uwot", quietly = TRUE))
Z <- as.matrix(read.csv(args[1], check.names = FALSE))
model_path <- args[2]
output <- args[3]

# Load the saved reference UMAP, or fit and save it
if (file.exists(model_path)) {
  reference <- uwot::load_uwot(model_path)
} else {
  set.seed(42)
  reference <- uwot::umap(
    Z,
    n_neighbors = min(15L, nrow(Z) - 1L), min_dist = .1, metric = "euclidean", n_components = 2L,
    n_epochs = 500L, scale = FALSE, n_threads = 1L, n_sgd_threads = 1L, ret_model = TRUE,
    verbose = FALSE
  )
  uwot::save_uwot(reference, model_path)
}

write.csv(
  data.frame(UMAP1 = reference$embedding[, 1], UMAP2 = reference$embedding[, 2]), output,
  row.names = FALSE
)

# Project each query score set onto the reference
if (length(args) > 3L) {
  for (i in seq.int(4L, length(args), 2L)) {
    query <- as.matrix(read.csv(args[i], check.names = FALSE))
    set.seed(42)
    xy <- uwot::umap_transform(
      query, reference, n_threads = 1L, n_sgd_threads = 1L, verbose = FALSE
    )
    write.csv(data.frame(UMAP1 = xy[, 1], UMAP2 = xy[, 2]), args[i + 1L], row.names = FALSE)
  }
}

# Record the UMAP settings
jsonlite::write_json(
  list(
    R = R.version.string, uwot = as.character(packageVersion("uwot")), seed = 42L,
    transform_seed = 42L, n_neighbors = min(15L, nrow(Z) - 1L), min_dist = .1, metric = "euclidean",
    n_epochs = 500L, n_threads = 1L, n_sgd_threads = 1L, reference_patients = nrow(Z),
    latent_dimensions = ncol(Z)
  ),
  paste0(model_path, ".settings.json"), auto_unbox = TRUE, pretty = TRUE
)
