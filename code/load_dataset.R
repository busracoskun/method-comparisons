# R session setup and dataset loader shared by every R analysis script
# Scripts are launched by run_benchmark.py, which sets CLEAN_ROOT and the input/output paths

CLEAN_ROOT <- Sys.getenv("CLEAN_ROOT")
if (!nzchar(CLEAN_ROOT)) stop("Use python run_benchmark.py so paths are configured.")
options(stringsAsFactors = FALSE)
.libPaths(c(Sys.getenv("R_LIBS_USER"), file.path(Sys.getenv("HOME"), "R", "4.5"), .libPaths()))

# Build an environment holding one prepared dataset and the method code. The sourcing order is
# (helpers, methods, inputs, then observed-data strategies) 
load_dataset <- function(dataset, observed = FALSE) {
  stopifnot(dataset %in% c("adni", "tcga"))
  suppressPackageStartupMessages(library(dplyr))

  # Dataset environment with input and output locations
  e <- new.env(parent = .GlobalEnv)
  e$CLEAN_ROOT <- CLEAN_ROOT
  e$PROJECT_DIR <- Sys.getenv("BENCH_WORK_DIR")
  e$ADNI_INPUT_DIR <- Sys.getenv("ADNI_INPUT_DIR")
  e$TCGA_INPUT_DIR <- Sys.getenv("TCGA_INPUT_DIR")

  # Input file, or its .gz version
  e$input_file <- function(folder, name) {
    path <- file.path(folder, name)
    if (!file.exists(path) && file.exists(paste0(path, ".gz"))) paste0(path, ".gz") else path
  }

  # Cached variable file, preferring CACHE_INPUT_DIR when it holds a copy
  e$cache_file <- function(name) {
    folder <- Sys.getenv("CACHE_INPUT_DIR")
    cached <- file.path(folder, name)
    if (nzchar(folder) && file.exists(cached)) cached else file.path(e$VARIABLE_DIR, name)
  }

  e$VARIABLE_DIR <- file.path(e$PROJECT_DIR, "variables")
  e$CLEAN_MODEL_DIR <- file.path(e$PROJECT_DIR, "fits")
  dir.create(e$VARIABLE_DIR, recursive = TRUE, showWarnings = FALSE)

  # Method code
  if ("bidifac" %in% strsplit(Sys.getenv("R_METHODS"), ",", fixed = TRUE)[[1]]) {
    # Prepared inputs and, for observed-data runs, the strategy-specific matrices
    source(
      Sys.getenv("BIDIFAC_FILE", file.path(CLEAN_ROOT, "external_sources", "bidifac", "BIDIFAC.R")),
      local = e
    )
  }
  source(file.path(CLEAN_ROOT, "code", "methods", "shared.R"), local = e)
  source(file.path(CLEAN_ROOT, "code", "methods", paste0(dataset, ".R")), local = e)

  # Some wrappers call package functions unqualified, as in the original library-loaded session
  packages <- c(icluster = "iClusterPlus", fabia = "fabia", intnmf = "IntNMF", nemo = "SNFtool")
  selected <- strsplit(Sys.getenv("R_METHODS"), ",", fixed = TRUE)[[1]]
  for (p in unique(unname(packages[intersect(names(packages), selected)]))) {
    suppressPackageStartupMessages(library(p, character.only = TRUE))
  }

  source(
    file.path(CLEAN_ROOT, "code", "preprocessing", paste0("load_prepared_", dataset, ".R")), local = e
  )
  if (observed) {
    prepare_file <- file.path(CLEAN_ROOT, "code", "observed_data", paste0("prepare_", dataset, ".R"))
    source(prepare_file, local = e)
  }
  e
}
