# Explicit optional setup. Never executed by the analyses or smoke tests
cran <- c(
  "missForest", "SIS", "dplyr", "tibble", "matrixStats", "clusterCrit", "mclust",
  "r.jive", "RGCCA", "IntNMF", "GFA", "OmicsPLS", "ade4", "SNFtool", "FactoMineR", "T4transport",
  "gplots", "denoiseR", "Matrix", "remotes", "uwot", "jsonlite"
)
install.packages(setdiff(cran, rownames(installed.packages())), repos = "https://cloud.r-project.org")
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
bioc <- c("iClusterPlus", "fabia", "omicade4", "mixOmics")
BiocManager::install(setdiff(bioc, rownames(installed.packages())), ask = FALSE, update = FALSE)
repos <- c(
  GPCCA = "Kaversoniano/GPCCA", BSFP = "sarahsamorodnitsky/BSFP",
  tensorICA = "jinghan1018/tensorICA", SLIDE = "irinagain/SLIDE", NEMO = "Shamir-Lab/NEMO/NEMO"
)
for (pkg in names(repos)) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    remotes::install_github(repos[[pkg]], dependencies = TRUE, upgrade = "never")
  }
}
# GIPCA requires the historical source snapshot; see docs/reproducibility.md
# These upstream locations are not version locks. Capture SHAs before a public release
