require("fastICA")
require("Rtsne")
require("scater")
require("SingleCellExperiment")
require("destiny")

set.seed(1328298)

# Feature selection


# Plot Clusters in reduced dimensions
ica <- fastICA(mat, n.comp=k)

pca <- prcomp(mat)

tsne <- Rtsne(t(mat), preplexity=K)

monocle_dims <-


# Colour cells by expression level
