# Plarania Atlas
# Plass et al. (2018) Cell Type Atlas and Lineage Tree of a Whole Complex Animalby Single-Cell Transcriptomics, Science
# doi : 10.1126/science.aaq1723

mat <- read.table("Planaria_dge.txt.gz")
ann <- read.table("Planaria_annotation.txt", sep=",")
pca_coords <- read.table("Planaria_PCA.txt")
pseudotime <- read.table("Planaria_pseudotime.txt")
cdat <- data.frame(cell_type1=ann, pseudo=pseudotime)
colnames(pca_coords) <- paste("PC", 1:50, sep="")
thing <- colnames(mat)
sample <- lapply(strsplit(thing, "_"), function(a){a[1]})
sample <- unlist(sample)
cdat$sample <- sample
cdat <- cbind(cdat, pca_coords)
colnames(cdat)[1:3] <- c("cell_type1", "pseudo", "sample");
require("SingleCellExperiment")
require("Matrix")
MAT <- Matrix(as.matrix(mat))
obj <- SingleCellExperiment(assays=list(counts=MAT), colData=cdat, rowData=data.frame(id=rownames(mat)))
saveRDS(obj, "Planaria.rds")
savehistory("Planaria_hist.R")
