# Mouse ESCs cell-cycle sorted
# Buettner et al. (2015). Computational analysis of cell-to-cell heterogeneity in single-cell RNA-sequencing data reveals hidden subpopulations of cells. Nature Biotechnology. 33:155-160
# PMID: 225599176
# ArrayExpress: E-MTAB-2805

g1_cells <- read.table("G1_singlecells_counts.txt", header=TRUE)
s_cells <- read.table("S_singlecells_counts.txt", header=TRUE)
g2m_cells <- read.table("G2M_singlecells_counts.txt", header=TRUE)

data <- cbind(g1_cells[,5:100], s_cells[,5:100], g2m_cells[,5:100])
rownames(data) <- g1_cells[,1]
anno <- rep(c("G1", "S", "G2M"), each=96)
gene_symbol <- g1_cells[,3]

ANN <- data.frame(Species = rep("Mus musculus", times=length(anno)), cell_type1 = anno, Source=rep("ESC", times=length(anno)))
rownames(ANN) <- colnames(data);

require("scater")
pd <- new("AnnotatedDataFrame", data=ANN)
buet <- SingleCellExperiment(assays=list(counts=as.matrix(data)), colData=ANN)
rowData(buet)$feature_symbol <- gene_symbol

saveRDS(buet, "buettner.rds")

