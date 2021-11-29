# Deng et al. Single-cell RNA-seq reveals dynamic, ranomd monoallelic gene expression in mammalian cells. Science (2014). 343(6167):193-6
# PMID: 24408435
# GEO: GSE45719


x <- read.table("Deng_counts.txt", header=T)
x <- x[,-grep("split", colnames(x))]
x <- x[,-grep("liver", colnames(x))]
x <- x[,-grep("fibro", colnames(x))]
x <- x[,-grep("smartseq2", colnames(x))]
x <- x[,-grep("C57", colnames(x))]
a <- colnames(x)
a <- strsplit(a, "_")

type <- unlist(lapply(a, function(b){b[[2]]}))
type <- sub("[1234]$", "", type)

require("SingleCellExperiment")
SingleCellExperiment(assays= list(counts=x), rowData=data.frame(cell_type1=type))
