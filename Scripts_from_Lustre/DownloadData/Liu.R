#Liu et al. (2016) Identification of key factors conquering developmental arrest of somatic cell cloned embryos by combining embryo biopsy and single-cell sequencing. Cell Discovery. 2 : 16010
# The Gene Expression Omnibus (GEO) accession numbers for the single-cell RNA-seq, ULI-NChIP-seq datasets are GSE70605, GSE70606, and GEO number for the single-cell RRBS is GSE70607.
# PMID: 27462457

data <- read.table("ExprMat.txt", stringsAsFactors=FALSE)
ann <- read.table("Liu_anno.txt", stringsAsFactors=FALSE)
