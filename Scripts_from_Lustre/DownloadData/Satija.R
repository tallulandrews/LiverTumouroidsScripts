# Zebrafish embryos (late blastula - 50% epiboly)
# Satija, R et al. (2015) Spatial reconstruction of single-cell gene expression data. Nature Biotechnology. 33(5):495-502
# GEO: GSE66688
# PMID: 25867923


x = read.table("GSE66688_zdata.expMatrix.txt.gz")


require("scater")
pd <- new("AnnotatedDataFrame", data=as.data.frame(A))
satija <- newSCESet(fpkmData=as.matrix(x), phenoData=pd, logExprsOffset=1)

saveRDS(satija, "satija.rds")

# ANXA2[A/B] and S100A10[A/B] hardly detected at all.
