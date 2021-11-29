# Human embryos and ESCs
# Yan et al. (2013) Single-cell RNA-Seq profiling of human preimplantation embryos and embryonic stem cells. Nature Structural Molecular Biology. 20(9):1131-9
# GEO: GSE36552
# PMID: 23934149



a = read.table("Anno.txt")
x = read.table("ExprMat.txt")

ann <- t(a)

ann[grep("stem cell", ann[,3]),3] <- "hESC"
ann[grep("2-cell", ann[,3]),3] <- "2cell"
ann[grep("Oocyte", ann[,3]),3] <- "oocyte"
ann[grep("ygote", ann[,3]),3] <- "zygote"
ann[grep("4-cell", ann[,1]),3] <- "4cell"
ann[grep("8-cell", ann[,1]),3] <- "8cell"
ann[grep("orulae", ann[,1]),3] <- "16cell"
ann[grep("blastocyst", ann[,1]),3] <- "blast"

A <- cbind(ann[,3], ann[,5])
rownames(A) <- ann[,6]
colnames(A) <- c("cell_type1", "source")



require("scater")
pd <- new("AnnotatedDataFrame", data=as.data.frame(A))
tang <- newSCESet(fpkmData=as.matrix(x), phenoData=pd, logExprsOffset=1)

saveRDS(tang, "tang.rds")


