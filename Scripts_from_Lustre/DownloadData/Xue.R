# Human & Mouse pre-implantation embryos
# Xue et al. (2013) Genetic programs in human and mouse early embryos revealed by single-cell RNA sequencing. Nature. 500(7464):593-7.
#GEO: GSE44183
#PMID: 23892778


# Human
h <- read.delim("GSE44183_human_expression_mat.txt.gz", header=TRUE)

rownames(h) <- h[,1]
h <- h[,-1]
type <- colnames(h)
type[grep("oocyte", type)] = "oocyte"
type[grep("zygote", type)] = "zygote"
type[grep("pronuc", type)] = "zygote"
type[grep("2.cell", type)] = "2cell"
type[grep("4.cell", type)] = "4cell"
type[grep("8.cell", type)] = "8cell"
type[grep("morula", type)] = "16cell"

require("scater")
P <- data.frame(cell_type1=type); rownames(P) <- colnames(h);
pd <- new("AnnotatedDataFrame", data=P)
XH <- newSCESet(fpkmData=as.matrix(h), phenoData=pd, logExprsOffset=1)
fData(XH)$feature_symbol <- rownames(exprs(XH));
saveRDS(XH, "Xue_human.rds")


# Mouse
m <- read.delim("GSE44183_mouse_expression_mat.txt.gz", header=TRUE)
rownames(m) <- m[,1]
m <- m[,-c(1, ncol(m))]
type <- colnames(m);
type[grep("Oocyte", type)] <- "oocyte";
type[grep("2cell", type)] <- "2cell";
type[grep("4cell", type)] <- "4cell";
type[grep("8cell", type)] <- "8cell";
type[grep("Morula", type)] <- "16cell";
type[grep("Pronucl", type)] <- "zygote";

require("scater")
P <- data.frame(cell_type1=type); rownames(P) <- colnames(m);
pd <- new("AnnotatedDataFrame", data=P)
XM <- newSCESet(fpkmData=as.matrix(m), phenoData=pd, logExprsOffset=1)
fData(XM)$feature_symbol <- rownames(exprs(XM));
saveRDS(XM, "Xue_mouse.rds")






