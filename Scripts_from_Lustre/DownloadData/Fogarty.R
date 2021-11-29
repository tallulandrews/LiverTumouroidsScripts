# Human embryos + Oct4 KO
# Fogarty et al. (2017) Genome editing reveals a role for OCT4 in human embryogenesis. Nature. 550(7674):67-73.
# GEO: GSE100118
# PMID: 28953884



x <- read.delim("GSE100118_scRNA_pou5f1crispr_rpkm_170603.csv.gz", header=TRUE, sep=",")

tmp <- x[1:2,]
colnames(x) <- as.character(unlist(as.vector(x[2,])))
x <- x[-c(1,2),]

F <- data.frame(ensembl_gene_id=x[,1], feature_symbol=x[,2])
rownames(F) <- F[,1]
x <- x[,-c(1,2)]
rownames(x) <- rownames(F)
tmp <- tmp[,-c(1,2)]

tmp <- t(tmp)

P <- data.frame(A=1:length(tmp[,1]), B= 1:length(tmp[,1]));
rownames(P) <- colnames(x);
colnames(P) <- c("cell_type1", "EmbryoID")
P <- t(P)
P[,grep("^2", tmp[,1])] <- c("blast", 2)
P[,grep("^7", tmp[,1])] <- c("blast", 7)
P[,grep("^8", tmp[,1])] <- c("blast", 8)
P[,grep("^5", tmp[,1])] <- c("blast", 5)
P[,grep("TE", tmp[,1])] <- c("TE", NA)
P[,grep("C16", tmp[,1])] <- c("KOblast", 16)
P[,grep("C8", tmp[,1])] <- c("KOblast", 16)
P[,grep("C12", tmp[,1])] <- c("KOblast", 16)
P[,grep("C9", tmp[,1])] <- c("KOblast", 16)
P[,grep("C24", tmp[,1])] <- c("KOblast", 16)

y <- apply(x, 1, as.numeric)
y <- y[,-22554]
F <- F[-22554,]
rownames(y) <- colnames(x);


require("scater")
pd <- new("AnnotatedDataFrame", data=as.data.frame(t(P)))
fd <- new("AnnotatedDataFrame", data=F)
fogarty <- newSCESet(fpkmData=as.matrix(t(y)), phenoData=pd, featureData=fd, logExprsOffset=1)

saveRDS(fogarty, "fogarty.rds")

