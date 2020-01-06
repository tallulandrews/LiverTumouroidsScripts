prot <- read.delim("/lustre/scratch117/cellgen/team218/TA/OtherDownloadedData/normal_tissue.tsv", header=T)
rna <- read.delim("/lustre/scratch117/cellgen/team218/TA/OtherDownloadedData/rna_tissue.tsv", header=T)


prot$Level <- factor(prot$Level, levels=c("Not detected", "Low", "Medium", "High"))
prot$Reliability <- factor(prot$Reliability, levels=c("Uncertain", "Supported", "Approved"))

gene_row_sets <- split(seq(nrow(prot)), factor(prot$Gene))
tot_gene_expr <- sapply(gene_row_sets, function(a) sum(as.numeric(prot[a,"Level"])-1))
chol <- prot[prot$Tissue == "liver" & prot$Cell.type == "bile duct cells",]
chol$TotalExpr <- tot_gene_expr[match(chol$Gene, names(tot_gene_expr))]
chol$Proportion <- (as.numeric(chol$Level)-1)/chol$TotalExpr


hep <- prot[prot$Tissue == "liver" & prot$Cell.type == "hepatocytes",]
hep$TotalExpr <- tot_gene_expr[match(hep$Gene, names(tot_gene_expr))]
hep$Proportion <- (as.numeric(hep$Level)-1)/hep$TotalExpr

genes <- hep$Gene[hep$Gene %in% chol$Gene]
genes <- sort(genes)
hep <- hep[match(genes, hep$Gene),]
chol <- chol[match(genes, chol$Gene),]

all <- cbind(chol, hep)
saveRDS(all, file="HPA_Hep_Chol_Markers.rds")

good <- all[all[,6] != "Uncertain" & all[,14] != "Uncertain",]
good <- good[abs(as.numeric(good[,5])-as.numeric(good[,13])) >= 2,]
good$RelativeProp <- good[,8]/good[,16]
good <- good[ good$RelativeProp >= 2 | good$RelativeProp <= 2 ,]

saveRDS(good, file="Good_HPA_Hep_Chol_Markers.rds")

