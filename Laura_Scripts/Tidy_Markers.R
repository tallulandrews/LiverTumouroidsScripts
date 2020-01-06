
#surface_terms = c("GO:0009986","GO:0004872","GO:0005102")
surface_genes <- read.table("/nfs/users/nfs_t/ta6/Collaborations/LiverOrganoids/Hsap_Surface_Genes.txt", header=T)

x <- readRDS("Laura_SC3_k6_Markers.rds")

CCA1 <- x$CCA1
CCA1$signif <- CCA1[,3] < 0.05/length(CCA1[,1])/6
CCA1$surface = CCA1$Gene %in% surface_genes[,1]
colnames(CCA1) <- c("AUC", "Cluster", "p.value", "Gene", "Symbol", "Significant", "GO Surface")
my_order = order(-CCA1$AUC)
write.table(CCA1[my_order,], file="CCA1_full_markers.txt", col.names=T, row.names=F, sep="\t")



HCC6 <- x$HCC6
HCC6$signif <- HCC6[,3] < 0.05/length(HCC6[,1])/6
HCC6$surface = HCC6$Gene %in% surface_genes[,1]
colnames(HCC6) <- c("AUC", "Cluster", "p.value", "Gene", "Symbol", "Significant", "GO Surface")
my_order = order(-HCC6$AUC)
write.table(HCC6[my_order,], file="HCC6_full_markers.txt", col.names=T, row.names=F, sep="\t")


HCC10 <- x$HCC10
HCC10$signif <- HCC10[,3] < 0.05/length(HCC10[,1])/6
HCC10$surface = HCC10$Gene %in% surface_genes[,1]
colnames(HCC10) <- c("AUC", "Cluster", "p.value", "Gene", "Symbol", "Significant", "GO Surface")
my_order = order(-HCC10$AUC)
write.table(HCC10[my_order,], file="HCC10_full_markers.txt", col.names=T, row.names=F, sep="\t")


