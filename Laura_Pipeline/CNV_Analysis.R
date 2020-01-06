locations <- read.delim("/lustre/scratch117/cellgen/team218/TA/genomebuilding/Homo_sapiens.GRCh38.79.genepos.txt", header=F)
locations <- locations[locations[,2] %in% 1:22,]
locations <- locations[order(as.numeric(locations[,2]), locations[,3], decreasing=FALSE),]

require("SingleCellExperiment")
require("scater")

SCE <- readRDS("Global_SCE.rds")
mat <- assays(SCE)[["lognorm"]]
mat <- mat[rowMeans(mat) > 1, ]

loc <- locations[locations[,1] %in% rownames(mat), ]
mat <- mat[match(loc[,1],rownames(mat)),]
toplot <- t(apply(mat, 1, scale))

c25 <- c("dodgerblue2","#E31A1C", # red
                "green4",
                "#6A3D9A", # purple
                "#FF7F00", # orange
                "black","gold1",
                "skyblue2","#FB9A99", # lt pink
                "palegreen2",
                "#CAB2D6", # lt purple
                "#FDBF6F", # lt orange
                "gray70", "khaki2",
                "maroon","orchid1","deeppink1","blue1","steelblue4",
                "darkturquoise","green1","yellow4","yellow3",
                "darkorange4","brown")

rowcols <- c25[as.numeric(loc[,2])];
Donor.col <- c("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", "#fdbf6f", "#ff7f00")

colcols <- Donor.col[SCE$Donor]
require("gplots")
#heatmap.2(toplot, Rowv=FALSE, Colv=FALSE, trace="none", RowSideColors=rowcols, ColSideColors=colcols, dendrogram="none", distfun=function(x) {as.dist(matrix(1, ncol=nrow(x), nrow=nrow(x)))})
png("InitCNV.png", width=8, height=6, units="in", res=300)
image(toplot, col=c("#2171b5", "#6baed6", "white", "#fb6a4a", "#cb181d"), breaks=c(-10000, -2,-1,1, 2, 10000), axes=FALSE)
# Row colors
res <- table(factor(SCE$Donor))
res <- res/sum(res)
res <- res[match(unique(SCE$Donor), names(res))]
for (i in 2:length(res)) {
	res[i] <- res[i-1]+res[i]
}
arrows(rep(-0.01,length(Donor.col)),
       c(0,res[-length(res)]), 
       rep(-0.01,length(Donor.col)), 
       res, 
       col=Donor.col, lwd=3, xpd=T, len=0)
abline(h=res, lwd=0.75, lty=2, col="black")
# Col colors
res <- table(factor(loc[,2]))
res <- res/sum(res)
res <- res[order(as.numeric(names(res)))]
chr.col <- c25[1:length(res)]
for (i in 2:length(res)) {
	res[i] <- res[i-1]+res[i]
}
arrows(c(0,res[-length(res)]), 
       rep(-0.01, max(as.numeric(loc[,2]))), 
       res,  
       rep(-0.01, max(as.numeric(loc[,2]))), 
       col=chr.col, lwd=3, xpd=T, len=0)
abline(v=res, lwd=0.75, lty=2, col="black")

title(xlab="Chromosome", ylab="Donor")
dev.off()


#### inferCNV (Broad) ####
require("scater")
SCE <- readRDS("Global_SCE.rds")
mat <- assays(SCE)[["lognorm"]]
ann <- cbind(as.character(colnames(SCE)), as.character(SCE$Donor))
#ann <- cbind(as.character(colnames(SCE)), paste(as.character(SCE$CC_state), as.character(SCE$Donor), sep=""))
tumour <- !(SCE$Donor %in% c("D3EM", "D3DM", "D9DM", "D9EM"))
ann[tumour,2] <- paste("malignant", ann[tumour,2], sep="_")
gene_order <- read.table("/lustre/scratch117/cellgen/team218/TA/genomebuilding/Homo_sapiens.GRCh38.79.genepos.txt")
gene_order[,2] <- paste("chr", as.character(gene_order[,2]), sep="")
gene_order <- gene_order[order(gene_order[,2], gene_order[,3], gene_order[,4]), ]
gene_order <- gene_order[gene_order[,1] %in% rownames(mat),]
mat <- mat[match( as.character(gene_order[,1]), rownames(mat) ), ]
write.table(mat, file="inferCNV_input_counts.tsv", sep="\t", col.names=TRUE, row.names=TRUE)
write.table(ann, file="inferCNV_input_ann.txt", sep="\t", col.names=FALSE, row.names=FALSE)
write.table(ann[!tumour,1], file="inferCNV_refcells.txt", sep="\t", col.names=FALSE, row.names=FALSE)
write.table(gene_order, file="inferCNV_input_geneorder.txt", sep="\t", col.names=FALSE, row.names=FALSE)


require(infercnv)
infercnv_obj <- CreateInfercnvObject(raw_counts_matrix="inferCNV_input_counts.tsv", annotations_file="inferCNV_input_ann.txt", delim="\t", gene_order_file="inferCNV_input_geneorder.txt", ref_group_names=c("D3DM", "D3EM", "D9DM", "D9EM"))
infercnv_obj = infercnv::run(infercnv_obj, cutoff=1, out_dir="inferCNVRes/", cluster_by_groups=T, include.spike=T)








# The below fails!
require("infercnv")
out <- infer_cnv(mat, gene_order, cutoff=1, reference_obs=ann[!tumour,1], transform_data=FALSE, window_length=11, max_centered_threshold=1, noise_threshold=0.1, num_ref_groups=1, out_path="/lustre/scratch117/cellgen/team218/TA/LiverOrganoids/Analysis_Pipeline_Output/inferCNVRes")


data = mat
cutoff = 1
reference_obs = ann[!tumour,1]
transform_data = FALSE
window_length = 11
max_centered_threshold = 1
noise_threshold = 0.1
num_ref_groups = 1
out_path = "/lustre/scratch117/cellgen/team218/TA/LiverOrganoids/Analysis_Pipeline_Output/inferCNVRes"
k_obs_groups = 1
plot_steps = FALSE
contig_tail = (window_length - 1) / 2
method_bound_vis = NA
lower_bound_vis = NA
upper_bound_vis = NA
ref_subtract_method = "by_mean"
hclust_method = 'complete'



