Camp <- readRDS("/lustre/scratch117/cellgen/team218/TA/scRNASeqDatasets/camp_human.rds")
require("scater")
require("RColorBrewer")

type_col_palette = colorRampPalette(c("grey50","red","orange","goldenrod","forestgreen","cornflowerblue","navy","purple"))
#type_col = brewer.pal(length(levels(pData(Camp)$cell_type1)), "Set2")
type_col = type_col_palette(length(levels(pData(Camp)$cell_type1)))
source_pch = c(18, 7, 25, 1, 15)

map = read.table("/lustre/scratch117/cellgen/team218/TA/Mirrored_Annotations/Hsap_Gene_Name_Mapping_Ensembl80.out", header=T)
ensg2symbol <- function(x) {
        new = as.character(map[match(x, map[,1]),2])
        new[is.na(new)] = as.character(x[is.na(new)])
        new[duplicated(new)] = x[duplicated(new)]
        return(new)
}

biotype <- read.table("/lustre/scratch117/cellgen/team218/TA/LiverOrganoids/Human_biotype.txt")
biotype[,1] <- ensg2symbol(biotype[,1])
pc <- biotype[biotype[,2] == "protein_coding",1]
mt <- biotype[grep("Mt_",biotype[,2]),1]
rRNA <- biotype[biotype[,2] == "rRNA",1]
sRNAs <- biotype[biotype[,2] == "sRNA" | biotype[,2] == "snRNA" | biotype[,2] == "snoRNA" | biotype[,2] == "scaRNA",1]

my_QC <- function(dat, type, plate) {
	par(mfrow=c(2,2))
	par(mar=c(4,4,1,1))
	plot(colSums(dat), colSums(dat > 0), col=type_col[type], pch=source_pch[plate], xlab="total counts", ylab="total genes")

	min_detect = 3500
	min_counts = 15000
	abline(h=min_detect)
	abline(v=min_counts) 
	filter_detect <- colSums(dat > 0) > min_detect
	filter_total <- colSums(dat) > min_counts

	plot(colSums(dat[rownames(dat) %in% pc,])/colSums(dat)*100, colSums(dat[rownames(dat) %in% mt,])/colSums(dat)*100, col=type_col[type], pch=source_pch[plate], xlab="%PC", ylab="%MT")

	plot(colSums(dat[rownames(dat) %in% pc,])/colSums(dat)*100, colSums(dat[rownames(dat) %in% rRNA,])/colSums(dat)*100, col=type_col[type], pch=source_pch[plate], xlab="%PC", ylab="%rRNA")

	plot(colSums(dat[rownames(dat) %in% pc,])/colSums(dat)*100, colSums(dat[rownames(dat) %in% sRNAs,])/colSums(dat)*100, col=type_col[type], pch=source_pch[plate], xlab="%PC", ylab="%sRNA")
	return(cbind(filter_detect, filter_total))
}

png("Camp_human_QC.png", width=7, height=7, units="in", res=300)
Filters <- my_QC(exprs(Camp), pData(Camp)$cell_type1, pData(Camp)$Source)
dev.off()
keep = Filters[,1] & Filters[,2]

Camp_QC <- Camp[,keep]
keep_genes <- rowSums(exprs(Camp_QC) > 0) > 2
Camp_QC <- Camp_QC[keep_genes,]


# Visualization
legend_plot <- function() {
	current = par("mar")
	names_col = c("definitive endo", "endothelial", "hepatic endo", "hepatoblast", "ipsc", "hepatocyte", "mesenchymal", "adult hepato", "erythroblast", "fetal hepato", "Kupffer", "lymphoblast", "stellate", "hepatic")
	cols=type_col
	names_pch = c("HUVEC","iPSC", "liver bud", "MSC", "Liver")
	pchs=source_pch
	par(mar=c(1,1,1,1))
	plot(1:10, col="white", axes=FALSE)
	legend("topleft", names_col, col=cols, pch=15, bty="n", title="Cell Type")
	legend("topright", names_pch, pch=pchs, col="black", bty="n", title="Source")
	par(mar=current)
}

# log2(CPM+1) normalization
renormed <- 2^exprs(Camp_QC)-1
sf <- colSums(renormed);
renormed <- t(t(renormed)/sf*median(sf))
unlogged <- renormed;
renormed <- log(renormed+1)/log(2);

#renormed <- exprs(Comp_QC)
# PCA

PCA <- prcomp(renormed)
# tSNE
require("Rtsne")
set.seed(123)
perp = 10
tmp <- renormed; rownames(tmp) <- 1:length(tmp[,1])
tmp <- t(tmp); tmp <- unique(tmp);
tsne_type <-pData(Camp_QC[,colnames(Camp_QC) %in% rownames(tmp)])$cell_type1
tsne_source <- pData(Camp_QC[,colnames(Camp_QC) %in% rownames(tmp)])$Source
tsne <- Rtsne(tmp, perplexity=perp)
tsne30 <- Rtsne(tmp, perplexity=30)
# DM
require("destiny")
set.seed(123)
DM <- DiffusionMap(t(tmp));

# cknn
source("~/ReviewArticle/functions.R")
require("igraph")
require("cccd")
cknn_mats <- cknn_matrix(t(renormed), 30)
cknn_clust <- C_agnostic_binary_search_clustering(cknn_mats$distance_matrix, cknn_mats$cknn_matrix,
                        Cs=c(1:min(dim(cknn_mats$cknn_matrix)[1]/20,50)), suppress.plot=FALSE)

cknn_graph <- cknn_clust$graph
coords <- layout_with_fr(cknn_graph)


png("Camp_human_Visualization.png", width=6*3/2, height=6, units="in", res=300)
par(mfrow=c(2,3))
par(mar=c(4,4,2,1))

plot(PCA$rotation[,1], PCA$rotation[,2], col=type_col[pData(Camp_QC)$cell_type1], pch=source_pch[pData(Camp_QC)$Source], xlab=paste("PC1 (",round(PCA$sdev[1]/sum(PCA$sdev)*100, digits=2)," %)",sep=""), ylab=paste("PC2 (",round(PCA$sdev[2]/sum(PCA$sdev)*100, digits=2)," %)",sep=""))
title(main="PCA")

plot(tsne$Y[,1], tsne$Y[,2], col=type_col[tsne_type], pch=source_pch[tsne_source], xlab="Component 1", ylab="Component 2")
title(main=paste("tSNE", perp))

plot(tsne30$Y[,1], tsne30$Y[,2], col=type_col[tsne_type], pch=source_pch[tsne_source], xlab="Component 1", ylab="Component 2")
title(main=paste("tSNE", 30))

plot(eigenvectors(DM)[,1], eigenvectors(DM)[,2], xlab = "Diffusion component 1", ylab = "Diffusion component 2", col=type_col[tsne_type], pch=source_pch[tsne_source])
title(main="DM")

plot(coords[,1], coords[,2], col=type_col[pData(Camp_QC)$cell_type1], pch=source_pch[pData(Camp_QC)$Source], xlab="", ylab="")
title(main="CkNN+FR")

legend_plot()

dev.off()

tsne_type <-pData(Camp_QC[,colnames(Camp_QC) %in% rownames(tmp)])$cell_type1

duplicates <-  which(!(colnames(Camp_QC) %in% rownames(tmp)))
reorder <- match(colnames(Camp_QC), rownames(tmp))

pData(Camp_QC)$PC1 <- PCA$rotation[,1]
pData(Camp_QC)$PC2 <- PCA$rotation[,2]
pData(Camp_QC)$tsne10_1 <- tsne$Y[reorder,1]
pData(Camp_QC)$tsne10_2 <- tsne$Y[reorder,2]
pData(Camp_QC)$tsne30_1 <- tsne30$Y[reorder,1]
pData(Camp_QC)$tsne30_2 <- tsne30$Y[reorder,2]
pData(Camp_QC)$DM1 <- eigenvectors(DM)[reorder,1]
pData(Camp_QC)$DM2 <- eigenvectors(DM)[reorder,2]
pData(Camp_QC)$DM3 <- eigenvectors(DM)[reorder,3]
pData(Camp_QC)$DM4 <- eigenvectors(DM)[reorder,4]
pData(Camp_QC)$CkNN_FR1 <- coords[,1]
pData(Camp_QC)$CkNN_FR2 <- coords[,2]
pData(Camp_QC)$CkNN_clusters <- components(cknn_graph)$membership

hasCol=function(A,a){colSums(a==A)==nrow(A)}
for( dub in duplicates) {
	same_cols <- which(hasCol(exprs(Camp_QC), exprs(Camp_QC)[,dub]))
	good_col <- same_cols[same_cols != dub]
	pData(Camp_QC)[dub,8:15] <- pData(Camp_QC)[good_col[1],8:length(pData(Camp_QC)[1,])]
}

# Feature Selection
require("M3Drop")
#unlogged <- 2^exprs(Camp_QC)-1

png("Camp_human_M3Drop_FS.png", width=6, height=6, units="in", res=300)
M3D_FS <- M3DropFeatureSelection(unlogged, mt_method="bonferroni", mt_threshold=0.05)
dev.off()

png("Camp_human_M3Drop_heatmap.png", width=8, height=8, units="in", res=300)
heatout <- M3DropExpressionHeatmap(M3D_FS, unlogged, cell_labels=pData(Camp_QC)$cell_type1)
dev.off()

M3D_table <- M3DropFeatureSelection(unlogged, mt_method="fdr", mt_threshold=2)
reorder <- match(rownames(fData(Camp_QC)), rownames(M3D_table))
fData(Camp_QC)$M3D_effect <- M3D_table[reorder,2]
fData(Camp_QC)$M3D_pval <- M3D_table[reorder,3]
fData(Camp_QC)$M3D_qval <- M3D_table[reorder,4]

new_counts <- round(2^exprs(Camp_QC)-1)
counts(Camp_QC) <- new_counts
exprs(Camp_QC) <- renormed


saveRDS(Camp_QC, file="Camp_QC.rds")

