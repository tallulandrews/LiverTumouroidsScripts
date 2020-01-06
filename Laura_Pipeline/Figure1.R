# Mega QCed SCE Object
require("SingleCellExperiment")
require("scater")
require("M3Drop")
CCA1 <- readRDS("CCA1_merged_SC3.rds") # Chol
CCA5 <- readRDS("CCA5_merged_SC3.rds") # Chol
HCC10 <- readRDS("HCC10_merged_SC3.rds") # Hepato
HCC23 <- readRDS("HCC23_merged_SC3.rds") # Btw
HCC24 <- readRDS("HCC24_merged_SC3.rds") # Hepato
HCC6 <- readRDS("HCC6_merged_SC3.rds") # Btw
D3DM <- readRDS("D3DM_merged_SC3.rds") # Btw
D3EM <- readRDS("D3EM_merged_SC3.rds") # Btw
D9DM <- readRDS("D9DM_merged_SC3.rds") # Btw
D9EM <- readRDS("D9EM_merged_SC3.rds") # Btw

SCEs <- list(CCA1=CCA1, CCA5=CCA5, HCC6=HCC6, HCC23=HCC23, HCC10=HCC10, HCC24=HCC24, D3DM=D3DM, D3EM=D3EM, D9DM=D9DM, D9EM=D9EM)
var_fs <- vector()
for (i in 1:length(SCEs)) {
	require("matrixStats")
	tmp <- SCEs[[i]]
	tmp <- tmp[rowData(tmp)$Symbol != "", ]
	tmp <- tmp[!is.na(rowData(tmp)$biotype) , ]
	tmp <- tmp[rowData(tmp)$biotype == "protein_coding", ]
	tmp <- tmp[order(rownames(tmp)), ]
	rowData(tmp)$vars <- rowVars(assays(tmp)[["lognorm"]])
	vs <- rownames(tmp)[order(rowData(tmp)$vars, decreasing=TRUE)]
	dat <- M3DropConvertData(tmp)
	fs <- M3DropFeatureSelection(dat, mt_method="fdr", mt_threshold=0.05, suppress.plot=TRUE)
	var_fs <- c(var_fs, fs$Gene[1:100])
	SCEs[[i]] <- tmp
}

consistent_genes <- rownames(SCEs[["CCA1"]])[rownames(SCEs[["CCA1"]]) %in% rownames(SCEs[["CCA5"]]) & rownames(SCEs[["CCA1"]]) %in% rownames(SCEs[["D3DM"]]) & rownames(SCEs[["CCA1"]]) %in% rownames(SCEs[["D9DM"]])]
consistent.symbol <- rowData(SCEs[["CCA1"]])$Symbol[rownames(SCEs[["CCA1"]]) %in% consistent_genes]

#Global <- cbind(assays(SCEs[["CCA1"]])[["lognorm"]][rownames(SCEs[["CCA1"]]) %in% consistent_genes,], 
#		assays(SCEs[["CCA5"]])[["lognorm"]][rownames(SCEs[["CCA5"]]) %in% consistent_genes,], 
#		assays(SCEs[["HCC6"]])[["lognorm"]][rownames(SCEs[["HCC6"]]) %in% consistent_genes,], 
#		assays(SCEs[["HCC23"]])[["lognorm"]][rownames(SCEs[["HCC23"]]) %in% consistent_genes,], 
#		assays(SCEs[["HCC10"]])[["lognorm"]][rownames(SCEs[["HCC10"]]) %in% consistent_genes,], 
#		assays(SCEs[["HCC24"]])[["lognorm"]][rownames(SCEs[["HCC24"]]) %in% consistent_genes,])

Global.counts <- cbind(assays(SCEs[["CCA1"]])[["counts"]][rownames(SCEs[["CCA1"]]) %in% consistent_genes,], 
		assays(SCEs[["CCA5"]])[["counts"]][rownames(SCEs[["CCA5"]]) %in% consistent_genes,], 
		assays(SCEs[["HCC6"]])[["counts"]][rownames(SCEs[["HCC6"]]) %in% consistent_genes,], 
		assays(SCEs[["HCC23"]])[["counts"]][rownames(SCEs[["HCC23"]]) %in% consistent_genes,], 
		assays(SCEs[["HCC10"]])[["counts"]][rownames(SCEs[["HCC10"]]) %in% consistent_genes,], 
		assays(SCEs[["HCC24"]])[["counts"]][rownames(SCEs[["HCC24"]]) %in% consistent_genes,],
		assays(SCEs[["D3DM"]])[["counts"]][rownames(SCEs[["D3DM"]]) %in% consistent_genes,],
		assays(SCEs[["D3EM"]])[["counts"]][rownames(SCEs[["D3EM"]]) %in% consistent_genes,],
		assays(SCEs[["D9DM"]])[["counts"]][rownames(SCEs[["D9DM"]]) %in% consistent_genes,],
		assays(SCEs[["D9EM"]])[["counts"]][rownames(SCEs[["D9EM"]]) %in% consistent_genes,]
		)
Global.means <- rowMeans(Global.counts)
Global.counts <- Global.counts[Global.means > 0,]
Global.symbol <- consistent.symbol[Global.means > 0]
Global.means <- Global.means[Global.means > 0]

size.factor <- colSums(Global.counts)
Global.norm <- t(t(Global.counts)/size.factor *median(size.factor))
Global.lognorm <- log2(Global.norm+1)

# Labels!
Donor <-  rep(c("CCA1", "CCA5", "HCC6", "HCC23", "HCC10", "HCC24", "D3DM", "D3EM", "D9DM", "D9EM"), times=c(ncol(CCA1), ncol(CCA5), ncol(HCC6), ncol(HCC23), ncol(HCC10), ncol(HCC24), ncol(D3DM), ncol(D3EM), ncol(D9DM), ncol(D9EM)))
Donor.col <- c("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", "#fdbf6f", "#ff7f00", "#cab2d6","#6a3d9a", "#fb9a99", "#e31a1c")
names(Donor.col) <- c("CCA1", "CCA5", "HCC6", "HCC23", "HCC10", "HCC24", "D3DM", "D3EM", "D9DM", "D9EM")
Type <-  rep(c("CC", "CC", "CHC", "CHC", "HCC", "HCC", "Ctrl", "Ctrl", "Ctrl", "Ctrl"), times=c(ncol(CCA1), ncol(CCA5), ncol(HCC6), ncol(HCC23), ncol(HCC10), ncol(HCC24), ncol(D3DM), ncol(D3EM), ncol(D9DM), ncol(D9EM)))
Type.col <- c("#1f78b4", "#33a02c", "#ff7f00", "#6a3d9a")
names(Type.col) <- c("CC", "CHC", "HCC", "Ctrl")

#Cycle <- c(as.character(CCA1$Cycle), as.character(CCA5$Cycle), 
#	   as.character(HCC6$Cycle), as.character(HCC23$Cycle), 
#	   as.character(HCC10$Cycle), as.character(HCC24$Cycle))
#Cycle[Cycle=="G0"] <- "None"
Cycle.col <- c("black", "#b15928", "#cab2d6")
names(Cycle.col) <- c("None", "G1S", "G2M")

Experiment <- rep(c("Exp1", "Exp2", "Exp1", "Exp2", "Exp1", "Exp2", "Exp3", "Exp3", "Exp4", "Exp4"), times=c(ncol(CCA1), ncol(CCA5), ncol(HCC6), ncol(HCC23), ncol(HCC10), ncol(HCC24), ncol(D3DM), ncol(D3EM), ncol(D9DM), ncol(D9EM)))
Experiment.col <- c("#f7f7f7", "#cccccc", "#969696", "#525252")
names(Experiment.col) <- c("Exp1", "Exp2", "Exp3", "Exp4")

Plate <- c(as.character(CCA1$Plate), as.character(CCA5$Plate), 
	   as.character(HCC6$Plate), as.character(HCC23$Plate), 
	   as.character(HCC10$Plate), as.character(HCC24$Plate), 
	   as.character(D3DM$Plate), as.character(D3EM$Plate), 
	   as.character(D9DM$Plate), as.character(D9EM$Plate))
Plate.col <- c("#fc8d62", "#66c2a5", "#8da0cb", "#a6d854", "#e78ac3", "#ffd92f", "#b3b3b3")
names(Plate.col) <- c("869", "870", "868", "3316", "3317", "3318", "3319")

# Visualize All together

require("M3Drop")
FS <- M3DropFeatureSelection(Global.norm, mt_method="fdr", mt_threshold=0.05, suppress.plot=FALSE)

# Global Cell-cycle
set.seed(1973)

require("CycleMix")
require("mclust")
source("~/NetworkInferencePipeline/Dropouts/My_R_packages/CycleMix/R/Analysis.R")
source("~/NetworkInferencePipeline/Dropouts/My_R_packages/CycleMix/R/Plotting.R")

#CC_genes <- rbind(HGeneSets$Tirosh, HGeneSets$Quiesc)
CC_genes <- rbind(HGeneSets$Tirosh)
SCE <- SingleCellExperiment(assays=list(lognorm=Global.lognorm, counts=Global.counts, norm=Global.norm))
source("~/R-Scripts/Ensembl_Stuff.R")
rowData(SCE)$Symbol <- General_Map(rownames(SCE), in.name="ensg", in.org="Hsap", out.name="symbol", out.org="Hsap")

out <- classifyCells(SCE, CC_genes, expr_name="lognorm", do.scale=FALSE, symbol_column="Symbol", allow.multi=FALSE)

colData(SCE)$CC_state <- factor(out$phase, levels=c("None", "G0", "G1S", "G2M"))
colData(SCE)$Proliferating <- ! (out$phase %in% c("G0", "None"))

colData(SCE)$Plate <- Plate
colData(SCE)$Type <- Type
colData(SCE)$Experiment <- Experiment
colData(SCE)$Donor <- Donor

saveRDS(SCE, "Global_SCE.rds")

# Make Plot
png("Supplementary_Fig1_ProlifMix.png", width=7, height=7, units="in", res=300)
par(mfrow=c(2,2))
par(mar=c(4,4,1,1))
for(i in out$fits) {
        plotMixture(i, BIC=FALSE)
}
dev.off()

Cycle <- as.character(colData(SCE)$CC_state)

# Lineage markers
#Chol_lineage <- read.table("/nfs/users/nfs_t/ta6/Collaborations/LiverOrganoids/Markers_130418_Chol.txt", header=TRUE)
#Hep_lineage <- read.table("/nfs/users/nfs_t/ta6/Collaborations/LiverOrganoids/Markers_130418_Hep.txt", header=TRUE)

#Hep_both <- Hep_lineage[ grepl("Prog", Hep_lineage[,2]) & grepl("Hep", Hep_lineage[,2]), 1]
#Chol_both <- Chol_lineage[ grepl("Prog", Chol_lineage[,2]) & grepl("Chol", Chol_lineage[,2]), 1]
#Prog_both <- Hep_lineage[ grepl("Prog", Hep_lineage[,2]) & Hep_lineage[,1] %in% Chol_lineage[Chol_lineage[,2] == "Prog",1], 1]

#Conflict1 <- Hep_lineage[ grepl("Hep", Hep_lineage[,2]) & Hep_lineage[,1] %in% Chol_lineage[Chol_lineage[,2] == "Chol",1], 1]
#Conflict2 <- Hep_lineage[ grepl("Prog", Hep_lineage[,2]) & Hep_lineage[,1] %in% Chol_lineage[Chol_lineage[,2] == "Chol",1], 1]
#Conflict3 <- Hep_lineage[ grepl("Hep", Hep_lineage[,2]) & Hep_lineage[,1] %in% Chol_lineage[Chol_lineage[,2] == "Prog",1], 1]
#Conflict4 <- Chol_lineage[ grepl("Prog", Chol_lineage[,2]) & Chol_lineage[,1] %in% Hep_lineage[Hep_lineage[,2] == "Hep",1], 1]
#
#Conflicts <- c(as.character(Conflict1), as.character(Conflict2), as.character(Conflict3), as.character(Conflicts4))

#Chol_lineage <- Chol_lineage[Chol_lineage[,1] %in% marker_genes & Chol_lineage[,1] %in% keep_genes,]
#Hep_lineage <- Hep_lineage[Hep_lineage[,1] %in% marker_genes & Hep_lineage[,1] %in% keep_genes,]

#Chol_lineage[,2] <- as.character(Chol_lineage[,2])
#Chol_lineage[Chol_lineage[,2] == "Prog",2] <- "Chol-Prog"
#Chol_lineage[Chol_lineage[,2] == "Chol",2] <- "Chol-Mature"
#Chol_lineage[Chol_lineage[,1] %in% Chol_both,2] <- "Chol-Both"
#Chol_lineage <- Chol_lineage[!(Chol_lineage[,1] %in% Conflicts),]

#Hep_lineage[,2] <- as.character(Hep_lineage[,2])
#Hep_lineage[Hep_lineage[,2] == "Prog",2] <- "Hep-Prog"
#Hep_lineage[Hep_lineage[,2] == "Hep",2] <- "Hep-Mature"
#Hep_lineage[Hep_lineage[,1] %in% Hep_both,2] <- "Hep-Both"
#Hep_lineage <- Hep_lineage[!(Hep_lineage[,1] %in% Conflicts),]

#Lineage <- rbind(Chol_lineage,Hep_lineage)
#Lineage[Lineage[,1] %in% Prog_both,2] <- "Common-Prog"
#Lineage <- Lineage[!duplicated(Lineage[,1]),]
#Lineage[,1] <- as.character(Lineage[,1])
#Lineage[Lineage[,1] == "05-Mar",1] <- "MARCH5"
#Lineage <- unique(Lineage)

#write.table(Lineage, file="Cleaned_Lineage.txt", row.names=F, col.names=F)

Lineage <- read.table("Cleaned_Lineage.txt", header=F, stringsAsFactors=FALSE)

# PCA
pca_data <- as.matrix(Global.lognorm[rownames(Global.lognorm) %in% FS$Gene[1:2000],])
rownames(pca_data) <- Global.symbol[rownames(Global.lognorm) %in% FS$Gene[1:2000]]


set.seed(3719869)
pca <- prcomp(pca_data)
png("Figure1_pcaSimple.png", width=5*2, height=5, units="in", res=300)
par(mfrow=c(1,2))
par(mar=c(4,4,1,1))
plot(pca$rotation[,1], pca$rotation[,2], pch=16, col=Donor.col[factor(Donor, levels=names(Donor.col))], 
	xlab=paste("PC1 (",round(pca$sdev[1]^2/sum(pca$sdev^2)*100, digits=1), "%)", sep=""),
	ylab=paste("PC2 (",round(pca$sdev[2]^2/sum(pca$sdev^2)*100, digits=1), "%)", sep=""))
legend("bottomleft", names(Donor.col), pch=16, col=Donor.col, bty="n", title="Donor")
plot(pca$rotation[,1], pca$rotation[,2], pch=16, col=Plate.col[factor(Plate, levels=names(Plate.col))], 
	xlab=paste("PC1 (",round(pca$sdev[1]^2/sum(pca$sdev^2)*100, digits=1), "%)", sep=""),
	ylab=paste("PC2 (",round(pca$sdev[2]^2/sum(pca$sdev^2)*100, digits=1), "%)", sep=""))
legend("bottomleft", as.character(1:length(Plate.col)), pch=16, col=Plate.col, bty="n", title="Plate")
dev.off()

plot(pca$rotation[,3], pca$rotation[,4], pch=16, col=Donor.col[factor(Donor, levels=names(Donor.col))], 
	xlab=paste("PC3 (",round(pca$sdev[3]^2/sum(pca$sdev^2)*100, digits=1), "%)", sep=""),
	ylab=paste("PC4 (",round(pca$sdev[4]^2/sum(pca$sdev^2)*100, digits=1), "%)", sep=""))


png("Supplementary_Fig1_PCAcutoff.png", width=5, height=5, units="in", res=300)
par(mar=c(4,4,1,1))
plot(pca$sdev^2/sum(pca$sdev^2)*100, type="b", xlab="Principal Component", ylab="Variance (%)", xlim=c(1,20))
cutoff <- 7
abline(v=cutoff+0.5, col="grey45", lty=3)
dev.off()

head(colMeans(pca$x[rownames(pca_data)%in% Lineage[Lineage[,2]=="Hep-Mature",1],]), cutoff)
head(colMeans(pca$x[rownames(pca_data)%in% Lineage[Lineage[,2]=="Chol-Mature",1],]), cutoff)
head(colMeans(pca$x[rownames(pca_data)%in% Lineage[grepl("Prog", Lineage[,2]),1],]), cutoff)
barplot_dat <-rbind(colMeans(pca$x[rownames(pca_data)%in% Lineage[Lineage[,2]=="Chol-Mature",1],1:cutoff]),
		    colMeans(pca$x[rownames(pca_data)%in% Lineage[grepl("Prog", Lineage[,2]),1],1:cutoff]),
		    colMeans(pca$x[rownames(pca_data)%in% Lineage[Lineage[,2]=="Hep-Mature",1],1:cutoff]))

png("Supplementary_Fig1_PCAlinScores.png", width=5, height=5, units="in", res=300)
par(mar=c(4,4,1,1))
barplot(barplot_dat, beside=T, col=Type.col[1:3], ylab="Avg Marker Weights")
abline(h=c(7.5,-7.5), col="grey50", lty=2)
legend("bottomright", c("Chol = -PC1-PC4", "Stem = -PC1 -PC2 +PC3 +PC4", "Hep = +PC1 -PC2 +PC4"), fill=Type.col, bty="n")
dev.off()

# from the top 10 components take those with average weights > 10
thing <- apply(pca$rotation, 2, scale)
chol_score <- (-thing[,1]-thing[,4])/2
hep_score <- (thing[,1]-thing[,2]-thing[,4])/3
stem_score <- (-thing[,1]-thing[,2]+thing[,4]+thing[,3])/4

xes <- hep_score-chol_score
#xes[hep_score>chol_score] <- hep_score[hep_score>chol_score]
yes <- stem_score


SCE$Hep_score <- hep_score
SCE$Chol_score <- chol_score
SCE$Stem_score <- stem_score


png("Figure1_Lineage.png", width=5*2, height=5, units="in", res=300)
par(mfrow=c(1,2))
par(mar=c(4,4,1,1))
plot(xes, yes, xlab="Cholangiocyte ---------- Hepatocyte", ylab="Stemness", pch=16, col=Donor.col[factor(Donor, levels=names(Donor.col))])
legend("topleft", names(Donor.col), col=Donor.col, pch=16, bty="n", ncol=2)
top_x <- mean(xes[Donor %in% c("HCC23", "HCC24")])
top_y <- mean(yes[Donor %in% c("HCC23", "HCC24")])
l_x <-  mean(xes[Donor %in% c("CCA5")])
l_y <-  mean(yes[Donor %in% c("CCA5")])
r_x <-  mean(xes[Donor %in% c("HCC10")])
r_y <-  mean(yes[Donor %in% c("HCC10")])
#arrows(top_x, top_y, l_x, l_y, len=0, lwd=2, col="black")
#arrows(top_x, top_y, r_x, r_y, len=0, lwd=2, col="black")


plot(xes, yes, xlab="Cholangiocyte ---------- Hepatocyte", ylab="Stemness", pch=16, col=Cycle.col[factor(Cycle, levels=names(Cycle.col))])
legend("topleft", names(Cycle.col), col=Cycle.col, pch=16, bty="n")
dev.off()


#plot(pca$rotation[,1], pca$rotation[,2], pch=16, col=Donor.col[factor(Donor, levels=names(Donor.col))])
# component 5 = stem-ness/cell-cycle, component 7 = semm-ness/cell-cycle, components 1 & 2 = CC vs HCC



# t-SNE
require("Rtsne")
set.seed(134)
TSNE <- Rtsne(t(Global.lognorm[rownames(Global.lognorm) %in% FS$Gene[1:2000],]), initial_dims=cutoff, dims=2, perplexity=30)
png("Figure1_tSNE.png", width=5*2, height=5, units="in", res=300)
par(mfrow=c(1,2))
par(mar=c(4,4,1,1))
plot(TSNE$Y[,1], TSNE$Y[,2],  pch=16, col=Donor.col[factor(Donor, levels=names(Donor.col))], xlab="Dim 1", ylab="Dim 2")
legend("bottomleft", names(Donor.col), pch=16, col=Donor.col, bty="n", ncol=2)
plot(TSNE$Y[,1], TSNE$Y[,2],  pch=16, col=Cycle.col[factor(Cycle, levels=names(Cycle.col))], xlab="Dim 1", ylab="Dim 2")
legend("bottomleft", names(Cycle.col), pch=16, col=Cycle.col, bty="n")
dev.off()




# Gene Heatmap
#M3DropExpressionHeatmap(fs$Gene[1:1000], Global.norm, cell_labels=Donor)
source("~/R-Scripts/heatmap.3.R")
require("RColorBrewer")
require("gplots")

heatcolours <- rev(brewer.pal(11, "RdBu"))
col_breaks <- c(-100, seq(-2, 2, length = 10), 100)
heat_data <- as.matrix(Global.lognorm[rownames(Global.lognorm) %in% FS$Gene[1:2000],])
rownames(heat_data) <- Global.symbol[rownames(Global.lognorm) %in% FS$Gene[1:2000]]
heat_data <- heat_data[rownames(heat_data) %in% Lineage[,1],]
ColColors1 <- Donor.col[factor(Donor, levels=names(Donor.col))]
ColColors2 <- Experiment.col[factor(Experiment, levels=names(Experiment.col))]
lin_label <- Lineage[match(rownames(heat_data), Lineage[,1]),2]
#lin_label[grep("Prog", lin_label)] <- "Prog"
RowColors <- c(Donor.col["CCA5"], Donor.col["CCA1"], Donor.col["HCC23"], Donor.col["HCC10"], Donor.col["HCC24"])[factor(lin_label, levels=c("Chol-Mature", "Chol-Prog", "Common-Prog", "Hep-Prog", "Hep-Mature"))]
RowColors[is.na(RowColors)] <- "white"

markers <- c("NDRG1", "CYP2E1", "ADH6", "ADH4", "KRT7", "CFTR", "ZWINT", "KRT19", "HNF4A", "CDK1", "TOP2A", "CCNB2")
gene_lab <- rep("", nrow(heat_data))
gene_lab[rownames(heat_data) %in% markers] <- rownames(heat_data)[rownames(heat_data) %in% markers]

lwid <- c(1, 0.2, 4)
lhei <- c(1, 0.2, 4)
lmat <- rbind(c(6, 0, 5), c(0, 0, 2), c(4, 1, 3))
require("gplots")
png("Figure1_heatmap.png", width=6*2, height=6, units="in", res=300)
heatmap_output <- suppressWarnings(heatmap.2(heat_data, 
            ColSideColors = ColColors1, RowSideColors = RowColors, 
            col = heatcolours, breaks = col_breaks, scale = "row", 
            symbreaks = TRUE, trace = "none", dendrogram = "both", #labRow = gene_lab,
            key = FALSE, Rowv = TRUE, Colv = TRUE, lwid = lwid, 
            lhei = lhei, lmat = lmat, hclustfun = function(x) {
                hclust(x, method = "ward.D2")
            }))
dev.off()

# Gene sets?
source("~/NetworkInferencePipeline/Dropouts/My_R_packages/M3D/R/Plotting_fxns.R")
g_groups <- M3DropGetHeatmapClusters(heatmap_output, k=14, type="gene")
RowColors<- rainbow(14)[g_groups]
g_groups[g_groups == 4] <- 3
g_groups[g_groups == 5] <- 6
g_groups[g_groups == 11] <- 12
RowColors<- rainbow(14)[g_groups]

heatmap_output <- suppressWarnings(heatmap.2(heat_data,
            ColSideColors = ColColors1, RowSideColors = RowColors,
            col = heatcolours, breaks = col_breaks, scale = "row",
            symbreaks = TRUE, trace = "none", dendrogram = "both", #labRow = gene_lab,
            key = FALSE, Rowv = TRUE, Colv = TRUE, lwid = lwid,
            lhei = lhei, lmat = lmat, hclustfun = function(x) {
                hclust(x, method = "ward.D2")
            }))




SCE <- regressCycle_partial(SCE, out, expr_name="lognorm", method="phase", phases=c("G1S", "G2M"))
colData(SCE)$tsne_x <- TSNE$Y[,1]
colData(SCE)$tsne_y <- TSNE$Y[,2]
pca_rot <- pca$rotation[,1:cutoff]
colnames(pca_rot) <- paste("PCA", 1:cutoff, sep="")
colData(SCE) <- cbind(colData(SCE), pca_rot);


pca_wei <- pca$x[,1:cutoff]
colnames(pca_wei) <- paste("PCA", 1:cutoff, sep="")
thing <- match(rowData(SCE)$Symbol, rownames(pca_wei))
rowData(SCE) <- cbind(rowData(SCE), pca_wei[thing,]);

saveRDS(SCE, "Global_SCE.rds")

# For Montreal 
SCE <- readRDS("Global_SCE.rds");
Lineage <- read.table("Cleaned_Lineage.txt", header=F, stringsAsFactors=FALSE)

Chol_Control <- SCE[,SCE$Donor %in% c("D3EM", "D3DM", "D9EM", "D9DM")]
SCE <- SCE[,!SCE$Donor %in% c("D3EM", "D3DM", "D9EM", "D9DM")]

# Colour Schemes
Donor.col <- c("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", "#fdbf6f", "#ff7f00", "#cab2d6","#6a3d9a", "#fb9a99", "#e31a1c")
names(Donor.col) <- c("CCA1", "CCA5", "HCC6", "HCC23", "HCC10", "HCC24", "D3DM", "D3EM", "D9DM", "D9EM")
Type.col <- c("#1f78b4", "#33a02c", "#ff7f00", "#6a3d9a")
names(Type.col) <- c("CC", "CHC", "HCC", "Ctrl")
Cycle.col <- c("black", "#b15928", "#cab2d6")
names(Cycle.col) <- c("None", "G1S", "G2M")
Experiment.col <- c("#f7f7f7", "#cccccc", "#969696", "#525252")
names(Experiment.col) <- c("Exp1", "Exp2", "Exp3", "Exp4")

require("M3Drop")
FS <- M3DropFeatureSelection(assays(SCE)[["norm"]], mt_method="fdr", mt_threshold=0.05, suppress.plot=FALSE)

pca_data <- as.matrix(assays(SCE)[["lognorm"]][rownames(SCE) %in% FS$Gene[1:2000],])
rownames(pca_data) <- rowData(SCE)$Symbol[rownames(SCE) %in% FS$Gene[1:2000]]


set.seed(3719869)
pca <- prcomp(pca_data)



plot(pca$rotation[,1], pca$rotation[,2], pch=16, col=Donor.col[factor(Donor, levels=names(Donor.col))],
        xlab=paste("PC1 (",round(pca$sdev[1]^2/sum(pca$sdev^2)*100, digits=1), "%)", sep=""),
        ylab=paste("PC2 (",round(pca$sdev[2]^2/sum(pca$sdev^2)*100, digits=1), "%)", sep=""))
legend("bottomleft", names(Donor.col), pch=16, col=Donor.col, bty="n", title="Donor")

thing <- apply(pca$rotation, 2, scale)
chol_score <- (-thing[,1])
stem_score <- (thing[,5]-thing[,2]-thing[,3])/3
hep_score <- (thing[,1]-thing[,2]+thing[,3])/3
