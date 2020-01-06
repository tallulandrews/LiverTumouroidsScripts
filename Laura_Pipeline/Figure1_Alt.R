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

SCEs <- list(CCA1=CCA1, CCA5=CCA5, HCC6=HCC6, HCC23=HCC23, HCC10=HCC10, HCC24=HCC24)
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

consistent_genes <- rownames(SCEs[["CCA1"]])[rownames(SCEs[["CCA1"]]) %in% rownames(SCEs[["CCA5"]])]
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
		assays(SCEs[["HCC24"]])[["counts"]][rownames(SCEs[["HCC24"]]) %in% consistent_genes,]
		)
Global.means <- rowMeans(Global.counts)
Global.counts <- Global.counts[Global.means > 0,]
Global.symbol <- consistent.symbol[Global.means > 0]
Global.means <- Global.means[Global.means > 0]

size.factor <- colSums(Global.counts)
Global.norm <- t(t(Global.counts)/size.factor *median(size.factor))
Global.lognorm <- log2(Global.norm+1)

# Labels!
Donor <-  rep(c("CCA1", "CCA5", "HCC6", "HCC23", "HCC10", "HCC24"), times=c(ncol(CCA1), ncol(CCA5), ncol(HCC6), ncol(HCC23), ncol(HCC10), ncol(HCC24)))
Donor.col <- c("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", "#fdbf6f", "#ff7f00")
names(Donor.col) <- c("CCA1", "CCA5", "HCC6", "HCC23", "HCC10", "HCC24")
Type <-  rep(c("CC", "CC", "CHC", "CHC", "HCC", "HCC"), times=c(ncol(CCA1), ncol(CCA5), ncol(HCC6), ncol(HCC23), ncol(HCC10), ncol(HCC24)))
Type.col <- c("#1f78b4", "#33a02c", "#ff7f00")
names(Type.col) <- c("CC", "CHC", "HCC")

#Cycle <- c(as.character(CCA1$Cycle), as.character(CCA5$Cycle), 
#	   as.character(HCC6$Cycle), as.character(HCC23$Cycle), 
#	   as.character(HCC10$Cycle), as.character(HCC24$Cycle))
#Cycle[Cycle=="G0"] <- "None"
Cycle.col <- c("black", "#b15928", "#cab2d6")
names(Cycle.col) <- c("None", "G1S", "G2M")

Experiment <- rep(c("Exp1", "Exp2", "Exp1", "Exp2", "Exp1", "Exp2"), times=c(ncol(CCA1), ncol(CCA5), ncol(HCC6), ncol(HCC23), ncol(HCC10), ncol(HCC24)))
Experiment.col <- c("#f7f7f7", "#cccccc")
names(Experiment.col) <- c("Exp1", "Exp2")

Plate <- c(as.character(CCA1$Plate), as.character(CCA5$Plate), 
	   as.character(HCC6$Plate), as.character(HCC23$Plate), 
	   as.character(HCC10$Plate), as.character(HCC24$Plate))
Plate.col <- c("#fc8d62", "#66c2a5", "#8da0cb", "#a6d854", "#e78ac3", "#ffd92f", "#b3b3b3")
names(Plate.col) <- c("869", "870", "868", "3316", "3317", "3318", "3319")

# Visualize All together

require("M3Drop")
FS <- M3DropFeatureSelection(Global.norm, mt_method="fdr", mt_threshold=0.05, suppress.plot=FALSE)
#FS2 <- BrenneckeGetVariableGenes(Global.norm)

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

rowData(SCE)$FStop2000 <- rownames(SCE) %in% FS$Gene[1:2000]

saveRDS(SCE, "Global_SCE_Alt.rds")

# Make Plot
png("Supplementary_Fig1_ProlifMix_Alt.png", width=7, height=7, units="in", res=300)
par(mfrow=c(2,2))
par(mar=c(4,4,1,1))
for(i in out$fits) {
        plotMixture(i, BIC=FALSE)
}
dev.off()

# Cycle by donor
Cycle <- as.character(colData(SCE)$CC_state)
Cycle <- factor(Cycle, levels=names(Cycle.col))
Donor <- as.character(colData(SCE)$Donor)
t <- table(Cycle, Donor);
t <- t( t(t)/colSums(t)*100 )
t <- t[,c(2,1,6,4,5,3)]

png("Supplementary_Fig1_ProlifDonor_Alt.png", width=5.5, height=5.5, units="in", res=300)
par(mar=c(4,4,1,1))
barplot(t, col=Cycle.col, las=2, ylab="Cells (%)")
dev.off()




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
pca_data <- as.matrix(Global.lognorm[rownames(Global.lognorm) %in% FS$Gene[1:1500],])
rownames(pca_data) <- Global.symbol[rownames(Global.lognorm) %in% FS$Gene[1:1500]]


set.seed(3719869)
pca <- prcomp(pca_data)
png("Figure1_pcaSimple_Alt.png", width=5*2, height=5, units="in", res=300)
par(mfrow=c(1,2))
par(mar=c(4,4,1,1))
plot(pca$rotation[,1], pca$rotation[,2], pch=16, col=Donor.col[factor(Donor, levels=names(Donor.col))], 
	xlab=paste("PC1 (",round(pca$sdev[1]^2/sum(pca$sdev^2)*100, digits=1), "%)", sep=""),
	ylab=paste("PC2 (",round(pca$sdev[2]^2/sum(pca$sdev^2)*100, digits=1), "%)", sep=""))
legend("bottomright", names(Donor.col), pch=16, col=Donor.col, bty="n", title="Donor")
plot(pca$rotation[,1], pca$rotation[,2], pch=16, col=Plate.col[factor(Plate, levels=names(Plate.col))], 
	xlab=paste("PC1 (",round(pca$sdev[1]^2/sum(pca$sdev^2)*100, digits=1), "%)", sep=""),
	ylab=paste("PC2 (",round(pca$sdev[2]^2/sum(pca$sdev^2)*100, digits=1), "%)", sep=""))
legend("bottomright", as.character(1:length(Plate.col)), pch=16, col=Plate.col, bty="n", title="Plate")
dev.off()

plot(pca$rotation[,3], pca$rotation[,4], pch=16, col=Donor.col[factor(Donor, levels=names(Donor.col))], 
	xlab=paste("PC3 (",round(pca$sdev[3]^2/sum(pca$sdev^2)*100, digits=1), "%)", sep=""),
	ylab=paste("PC4 (",round(pca$sdev[4]^2/sum(pca$sdev^2)*100, digits=1), "%)", sep=""))


png("Supplementary_Fig1_PCAcutoff_Alt.png", width=5, height=5, units="in", res=300)
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

png("Supplementary_Fig1_PCAlinScores_Alt.png", width=5, height=5, units="in", res=300)
par(mar=c(4,4,1,1))
barplot(barplot_dat, beside=T, col=Type.col[1:3], ylab="Avg Marker Weights")
abline(h=c(7.5,-7.5), col="grey50", lty=2)
legend("bottomright", c("Chol = PC1", "Stem = -PC2 -PC3 +PC5", "Hep = -PC2 +PC3"), fill=Type.col, bty="n")
dev.off()

# from the top 10 components take those with average weights > 10
thing <- apply(pca$rotation, 2, scale)
chol_score <- (thing[,1])/1
hep_score <- (-thing[,2]+thing[,3])/2
stem_score <- (thing[,5]-thing[,2]-thing[,3])/3

xes <- hep_score-chol_score
#xes[hep_score>chol_score] <- hep_score[hep_score>chol_score]
yes <- stem_score


SCE$Hep_score <- hep_score
SCE$Chol_score <- chol_score
SCE$Stem_score <- stem_score


png("Figure1_Lineage_Alt.png", width=5*2, height=5, units="in", res=300)
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
png("Figure1_tSNE_Alt.png", width=5*2, height=5, units="in", res=300)
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
png("Figure1_heatmap_Alt.png", width=6*2, height=6, units="in", res=300)
heatmap_output <- suppressWarnings(heatmap.2(heat_data, 
            ColSideColors = ColColors1, RowSideColors = RowColors, 
            col = heatcolours, breaks = col_breaks, scale = "row", 
            symbreaks = TRUE, trace = "none", dendrogram = "both", #labRow = gene_lab,
            key = FALSE, Rowv = TRUE, Colv = TRUE, lwid = lwid, 
            lhei = lhei, lmat = lmat, hclustfun = function(x) {
                hclust(x, method = "ward.D2")
            }))
dev.off()


require("CellTypeProfiles")
M <- complex_markers(heat_data, Donor)
M[M$CCA1==1& M$CCA5==1& M$HCC6==1 & RowColors=="#1f78b4" & M$HCC23==0 & M$HCC24==0,]
M[M$HCC10==1& M$CCA5==0& M$HCC6==0 & M$CCA1==0 & M$HCC23==0 & M$HCC24==0 & RowColors=="#ff7f00",]
M[M$HCC10==0& M$CCA5==0& M$HCC6==0 & M$CCA1==0 & M$HCC23==1 & M$HCC24==1 & RowColors%in% c("#33a02c", "#a6cee3", "#fdbf6f") ,]

#Known_markers_Chol <- c("GPX2", "CYP2S1", "ELF3", "LAD1", "ITGB42","TJP3","TSPO")
Known_markers_Chol <- rownames(M[RowColors==Donor.col["CCA5"] & M$AUC > 0.9 & M$CCA1+M$CCA5+M$HCC6 > 1 & M$HCC10==0,])
#Known_markers_Hep <- c("AQP9", "APOA5", "ADH4", "PAH", "ADH6", "CYP2E1")
Known_markers_Hep <- rownames(M[RowColors==Donor.col["HCC24"] & M$AUC > 0.9 & M$HCC10== 1 & M$CCA1+M$CCA5+M$HCC6 < 2,])
#Known_markers_Stem <- c("LAMB1", "TYRO3", "SALL4", "SPARC", "ROR2", "IFI30")
Known_markers_Stem <- rownames(M[RowColors%in% c(Donor.col["CCA1"], Donor.col["HCC23"], Donor.col["HCC10"]) & M$AUC > 0.9 & M$HCC23+M$HCC24>= 1,])

G_to_use <- c(Known_markers_Chol, Known_markers_Hep, Known_markers_Stem)
G_to_use <- G_to_use[-grep("-", G_to_use)]

heat_data2 <- as.matrix(Global.lognorm[Global.symbol %in% G_to_use,])
rownames(heat_data2) <- Global.symbol[Global.symbol %in% G_to_use]
ColColors1 <- Donor.col[factor(Donor, levels=names(Donor.col))]
ColColors2 <- Experiment.col[factor(Experiment, levels=names(Experiment.col))]
lin_label <- Lineage[match(rownames(heat_data2), Lineage[,1]),2]
RowColors2 <- c("forestgreen", "darkgoldenrod", "darkgoldenrod", "darkgoldenrod", "firebrick")[factor(lin_label, levels=c("Chol-Mature", "Chol-Prog", "Common-Prog", "Hep-Prog", "Hep-Mature"))]
RowColors2[is.na(RowColors2)] <- "white"

png("Figure1_heatmap2_Alt.png", width=6*2, height=6, units="in", res=300)
heatmap_output <- suppressWarnings(heatmap.2(t(heat_data2), 
            ColSideColors = RowColors2, RowSideColors = ColColors1, 
            col = heatcolours, breaks = col_breaks, scale = "row", 
            symbreaks = TRUE, trace = "none", dendrogram = "both", #labRow = gene_lab,
            key = FALSE, Rowv = TRUE, Colv = TRUE, lwid = lwid, 
            lhei = lhei, lmat = lmat, hclustfun = function(x) {
                hclust(x, method = "complete")
            }))
dev.off()

### Novel CC markers ###
M2 <- complex_markers(Global.lognorm, Donor)
novel_chol <- M2[ M2$CCA1 == 1 & M2$CCA5 == 1 & M2$HCC6 == 1 & M2$HCC10==0 & M2$HCC23==0 & M2$HCC24==0 & !Global.symbol %in% Lineage[,1],]

D3DM <- readRDS("D3DM_merged_SC3.rds") # Btw
D3EM <- readRDS("D3EM_merged_SC3.rds") # Btw
D9DM <- readRDS("D9DM_merged_SC3.rds") # Btw
D9EM <- readRDS("D9EM_merged_SC3.rds") # Btw


d3dm_pbulk <- rowMeans(assays(D3DM)[["lognorm"]])
d3em_pbulk <- rowMeans(assays(D3EM)[["lognorm"]])
d9dm_pbulk <- rowMeans(assays(D9DM)[["lognorm"]])
d9em_pbulk <- rowMeans(assays(D9EM)[["lognorm"]])
cca1_pbulk <- rowMeans(assays(CCA1)[["lognorm"]])
cca5_pbulk <- rowMeans(assays(CCA5)[["lognorm"]])
hcc6_pbulk <- rowMeans(assays(HCC6)[["lognorm"]])
genes <- names(cca1_pbulk)[names(cca1_pbulk) %in% names(cca5_pbulk) & names(cca1_pbulk) %in% names(d3dm_pbulk) & names(cca1_pbulk) %in% names(d9dm_pbulk)]

tab <- cbind(d3dm_pbulk[names(d3dm_pbulk) %in% genes], 
	     d3em_pbulk[names(d3em_pbulk) %in% genes], 
	     d9em_pbulk[names(d9em_pbulk) %in% genes], 
	     d9dm_pbulk[names(d9dm_pbulk) %in% genes], 
	     cca1_pbulk[names(cca1_pbulk) %in% genes], 
	     cca5_pbulk[names(cca5_pbulk) %in% genes], 
	     hcc6_pbulk[names(hcc6_pbulk) %in% genes])

require("CycleMix")
out <- vector()
pro_mean <- vector()
qui_mean <- vector()
for (d in unique(Donor)) {
	subset <- SCE[,Donor==d]
	de <- apply(assays(subset)[["lognorm"]], 1, function(x) {wilcox.test(x[subset$Proliferating], x[!subset$Proliferating])$p.value})
	out <- cbind(out, de);
	pro_mean <- cbind(pro_mean, rowMeans(assays(subset)[["lognorm"]][,subset$Proliferating]))
	qui_mean <- cbind(qui_mean, rowMeans(assays(subset)[["lognorm"]][,!subset$Proliferating]))
}
colnames(out) <- unique(Donor);
colnames(pro_mean) <- unique(Donor);
colnames(qui_mean) <- unique(Donor);
sig_cancer <- apply(out, 2, function(y){p.adjust(y, method="fdr") < 0.05})
consistent <- rowSums(sig_cancer, na.rm=T)==6

normal <- list(D3DM, D3EM, D9DM, D9EM)
ctrl_p <- vector()
ctrl_pro <- vector()
ctrl_qui <- vector()
for (i in normal) {
	de <- apply(assays(i)[["lognorm"]], 1, function(x) {wilcox.test(x[i$Proliferating], x[!i$Proliferating])$p.value})
	ctrl_p <- cbind(ctrl_p, de);
	ctrl_pro <- cbind(ctrl_pro, rowMeans(assays(i)[["lognorm"]][,i$Proliferating]))
	ctrl_qui <- cbind(ctrl_qui, rowMeans(assays(i)[["lognorm"]][,!i$Proliferating]))
}

colnames(ctrl_p) <- c("D3DM", "D3EM", "D9DM", "D9EM")
colnames(ctrl_pro) <- c("D3DM", "D3EM", "D9DM", "D9EM")
colnames(ctrl_qui) <- c("D3DM", "D3EM", "D9DM", "D9EM")

sig_norm <- apply(ctrl_p, 2, function(y){p.adjust(y, method="fdr") < 0.05})
multi <- rowSums(sig_norm, na.rm=T)>2
sig_norm <- sig_norm[match(rownames(sig_cancer), rownames(sig_norm)),]


sum( rowSums(sig_cancer, na.rm=T)>=6 & rowSums(sig_norm, na.rm=T) < 2 )

qui_mean2<-qui_mean
qui_mean2[qui_mean2==0] <- 1/400
ctrl_qui2<-ctrl_qui
ctrl_qui2[ctrl_qui2==0] <- 1/400
l2fc_cancer <- pro_mean/qui_mean2
l2fc_norm <- ctrl_pro/ctrl_qui2
l2fc_norm<- l2fc_norm[match(rownames(l2fc_cancer), rownames(l2fc_norm)),]

rownames(sig_cancer) <- Global.symbol
rownames(sig_norm) <- Global.symbol
rownames(l2fc_cancer) <- Global.symbol
rownames(l2fc_norm) <- Global.symbol


require("RColorBrewer")
GO <- read.delim("/lustre/scratch117/cellgen/team218/TA/Mirrored_Annotations/GO/hsapiens_80_GO_Annoations_Emsembl.out", sep="\t", header=F)
GO_cc <- as.character(GO[GO[,3] == "cell cycle",1])
source("~/R-Scripts/Ensembl_Stuff.R")
GO_cc<-General_Map(GO_cc, in.org="Hsap", out.org="Hsap", in.name="ensg", out.name="symbol")

heat_cols <- brewer.pal(8,"Reds")
top <- rownames(sig_cancer)[rowSums(sig_cancer, na.rm=T)>=6]# & rowSums(sig_norm, na.rm=T) < 2]
top <- top[!top %in% c(HGeneSets$Tirosh$Gene, HGeneSets$Macosko$Gene, HGeneSets$Whitfield$Gene, GO_cc)]
effect <- cbind(l2fc_cancer)#, l2fc_norm)
good <- rownames(effect)[rowMeans(effect) > 2]
require("gplots")
png("Consistent_Prolif_Alt.png", width=6, height=9, units="in", res=300)
heatmap.2(log2(effect[rownames(effect) %in% top & rownames(effect) %in% good,]), trace="none", col=heat_cols)
dev.off()
write.table(effect[rownames(effect) %in% top & rownames(effect) %in% good,], file="Consistent_Prolif_notCC_Alt.txt", row.names=TRUE, col.names=TRUE)



# Save results
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

saveRDS(SCE, "Global_SCE_Alt.rds")

