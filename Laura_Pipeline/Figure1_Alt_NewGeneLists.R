# Mega QCed SCE Object
require("SingleCellExperiment")
require("scater")
require("M3Drop")

SCE <- readRDS("Global_SCE_Alt.rds")
Lineage <- read.table("~/Collaborations/LiverOrganoids/Cleaned_PatLineage.txt", header=F, stringsAsFactors=FALSE)
# Labels!
Donor.col <- c("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", "#fdbf6f", "#ff7f00")
names(Donor.col) <- c("CCA1", "CCA5", "HCC6", "HCC23", "HCC10", "HCC24")
Type.col <- c("#1f78b4", "#33a02c", "#ff7f00")
names(Type.col) <- c("CC", "CHC", "HCC")
Cycle.col <- c("black", "#b15928", "#cab2d6")
names(Cycle.col) <- c("None", "G1S", "G2M")
Experiment.col <- c("#f7f7f7", "#cccccc")
names(Experiment.col) <- c("Exp1", "Exp2")
Plate.col <- c("#fc8d62", "#66c2a5", "#8da0cb", "#a6d854", "#e78ac3", "#ffd92f", "#b3b3b3")
names(Plate.col) <- c("869", "870", "868", "3316", "3317", "3318", "3319")

# Visualize All together
require("M3Drop")
png("Figure1_M3Drop_FS.png", width=5, height=5, units="in", res=300)
FS <- M3DropFeatureSelection(assays(SCE)[["norm"]], mt_method="fdr", mt_threshold=0.05, suppress.plot=FALSE)
dev.off();


# PCA
pca_data <- as.matrix(assays(SCE)[["lognorm"]][rownames(SCE) %in% FS$Gene[1:1500],])
rownames(pca_data) <- rowData(SCE)$Symbol[rownames(SCE) %in% FS$Gene[1:1500]]


set.seed(3719869)
pca <- prcomp(pca_data)
#png("Supplementary_Fig1_PCAcutoff_Alt.png", width=5, height=5, units="in", res=300)
#par(mar=c(4,4,1,1))
plot(pca$sdev^2/sum(pca$sdev^2)*100, type="b", xlab="Principal Component", ylab="Variance (%)", xlim=c(1,20))
cutoff <- 7
abline(v=cutoff+0.5, col="grey45", lty=3)
#dev.off()

head(colMeans(pca$x[rownames(pca_data)%in% Lineage[Lineage[,2]=="Hep-Mature",1],]), cutoff)
head(colMeans(pca$x[rownames(pca_data)%in% Lineage[Lineage[,2]=="Chol-Mature",1],]), cutoff)
head(colMeans(pca$x[rownames(pca_data)%in% Lineage[grepl("Prog", Lineage[,2]),1],]), cutoff)
barplot_dat <-rbind(colMeans(pca$x[rownames(pca_data)%in% Lineage[Lineage[,2]=="Chol-Mature",1],1:cutoff]),
		    colMeans(pca$x[rownames(pca_data)%in% Lineage[grepl("Prog", Lineage[,2]),1],1:cutoff]),
		    colMeans(pca$x[rownames(pca_data)%in% Lineage[Lineage[,2]=="Hep-Mature",1],1:cutoff]))

png("Supplementary_Fig1_PCAlinScores_NewGeneListAlt.png", width=5, height=5, units="in", res=300)
par(mar=c(4,4,1,1))
barplot(barplot_dat, beside=T, col=Type.col[1:3], ylab="Avg Marker Weights")
abline(h=c(10,-10), col="grey50", lty=2) # previously 7.5, -7.5
legend("bottomright", c("Chol = PC1 + PC2", "Stem = -PC3 +PC5", "Hep = -PC2 +PC3"), fill=Type.col, bty="n")
dev.off()

# from the top 10 components take those with average weights > 10
thing <- apply(pca$rotation, 2, scale)
chol_score <- (thing[,1]+thing[,2])/2
hep_score <- (-thing[,2]+thing[,3])/2
stem_score <- (thing[,5]-thing[,3])/2

xes <- hep_score-chol_score
yes <- stem_score


SCE$Hep_score <- hep_score
SCE$Chol_score <- chol_score
SCE$Stem_score <- stem_score

png("Supplementary_Hep_Vs_Chol.png", width=5, height=5, units="in", res=300)
plot(SCE$Hep_score, SCE$Chol_score, xlab="Hepatocyte Score", ylab="Cholangiocyte Score", pch=16, col=Donor.col[factor(SCE$Donor, levels=names(Donor.col))])
dev.off()



png("Figure1_Lineage_NewGeneListAlt.png", width=5*2, height=5, units="in", res=300)
par(mfrow=c(1,2))
par(mar=c(4,4,1,1))
plot(xes, yes, xlab="Cholangiocyte ---------- Hepatocyte", ylab="Stemness", pch=16, col=Donor.col[factor(SCE$Donor, levels=names(Donor.col))])
legend("topright", names(Donor.col), col=Donor.col, pch=16, bty="n", ncol=2)
top_x <- mean(xes[SCE$Donor %in% c("HCC23", "HCC24")])
top_y <- mean(yes[SCE$Donor %in% c("HCC23", "HCC24")])
l_x <-  mean(xes[SCE$Donor %in% c("CCA5")])
l_y <-  mean(yes[SCE$Donor %in% c("CCA5")])
r_x <-  mean(xes[SCE$Donor %in% c("HCC10")])
r_y <-  mean(yes[SCE$Donor %in% c("HCC10")])
#arrows(top_x, top_y, l_x, l_y, len=0, lwd=2, col="black")
#arrows(top_x, top_y, r_x, r_y, len=0, lwd=2, col="black")


plot(xes, yes, xlab="Cholangiocyte ---------- Hepatocyte", ylab="Stemness", pch=16, col=Cycle.col[factor(SCE$CC_state, levels=names(Cycle.col))])
legend("topright", names(Cycle.col), col=Cycle.col, pch=16, bty="n")
dev.off()


# Gene Heatmap
#M3DropExpressionHeatmap(fs$Gene[1:1000], Global.norm, cell_labels=Donor)
source("~/R-Scripts/heatmap.3.R")
require("RColorBrewer")
require("gplots")

heatcolours <- rev(brewer.pal(11, "RdBu"))
col_breaks <- c(-100, seq(-2, 2, length = 10), 100)
heat_data <- as.matrix(assays(SCE)[["lognorm"]][rownames(SCE) %in% FS$Gene[1:2000],])
rownames(heat_data) <- rowData(SCE)$Symbol[rownames(SCE) %in% FS$Gene[1:2000]]
heat_data <- heat_data[rownames(heat_data) %in% Lineage[,1],]
ColColors1 <- Donor.col[factor(SCE$Donor, levels=names(Donor.col))]
ColColors2 <- Experiment.col[factor(SCE$Experiment, levels=names(Experiment.col))]
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
png("Figure1_heatmap_NewGeneListAlt.png", width=6*2, height=6, units="in", res=300)
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
Global.symbol <- rowData(SCE)$Symbol
heat_data <- as.matrix(assays(SCE)[["lognorm"]][rowSums(assays(SCE)[["lognorm"]]>0) > 50 & Global.symbol %in% Lineage[,1], ])
rownames(heat_data) <- Global.symbol[rowSums(assays(SCE)[["lognorm"]]>0) > 50 & Global.symbol %in% Lineage[,1]]
M <- complex_markers(heat_data, SCE$Donor)

lin_label <- Lineage[match(rownames(heat_data), Lineage[,1]),2]
RowColors <- c(Donor.col["CCA5"], Donor.col["CCA1"], Donor.col["HCC23"], Donor.col["HCC10"], Donor.col["HCC24"])[factor(lin_label, levels=c("Chol-Mature", "Chol-Prog", "Common-Prog", "Hep-Prog", "Hep-Mature"))]
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


heat_data2 <- as.matrix(assays(SCE)[["lognorm"]][Global.symbol %in% G_to_use,])
rownames(heat_data2) <- Global.symbol[Global.symbol %in% G_to_use]
ColColors1 <- Donor.col[factor(SCE$Donor, levels=names(Donor.col))]
ColColors2 <- Experiment.col[factor(SCE$Experiment, levels=names(Experiment.col))]
lin_label <- Lineage[match(rownames(heat_data2), Lineage[,1]),2]
RowColors2 <- c("forestgreen", "darkgoldenrod", "darkgoldenrod", "darkgoldenrod", "firebrick")[factor(lin_label, levels=c("Chol-Mature", "Chol-Prog", "Common-Prog", "Hep-Prog", "Hep-Mature"))]
RowColors2[is.na(RowColors2)] <- "white"

png("Figure1_heatmap2_NewGeneListAlt.png", width=6*2, height=6, units="in", res=300)
heatmap_output <- suppressWarnings(heatmap.2(t(heat_data2), 
            ColSideColors = RowColors2, RowSideColors = ColColors1, 
            col = heatcolours, breaks = col_breaks, scale = "row", 
            symbreaks = TRUE, trace = "none", dendrogram = "both", #labRow = gene_lab,
            key = FALSE, Rowv = TRUE, Colv = TRUE, lwid = lwid, 
            lhei = lhei, lmat = lmat, hclustfun = function(x) {
                hclust(x, method = "ward.D")
            }))
dev.off()

s <- t(heat_data2)[,heatmap_output$colInd]
c_col <- RowColors2[heatmap_output$colInd]
r_col <- ColColors1
reorder <- order(c_col);
png("Figure1_heatmap2_NewGeneListAlt_fix.png", width=6*2, height=6, units="in", res=300)
heatmap_output2 <- suppressWarnings(heatmap.2(s[,reorder],
            ColSideColors = c_col[reorder], RowSideColors = r_col,
            col = heatcolours, breaks = col_breaks, scale = "col",
            symbreaks = TRUE, trace = "none", dendrogram = "both", #labRow = gene_lab,
            key = FALSE, Rowv = TRUE, Colv = FALSE, lwid = lwid,
            lhei = lhei, lmat = lmat)
            )
dev.off()

