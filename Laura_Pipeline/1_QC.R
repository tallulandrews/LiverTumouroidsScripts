
args <- commandArgs(trailingOnly=TRUE)

require("RColorBrewer")
# countfile, annotationfile, prefix
#args <- c("Correct_Raw_FC_NoMulti.txt","Annotation_Table.txt","Test")


source("/nfs/users/nfs_t/ta6/Collaborations/LiverOrganoids/Laura_Pipeline/0_ColourScheme.R")
#plate_pch = c(16, 7, 25, 3, 4, 23, 1)
#type_col = brewer.pal(8, "Set1")


dat_matrix <- read.table(args[1], header=T)
ann <- read.delim(args[2], header=T)
ann$CellID <- paste("X", ann$CellID, sep="")

ann$Plate <- factor(ann$Plate)
plate_pch <- plate_pch[1:length(levels(ann$Plate))]

#type_col = c("purple", "blue", "forestgreen", "black")
#plate_pch = c(16, 7, 25)


tot_reads <- colSums(dat_matrix);
dat_matrix <- dat_matrix[rownames(dat_matrix) != "__Unassigned_Various",]
gene_len = dat_matrix[,1]
dat_matrix <- dat_matrix[,-1]
tot_reads <- tot_reads[-1]

empty <- ann$CellID[grep("no cells", ann$Type) ]
pools <- ann$CellID[ grep("cells", ann$Type) ]
pools <- pools[! (pools %in% empty) ]
single_cells <- ann$CellID[ !(ann$CellID %in% pools | ann$CellID %in% empty)]

ann$Type <- unlist(lapply(strsplit(as.character(ann$Type), " "), function(dat_matrix){dat_matrix[length(dat_matrix)]}))
ann$Type[ann$Type == "cells"] <- "Empty"
lvls <- sort(unique(ann$Type)); lvls <- lvls[lvls != "Empty"]; 
ann$Type <- factor(ann$Type, levels=c(lvls, "Empty"))

type_col <- c(type_col[1:(length(levels(ann$Type))-1)], "black")

# Error Checking
dat_matrix <- dat_matrix[, match(ann$CellID, colnames(dat_matrix))]
print(paste("Cell order correct?", identical(ann$CellID, colnames(dat_matrix))))


edat <- dat_matrix[,match(empty, colnames(dat_matrix))]
eann <- ann[ann$CellID %in% empty,]
#eplate <- ann$Plate[ann$CellID %in% empty]
#eplate <- factor(eplate, levels=c("868", "869", "870"))
#etype <- ann$Type[ann$CellID %in% empty]
#etype <- factor(etype, levels=c("CCA1","HCC6","HCC10","no cells"))

pdat <- dat_matrix[,match(pools, colnames(dat_matrix))]
pann <- ann[ann$CellID %in% pools,]
#pplate <- ann$Plate[ann$CellID %in% pools]
#pplate <- factor(pplate, levels=c("868", "869", "870"))
#ptype <- ann$Type[ann$CellID %in% pools]
#ptype[grep("HCC6", ptype)] = "HCC6"
#ptype[grep("CCA1", ptype)] = "CCA1"
#ptype[grep("HCC10", ptype)] = "HCC10"
#ptype <- factor(ptype, levels=c("CCA1","HCC6","HCC10","no cells"))

scdat <- dat_matrix[,match(single_cells, colnames(dat_matrix))]
scann <- ann[ann$CellID %in% single_cells,]
#scplate <- ann$Plate[ann$CellID %in% single_cells]
#scplate <- factor(scplate, levels=c("868", "869", "870"))
#sctype <- ann$Type[ann$CellID %in% single_cells]
#sctype <- factor(sctype, levels=c("CCA1","HCC6","HCC10","no cells"))

my_legend <- function(loc) {
	legend(loc,c(levels(ann$Type), "", levels(ann$Plate)), 
		col=c(type_col, "white", rep("black", times=length(plate_pch))), 
		pch=c(rep(16, times=length(type_col)+1), plate_pch), ncol=2, bty="n")
}

png(paste(args[3], "All_genes_vs_counts.png", sep="_"), width=7, height=7, units="in", res=300)
plot(colSums(cbind(scdat, edat)), colSums(cbind(scdat,edat) > 0), col=type_col[c(scann$Type, eann$Type)], pch=plate_pch[c(scann$Plate, eann$Plate)], xlab="Total Counts", ylab="Total Genes")
points(colSums(pdat), colSums(pdat > 0), col=type_col[pann$Type], pch=plate_pch[pann$Plate], cex=3)
my_legend("bottomright")
dev.off()

biotype <- read.table("/lustre/scratch117/cellgen/team218/TA/LiverOrganoids/Human_biotype.txt")

require("scater")
source("~/R-Scripts/Ensembl_Stuff.R")
g_d <- data.frame(Length=gene_len, Symbol=map_symbol_ensg(rownames(dat_matrix), is.org="Hsap", is.name="ensg"), biotype=biotype[match(rownames(scdat),biotype[,1]), 2]);
rownames(g_d) = rownames(dat_matrix);
gd <- new("AnnotatedDataFrame", data=g_d)
# Pools
pd <- new("AnnotatedDataFrame", data=data.frame(pann))
rownames(pd) <- colnames(pdat)
pSCE <- newSCESet(countData = pdat, phenoData = pd, featureData = gd)
pSCE <- calculateQCMetrics(pSCE)
# single-cells
pd <- new("AnnotatedDataFrame", data=data.frame(scann))
rownames(pd) <- colnames(scdat)
scSCE <- newSCESet(countData = scdat, phenoData = pd, featureData = gd)
scSCE <- calculateQCMetrics(scSCE)
# empty
pd <- new("AnnotatedDataFrame", data=data.frame(eann))
rownames(pd) <- colnames(edat)
eSCE <- newSCESet(countData = edat, phenoData = pd, featureData = gd)
eSCE <- calculateQCMetrics(eSCE)



# Thresholds

calculate_thresholds <- function(scSCE, pSCE, eSCE, biotype) {
	pc <- biotype[biotype[,2] == "protein_coding",1]
	mt <- biotype[grep("Mt_",biotype[,2]),1]
	rRNA <- biotype[biotype[,2] == "rRNA",1]
	sRNAs <- biotype[biotype[,2] == "sRNA" | biotype[,2] == "snRNA" | biotype[,2] == "snoRNA" | biotype[,2] == "scaRNA",1]


	Pos_ctrl <- pSCE$total_features;
	Pos_ctrl <- c(min(Pos_ctrl), max(Pos_ctrl))
	if (median(Pos_ctrl) < 500 | min(Pos_ctrl) < 100) { # Pos_ctrl failed
		#Pos_ctrl <- Pos_ctrl[Pos_ctrl > 1000];
		#Pos_ctrl <- c(min(Pos_ctrl), max(scSCE$total_features))
		require("mclust")
		f <- Mclust(cbind(log(scSCE$total_counts+1), log(scSCE$total_features+1)), G=2)
		good <- f$classification == 2
		Pos_ctrl <- c(min(scSCE$total_features[good]), max(scSCE$total_features[good]))
	}
	scmean <- mean(scSCE$total_features[scSCE$total_features > Pos_ctrl[1] & scSCE$total_features < Pos_ctrl[2]])
	scsd <- sd(scSCE$total_features[scSCE$total_features > Pos_ctrl[1] & scSCE$total_features < Pos_ctrl[2]])
	#gLo <- round(max(min(Pos_ctrl, scmean-2.5*scsd), quantile(eSCE$total_features, probs=0.75)), digits=-3) # changed from 3 to 2.5 on 19 April 2018
	gLo <- round(min(Pos_ctrl, scmean-2.5*scsd), digits=-3) # changed from 3 to 2.5 on 19 April 2018
	#if (max(eSCE$total_features) > gLo) {
	#	gLo <- round(max(eSCE$total_features), digits=-3)
	#}
	gHi <- round(max(Pos_ctrl, scmean+3*scsd), digits=-3)


	pc_prop <- colSums(counts(scSCE)[rownames(scSCE) %in% pc,])/scSCE$total_counts*100
	m <- mean(pc_prop, na.rm=T); s <- sd(pc_prop, na.rm=T)+0.001;
	pcLo <- signif(m-3*s, digits=1)

	mt_prop <- colSums(counts(scSCE)[rownames(scSCE) %in% mt,])/scSCE$total_counts*100
	m <- mean(mt_prop, na.rm=T); s <- sd(mt_prop, na.rm=T)+0.001;
	mtHi <- signif(m+3*s, digits=1)
	r_prop <- colSums(counts(scSCE)[rownames(scSCE) %in% rRNA,])/scSCE$total_counts*100
	m <- mean(r_prop, na.rm=T); s <- sd(r_prop, na.rm=T)+0.001;
	rHi <- signif(m+3*s, digits=1)
	s_prop <- colSums(counts(scSCE)[rownames(scSCE) %in% sRNAs,])/scSCE$total_counts*100
	m <- mean(s_prop, na.rm=T); s <- sd(s_prop, na.rm=T)+0.001;
	sHi <- signif(m+3*s, digits=1)

	thresholds <- list(gHi = gHi, gLo=gLo, pc=pcLo, mt=mtHi, rRNA=rHi, sRNA=sHi)
	return(thresholds)
}

thresholds <- calculate_thresholds(scSCE, pSCE, eSCE, biotype)

my_QC <- function(SCE, thresholds) {
	pc <- rownames(fData(SCE))[fData(SCE)$biotype == "protein_coding"]
	mt <- rownames(fData(SCE))[fData(SCE)$biotype == "Mt_"]
	rRNA <- rownames(fData(SCE))[fData(SCE)$biotype == "rRNA"]
	sRNAs <- rownames(fData(SCE))[fData(SCE)$biotype %in% c("sRNA","snRNA", "snoRNA", "scaRNA")]

	par(mfrow=c(2,2))
	par(mar=c(4,4,1,1))
	dat <- counts(SCE)
	plot(SCE$total_counts, SCE$total_features, col=type_col[SCE$Type], pch=plate_pch[SCE$Plate], xlab="total counts", ylab="total genes")

	abline(h=c(thresholds["gHi"], thresholds["gLo"]))
	filter_detect <- SCE$total_features > thresholds["gLo"] & SCE$total_features < thresholds["gHi"]

	plot(colSums(dat[rownames(dat) %in% pc,])/colSums(dat)*100, colSums(dat[rownames(dat) %in% mt,])/colSums(dat)*100, col=type_col[SCE$Type], pch=plate_pch[SCE$Plate], xlab="%PC", ylab="%MT")
	abline(v=thresholds["pc"])
	filter_pc <- colSums(dat[rownames(dat) %in% pc,])/colSums(dat)*100 > thresholds["pc"]
	abline(h=thresholds["mt"])
	filter_mt <- colSums(dat[rownames(dat) %in% mt,])/colSums(dat)*100 < thresholds["mt"]

	plot(colSums(dat[rownames(dat) %in% pc,])/colSums(dat)*100, colSums(dat[rownames(dat) %in% rRNA,])/colSums(dat)*100, col=type_col[SCE$Type], pch=plate_pch[SCE$Plate], xlab="%PC", ylab="%rRNA")
	abline(v=thresholds["pc"])
	abline(h=thresholds["rRNA"])
	filter_rRNA <- colSums(dat[rownames(dat) %in% rRNA,])/colSums(dat)*100 < thresholds["rRNA"]

	plot(colSums(dat[rownames(dat) %in% pc,])/colSums(dat)*100, colSums(dat[rownames(dat) %in% sRNAs,])/colSums(dat)*100, col=type_col[SCE$Type], pch=plate_pch[SCE$Plate], xlab="%PC", ylab="%sRNA")
	abline(v=thresholds["pc"])
	abline(h=thresholds["sRNA"])
	filter_sRNA <- colSums(dat[rownames(dat) %in% sRNAs,])/colSums(dat)*100 < thresholds["sRNA"]
	pData(SCE)$filter_pc <- filter_pc
	pData(SCE)$filter_detect <- filter_detect
	pData(SCE)$filter_mtRNA <- filter_mt
	pData(SCE)$filter_rRNA <- filter_rRNA
	pData(SCE)$filter_sRNA <- filter_sRNA
	pData(SCE)$Good <- filter_pc & filter_detect & filter_mt & filter_rRNA & filter_sRNA
	return(SCE)
}


png(paste(args[3],"SC_QC.png", sep="_"), width=7, height=7, units="in", res=300)
sc_qc <- my_QC(scSCE, thresholds)
dev.off()
png(paste(args[3],"Pool_QC.png", sep="_"), width=7, height=7, units="in", res=300)
pool_qc <- my_QC(pSCE, thresholds)
dev.off()
empty_qc <- my_QC(eSCE, thresholds)

sum(pData(sc_qc)$Good)
sum(pData(pool_qc)$Good)
sum(pData(empty_qc)$Good)

sink(paste(args[3], "QC_Thresholds.txt", sep="_"))
print(thresholds)
sink()

# Plate position

for(p in levels(scSCE$Plate)) {
	file=paste(args[3], "SC",p,"PlateQC.png", sep="_")
	print(file)
	png(file, width=8*2, height=8, units="in", res=300)
	scater::plotPlatePosition(scSCE[,scSCE$Plate==p], plate_position=scSCE$Well[scSCE$Plate==p], colour_by="total_features")
	dev.off()
	
	file=paste(args[3], "SC",p,"PlateFiltered.png", sep="_")
	png(file, width=8*2, height=8, units="in", res=300)
        scater::plotPlatePosition(sc_qc[,sc_qc$Plate==p], plate_position=sc_qc$Well[sc_qc$Plate==p], colour_by="Good")
        dev.off()
}

table(pData(sc_qc)$Good, sc_qc$Plate)
table(pData(sc_qc)$Good, sc_qc$Type)

if (sum(sc_qc$Plate == 3313) > 0) {
sc_qc$Good[sc_qc$Plate == 3313] <- FALSE
}



# Normalization
scale_factor <- median(pData(sc_qc)$total_counts[pData(sc_qc)$Good])
set_exprs(sc_qc, "norm") <- t(t(counts(sc_qc))/pData(sc_qc)$total_counts*scale_factor)
set_exprs(sc_qc, "lognorm") <- log(get_exprs(sc_qc, "norm")+1)/log(2)
saveRDS(sc_qc, file=paste(args[3], "Initial.rds", sep="_"))


scQC <- sc_qc[,pData(sc_qc)$Good]

#require("scran"); require("limSolve");
#qclust <- quickCluster(scQC, min.size = 50)
#scQC <- computeSumFactors(scQC, sizes=c(10,15,20,25), clusters = qclust, positive=TRUE)
#scQC <- normalize(scQC) # get Inf?
#
#set_exprs(scQC, "scran") <- exprs(scQC)
#exprs(scQC) <- get_exprs(scQC, "lognorm")
#
#plot(pData(scQC)$total_counts/scale_factor, pData(scQC)$size_factor, xlab="cpm", ylab="scran");
#abline(b=1, a=0, col="red")

# Batch Effect Removal
factor_counts <- function (vec) {
    vec <- factor(vec)
    x <- split(seq(length(vec)), vec)
    result <- sapply(x, function(a) length(a))
    return(result)
}


glm_fun <- function(x) {
	res <- glm(x~pData(scQC)$Plate+pData(scQC)$Type)
	ncells <- factor_counts(pData(scQC)$Plate)
	n_nozero <- sapply(split(seq(length(x)), pData(scQC)$Plate), function(a) sum(x[a] > 0) )
	effect <- c(0,res$coef[2:length(levels(pData(scQC)$Plate))]) * ncells/n_nozero
	effect[!is.finite(effect) | is.na(effect)] <- 0;

	correction_factor <- rep(effect, times=ncells)
	corrected <- x-correction_factor
	corrected[x == 0] <- 0
	corrected[corrected < 0] <- 0
	return(corrected)
}

batch_corrected <- apply(get_exprs(scQC, "lognorm"), 1, glm_fun)
batch_corrected[is.na(batch_corrected)] <- 0
batch_corrected[!is.finite(batch_corrected)] <- 0
set_exprs(scQC, "batchcor") <- t(batch_corrected)

saveRDS(scQC, file=paste(args[3], "QCed_BatchCor.rds", sep="_"))


PCA_raw <- prcomp(get_exprs(scQC, "lognorm"))
PCA_cor <- prcomp(get_exprs(scQC, "batchcor"))
png(paste(args[3],"QCed_Norm_BatchCor_PCA.png", sep="_"), width=6*2, height=6, units="in", res=300)
par(mfrow=c(1,2))
plot(PCA_raw$rotation[,1], PCA_raw$rotation[,2], col=type_col[pData(scQC)$Type], pch=plate_pch[scQC$Plate], xlab=paste("PC1 (",round(PCA_raw$sdev[1]/sum(PCA_raw$sdev)*100, digits=2)," %)",sep=""), ylab=paste("PC2 (",round(PCA_raw$sdev[2]/sum(PCA_raw$sdev)*100, digits=2)," %)",sep=""), main="logFPKM")
plot(PCA_cor$rotation[,1], PCA_cor$rotation[,2], col=type_col[pData(scQC)$Type], pch=plate_pch[scQC$Plate], xlab=paste("PC1 (",round(PCA_cor$sdev[1]/sum(PCA_cor$sdev)*100, digits=2)," %)",sep=""), ylab=paste("PC2 (",round(PCA_cor$sdev[2]/sum(PCA_cor$sdev)*100, digits=2)," %)",sep=""), main="Batch Corrected")
dev.off()

types <- sort(unique(scQC$Type))
for(ty in types) {
	subset <- scQC[,pData(scQC)$Type==as.character(ty)]
	subset <- calculateQCMetrics(subset);	

	# Feature Selection
	require("M3Drop")
	# 13 Feb 2018 modification
	FS <- M3DropFeatureSelection(get_exprs(subset,"norm"), mt_method="fdr", mt_threshold=1)
	fData(subset)$M3Drop_q.value <- FS[match(rownames(subset), FS$Gene),]$q.value
	HVG <- BrenneckeGetVariableGenes(get_exprs(subset,"norm"), fdr=2)
	fData(subset)$HVG_q.value <- HVG[match(rownames(subset), HVG$Gene),]$q.value
	#fData(subset)$is.Feature <- rownames(fData(subset)) %in% FS$Gene

	mat <- get_exprs(subset, "lognorm")[fData(subset)$M3Drop_q.value < 0.05,]

	# Dim Reductions
	set.seed(1091)

	#PCA <- prcomp(mat)
	#pData(subset)$PC1 <- PCA$rotation[,1] 
	#pData(subset)$PC2 <- PCA$rotation[,2]

	#require("Rtsne")
	#set.seed(123)
	#tsne <- Rtsne(t(mat), perplexity=10)
	#pData(subset)$Tsne1 <- tsne$Y[,1]
	#pData(subset)$Tsne2 <- tsne$Y[,2]

	#require("destiny")
	#set.seed(1)
	#dm <- DiffusionMap(t(mat))
	#pData(subset)$DM1 <- eigenvectors(dm)[,1]
	#pData(subset)$DM2 <- eigenvectors(dm)[,2]
	#pData(subset)$DM3 <- eigenvectors(dm)[,3]
	saveRDS(subset, file=paste(as.character(ty), args[3], "QCed_BatchCor_DimRed.rds", sep="_"))
}
