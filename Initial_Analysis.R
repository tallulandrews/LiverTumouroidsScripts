x <- read.table("Correct_Raw_FC_NoMulti.txt", header=T)
ann <- read.delim("Annotation_Table.txt", header=T)

type_col = c("purple", "blue", "forestgreen", "black")
plate_pch = c(16, 7, 25)

biotype <- read.table("Human_biotype.txt")
pc <- biotype[biotype[,2] == "protein_coding",1]
mt <- biotype[grep("Mt_",biotype[,2]),1]
rRNA <- biotype[biotype[,2] == "rRNA",1]
sRNAs <- biotype[biotype[,2] == "sRNA" | biotype[,2] == "snRNA" | biotype[,2] == "snoRNA" | biotype[,2] == "scaRNA",1]

gene_len = x[,1]
x <- x[,-1]
tot_reads <- colSums(x);
x <- x[rownames(x) != "__Unassigned_Various",]
exon_reads <- colSums(x);

my_legend <- function(loc) {
	legend(loc,c("CCA1", "HCC6", "HCC10", "empty", "", "Plate1","Plate2","Plate3"), col=c(type_col, "white", rep("black", times=length(plate_pch))), pch=c(rep(16, times=length(type_col)+1), plate_pch), ncol=2)
}
empty <- ann$CellID[ann$Well %in% c("P21", "P22", "P23", "P24")]
pools <- ann$CellID[ann$Well %in% c("A01", "A02", "B01", "B02", "C01", "C02")]
single_cells <- ann$CellID[ !(ann$CellID %in% pools | ann$CellID %in% empty)]

edat <- x[,match(paste("X",empty, sep=""), colnames(x))]
eplate <- ann$Plate[ann$CellID %in% empty]
eplate <- factor(eplate, levels=c("868", "869", "870"))
etype <- ann$Type[ann$CellID %in% empty]
etype <- factor(etype, levels=c("CCA1","HCC6","HCC10","no cells"))

pdat <- x[,match(paste("X", pools, sep=""), colnames(x))]
pplate <- ann$Plate[ann$CellID %in% pools]
pplate <- factor(pplate, levels=c("868", "869", "870"))
ptype <- ann$Type[ann$CellID %in% pools]
ptype[grep("HCC6", ptype)] = "HCC6"
ptype[grep("CCA1", ptype)] = "CCA1"
ptype[grep("HCC10", ptype)] = "HCC10"
ptype <- factor(ptype, levels=c("CCA1","HCC6","HCC10","no cells"))
plot(colSums(pdat), colSums(pdat > 0), col=type_col[ptype], pch=plate_pch[pplate])
my_legend("bottomright")

scdat <- x[,match(paste("X", single_cells, sep=""), colnames(x))]
scplate <- ann$Plate[ann$CellID %in% single_cells]
scplate <- factor(scplate, levels=c("868", "869", "870"))
sctype <- ann$Type[ann$CellID %in% single_cells]
sctype <- factor(sctype, levels=c("CCA1","HCC6","HCC10","no cells"))


plot(colSums(cbind(scdat, edat)), colSums(cbind(scdat,edat) > 0), col=type_col[c(sctype, etype)], pch=plate_pch[c(scplate, eplate)])
my_legend("bottomright")

my_QC <- function(dat, type, plate) {
	par(mfrow=c(2,2))
	par(mar=c(4,4,1,1))
	plot(colSums(dat), colSums(dat > 0), col=type_col[type], pch=plate_pch[plate], xlab="total counts", ylab="total genes")
	abline(h=6000)
	filter_detect <- colSums(dat > 0) > 6000

	plot(colSums(dat[rownames(dat) %in% pc,])/colSums(dat)*100, colSums(dat[rownames(dat) %in% mt,])/colSums(dat)*100, col=type_col[type], pch=plate_pch[plate], xlab="%PC", ylab="%MT")
	abline(v=70)
	filter_pc <- colSums(dat[rownames(dat) %in% pc,])/colSums(dat)*100 > 70

	plot(colSums(dat[rownames(dat) %in% pc,])/colSums(dat)*100, colSums(dat[rownames(dat) %in% rRNA,])/colSums(dat)*100, col=type_col[type], pch=plate_pch[plate], xlab="%PC", ylab="%rRNA")
	abline(v=70)
	abline(h=0.01)
	filter_rRNA <- colSums(dat[rownames(dat) %in% rRNA,])/colSums(dat)*100 < 0.01

	plot(colSums(dat[rownames(dat) %in% pc,])/colSums(dat)*100, colSums(dat[rownames(dat) %in% sRNAs,])/colSums(dat)*100, col=type_col[type], pch=plate_pch[plate], xlab="%PC", ylab="%sRNA")
	abline(v=70)
	abline(h=0.04)
	filter_sRNA <- colSums(dat[rownames(dat) %in% sRNAs,])/colSums(dat)*100 < 0.04
	return(cbind(filter_detect, filter_pc, filter_rRNA, filter_sRNA))
}

png("SC_QC.png", width=7, height=7, units="in", res=300)
sc_qc <- my_QC(cbind(scdat, edat), c(sctype, etype), c(scplate, eplate))
dev.off()
png("Pool_QC.png", width=7, height=7, units="in", res=300)
pool_qc <- my_QC(pdat, ptype, pplate)
dev.off()

# Do QC
sc_good <- cbind(scdat, edat); sc_good <- sc_good[,rowSums(sc_qc) == 4];
sc_plate <-  c(scplate, eplate); sc_plate <- sc_plate[rowSums(sc_qc) == 4];
sc_type <- c(as.character(sctype), as.character(etype)); sc_type<- sc_type[rowSums(sc_qc) == 4]; sc_type <- factor(sc_type,  levels=c("CCA1","HCC6","HCC10"))


# PCA
sc_good <- sc_good[rowSums(sc_good > 5) > 2,]
size_factors <- colSums(sc_good)
norm <- sc_good/size_factors*median(size_factors)

PCA <- prcomp(log(norm+1)/log(2))
png("SC_QCed_lognorm_PCA.png", width=6, height=6, units="in", res=300)
plot(PCA$rotation[,1], PCA$rotation[,2], col=type_col[sc_type], pch=plate_pch[sc_plate], xlab=paste("PC1 (",round(PCA$sdev[1]/sum(PCA$sdev)*100, digits=2)," %)",sep=""), ylab=paste("PC2 (",round(PCA$sdev[2]/sum(PCA$sdev)*100, digits=2)," %)",sep=""))
dev.off()

# Feature Selection
require("M3Drop")

png("M3Drop_FS.png", width=6, height=6, units="in", res=300)
M3D_FS <- M3DropFeatureSelection(norm, mt_method="bonferroni", mt_threshold=0.05)
dev.off()

png("M3Drop_heatmap.png", width=8, height=8, units="in", res=300)
heatout <- M3DropExpressionHeatmap(M3D_FS, norm, cell_labels=sc_type)
dev.off()

# Digging out Substructure
cell_groups <- M3DropGetHeatmapClusters(heatout, k=4, type="cell") # type="gene" is not working!
markers <- M3DropGetMarkers(norm, cell_groups)
markers2 <- M3DropGetMarkers(norm[,cell_groups==3 | cell_groups==4], cell_groups[cell_groups==3 | cell_groups==4])


# DE 
#source("/nfs/users/nfs_t/ta6/NetworkInferencePipeline/Dropouts/Results_Git/DiffExpr/DE_functions.R")

# Monocle DE
require("monocle")
cell_data <- data.frame(group=sc_type, batch=sc_plate)
rownames(cell_data) = colnames(sc_good)
pd <- new("AnnotatedDataFrame", cell_data)
Obj <- newCellDataSet(as.matrix(sc_good), phenoData=pd, expressionFamily=negbinomial.size())
Obj <- estimateSizeFactors(Obj)
Obj <- estimateDispersions(Obj)
diff_test_res <- differentialGeneTest(Obj,fullModelFormulaStr="~group+batch")
## Need to edit the below
diff_test_res[diff_test_res[,1] != "OK",4] = 1;
diff_test_res[diff_test_res[,1] != "OK",5] = 1;
pvals <- diff_test_res[,4]
qvals <- diff_test_res[,5]

norm <- t(t(exprs(Obj))/pData(Obj)$Size_Factor)
Means = aggregate(t(norm), by=list(pData(Obj)$group),mean);
rownames(Means) = Means[,1]; Means = Means[,-1];
TABLE <- cbind(rownames(norm),t(Means),pvals,p.adjust(pvals, method="fdr"));
colnames(TABLE) <- c("Gene",rownames(Means),"p.value","q.value");




# Clustering & Markers
source("/nfs/users/nfs_t/ta6/R-Scripts/heatboxplot.R")
