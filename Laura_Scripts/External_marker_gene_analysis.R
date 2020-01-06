require("scater")
require("RColorBrewer")


CCA1 <- readRDS("/lustre/scratch117/cellgen/team218/TA/LiverOrganoids/CCA1_SC3.rds")
HCC6 <- readRDS("/lustre/scratch117/cellgen/team218/TA/LiverOrganoids/HCC6_SC3.rds")
HCC10 <- readRDS("/lustre/scratch117/cellgen/team218/TA/LiverOrganoids/HCC10_SC3.rds")

source("/nfs/users/nfs_t/ta6/R-Scripts/Ensembl_Stuff.R")

fData(CCA1)$feature_symbol <- map_symbol_ensg(rownames(CCA1), is.org="Hsap", is.name="ensg")
fData(HCC6)$feature_symbol <- map_symbol_ensg(rownames(HCC6), is.org="Hsap", is.name="ensg")
fData(HCC10)$feature_symbol <- map_symbol_ensg(rownames(HCC10), is.org="Hsap", is.name="ensg")

require("CellTypeProfiles")

CCA1_markers <- complex_markers(exprs(CCA1), CCA1$sc3_6_clusters)
fData(CCA1)$is.marker <- rownames(fData(CCA1)) %in% rownames(CCA1_markers[CCA1_markers$q.value<0.05,])
HCC6_markers <- complex_markers(exprs(HCC6), HCC6$sc3_6_clusters)
fData(HCC6)$is.marker <- rownames(fData(HCC6)) %in% rownames(HCC6_markers[HCC6_markers$q.value<0.05,])
HCC10_markers <- complex_markers(exprs(HCC10), HCC10$sc3_6_clusters)
fData(HCC10)$is.marker <- rownames(fData(HCC10)) %in% rownames(HCC10_markers[HCC10_markers$q.value<0.05,])

get_dm_coords <- function(SCE) {
        require("destiny")
        set.seed(1)
        size_factors <- colSums(counts(SCE))
        norm <- t(t(counts(SCE))/size_factors*median(size_factors))
        norm <- log(norm+1)/log(2)
	norm <- norm[fData(SCE)$is.marker,]
#       norm <- exprs(SCE)

        plate <- factor(pData(SCE)[,2], levels=c("868", "869", "870"))
#       type  <- factor(pData(SCE)[,4], levels=c("CCA1","HCC6","HCC10","no cells"))
        type  <- pData(SCE)$sc3_6_clusters

        dm <- DiffusionMap(t(norm))
        dm_dims <- eigenvectors(dm)
#        plot(dm_dims[,1], dm_dims[,2], col=group_cols[type], pch=plate_pch[plate], xlab="Dimension 1", ylab="Dimension 2", main="")
#        legend("topright",as.character(1:6),fill=group_cols, bty="n")
	pData(SCE)$dm1 <- dm_dims[,1]
	pData(SCE)$dm2 <- dm_dims[,2]
        return(SCE)
}

CCA1 <- get_dm_coords(CCA1)
HCC6 <- get_dm_coords(HCC6)
HCC10 <- get_dm_coords(HCC10)


### Load External Markers ###
SC <- read.table("/nfs/users/nfs_t/ta6/Collaborations/LiverOrganoids/FromLaura_Camp_Hepatoblast_Markers.txt")
MSC <- read.table("/nfs/users/nfs_t/ta6/Collaborations/LiverOrganoids/FromLaura_Camp_Mesenchymal_Markers.txt")
HC <- read.table("/nfs/users/nfs_t/ta6/Collaborations/LiverOrganoids/FromLaura_Camp_Hepatocyte_Markers.txt")
Chol_arrays <- readRDS("/lustre/scratch117/cellgen/team218/TA/LiverOrganoids/ExternalData/Sampaziotis_CholMarkers.rds")
HPA_Chol <- readRDS("/lustre/scratch117/cellgen/team218/TA/OtherDownloadedData/HPA_Markers_liver_bile_duct_cells.rds")
HPA_Chol <- HPA_Chol[HPA_Chol$Level %in% c("Medium", "High") & HPA_Chol$Reliability %in% c("Supported", "Approved") & HPA_Chol$Proportion > 0.04,]
HPA_Hep <- readRDS("/lustre/scratch117/cellgen/team218/TA/OtherDownloadedData/HPA_Markers_liver_hepatocytes.rds")
HPA_Hep <- HPA_Hep[HPA_Hep$Level %in% c("Medium", "High") & HPA_Hep$Reliability %in% c("Supported", "Approved") & HPA_Hep$Proportion > 0.04,]


Extern_table <- data.frame(Gene=c(as.character(SC[,1]), as.character(MSC[,1]), as.character(HC[,1]), as.character(HPA_Chol$Gene.name), as.character(HPA_Hep$Gene.name), names(Chol_arrays[[1]]), names(Chol_arrays[[2]]), names(Chol_arrays[[3]]) ), Type=rep(c("SC", "MSC", "Hep", "Chol_HPA", "Hep_HPA", names(Chol_arrays)), times=c(nrow(SC), nrow(MSC), nrow(HC), nrow(HPA_Chol), nrow(HPA_Hep), length(Chol_arrays[[1]]), length(Chol_arrays[[2]]), length(Chol_arrays[[3]]))))
#CC <- ?

check_assignment <- function(marker_obj, ext_markers, SCE) {
	assign <- marker_obj[,-c(1,ncol(marker_obj), ncol(marker_obj)-1)]
	assign <- assign[marker_obj$q.value < 0.05,]
	shared <- colSums(assign[rownames(assign) %in% rownames(fData(SCE))[fData(SCE)$feature_symbol %in% ext_markers],])
	tot <- colSums(assign)
#	shared/tot
	return(shared)
}


N_types <- length(unique(Extern_table[,2]))
Types <- unique(Extern_table[,2])

require("CellTypeProfiles")
# Boxplots
png("CCA1_external_markers.png", width=2.75*4, height=2*2, units="in", res=300)
par(mfrow=c(2,4))
par(mar=c(2,2,1,1))
for(t in Types) {
	thing <- max(my_row_mean_aggregate(exprs(CCA1)[fData(CCA1)$feature_symbol %in% Extern_table[,1],], CCA1$sc3_6_clusters));
	boxplot(my_row_mean_aggregate(exprs(CCA1)[fData(CCA1)$feature_symbol %in% Extern_table[Extern_table[,2]==t,1],], CCA1$sc3_6_clusters), 
		notch=TRUE, col=brewer.pal(6,"Set3"), main=paste("CCA1 - ", t), ylim=c(0,thing))
}
dev.off()

png("HCC6_external_markers.png", width=2.75*4, height=2*2, units="in", res=300)
par(mfrow=c(2,4))
par(mar=c(2,2,1,1))
for(t in Types) {
	thing <- max(my_row_mean_aggregate(exprs(HCC6)[fData(HCC6)$feature_symbol %in% Extern_table[,1],], HCC6$sc3_6_clusters));
	boxplot(my_row_mean_aggregate(exprs(HCC6)[fData(HCC6)$feature_symbol %in% Extern_table[Extern_table[,2]==t,1],], HCC6$sc3_6_clusters), 
		notch=TRUE, col=brewer.pal(6,"Set3"), main=paste("HCC6 - ", t), ylim=c(0,thing))
}
dev.off()

png("HCC10_external_markers.png", width=2.75*4, height=2*2, units="in", res=300)
par(mfrow=c(2,4))
par(mar=c(2,2,1,1))
for(t in Types) {
	thing <- max(my_row_mean_aggregate(exprs(HCC10)[fData(HCC10)$feature_symbol %in% Extern_table[,1],], HCC10$sc3_6_clusters));
	boxplot(my_row_mean_aggregate(exprs(HCC10)[fData(HCC10)$feature_symbol %in% Extern_table[Extern_table[,2]==t,1],], HCC10$sc3_6_clusters), 
		notch=TRUE, col=brewer.pal(6,"Set3"), main=paste("HCC10 - ", t), ylim=c(0,thing))
}
dev.off()

# Markerplots
overall_score_plot <- function(SCE, ext_markers) {
	is.marker <- fData(SCE)$feature_symbol %in% as.character(ext_markers);
	score <- colMeans(exprs(SCE)[is.marker,])

	n = 6
	colours = brewer.pal(n, "Blues");
	bins <- cut(score, breaks=quantile(score, probs=seq(from=0, to=1, length=n)))
	plot(SCE$dm1, SCE$dm2, col=colours[bins]);
	return(score)
}

intensity_plot <- function(SCE, gene) {
	vals <- exprs(SCE)[fData(SCE)$feature_symbol == gene,];
		
	splits <- quantile(vals, probs=seq(from=0, to=1, length=8))
	splits <- unique(splits);
	splits <- seq(from=0, to=max(splits), length=5)
	bins <- cut(vals, breaks=splits, include.lowest=TRUE);
	my_cols <- colorRampPalette(brewer.pal(6, "Blues"))(length(levels(bins)));
	plot(SCE$dm1, SCE$dm2, col=my_cols[bins], pch=16, main=gene, xaxt="n", yaxt="n", xlab="", ylab="", bty="n")
}
Common_markers <- function(SCE, SCE_markers, ext_markers, file=NA) {
	common_markers <- SCE_markers[fData(SCE)$feature_symbol %in% ext_markers,]
	common_markers <- common_markers[order(-common_markers$AUC),]
	sig <- common_markers[common_markers$q.value < 0.05 & common_markers$AUC > 0.7,]
	assign <- sig[,c(-1, -ncol(sig), -ncol(sig)+1)]
	sig <- sig[rowSums(assign) < ncol(assign)/2,]
	sig$name <- fData(SCE[match(rownames(sig), rownames(SCE)),])$feature_symbol
	sig <- sig[order(-sig$AUC),]
	max_panels <- 15;
	if (!is.na(file)) {
		png(file, width=1*min(nrow(sig), max_panels)/2, height=4, units="in", res=300)
	}

	par(cex=0.95)
	par(mfrow=c(3, min(ceiling(nrow(sig)/3), max_panels/3)))
	par(mar=c(0,0,2,0))
	if (nrow(sig) > max_panels) {
		for(g in sig$name[1:max_panels]) {intensity_plot(SCE, g)}
	} else {
		for(g in sig$name) {intensity_plot(SCE, g)}
	}
	if (!is.na(file)) {
		dev.off()
	}
	return(sig)
}

#for(t in Types) {
#	out <- Common_markers(CCA1, CCA1_markers, Extern_table[Extern_table[,2]==t,1])
#	print(colSums(out[,c(-1, -ncol(out), -ncol(out)+1)]))
#}

out <- Common_markers(CCA1, CCA1_markers, Extern_table[Extern_table[,2]=="SC",1], "Extern_SC_in_CCA1.png")
CCA1_sc_agree <- rownames(out)[out[,5]==1 | out[,2]==1 ]
nrow(out); colSums(out[,2:7])
out <- Common_markers(CCA1, CCA1_markers, Extern_table[Extern_table[,2]=="Hep",1], "Extern_Hep_in_CCA1.png")
nrow(out); colSums(out[,2:7])
out <- Common_markers(CCA1, CCA1_markers, Extern_table[Extern_table[,2]=="MSC",1], "Extern_MSC_in_CCA1.png")
nrow(out); colSums(out[,2:7])
out <- Common_markers(CCA1, CCA1_markers, Extern_table[Extern_table[,2]=="Chol",1], "Extern_Chol_in_CCA1.png")
nrow(out); colSums(out[,2:7])
out <- Common_markers(CCA1, CCA1_markers, Extern_table[Extern_table[,2]=="CLC",1], "Extern_CLC_in_CCA1.png")
nrow(out); colSums(out[,2:7])
out <- Common_markers(CCA1, CCA1_markers, Extern_table[Extern_table[,2]=="Mature",1], "Extern_C-Mature_in_CCA1.png")
nrow(out); colSums(out[,2:7])
out <- Common_markers(CCA1, CCA1_markers, Extern_table[Extern_table[,2]=="Chol_HPA",1], "Extern_Chol-HPA_in_CCA1.png")
nrow(out); colSums(out[,2:7])
out <- Common_markers(CCA1, CCA1_markers, Extern_table[Extern_table[,2]=="Hep_HPA",1], "Extern_Hep-HPA_in_CCA1.png")
nrow(out); colSums(out[,2:7])

out <- Common_markers(HCC6, HCC6_markers, Extern_table[Extern_table[,2]=="SC",1], "Extern_SC_in_HCC6.png")
HCC6_sc_agree <- rownames(out)[out[,5]==1 ]
nrow(out); colSums(out[,2:7])
out <- Common_markers(HCC6, HCC6_markers, Extern_table[Extern_table[,2]=="Hep",1], "Extern_Hep_in_HCC6.png")
HCC6_hep_agree <- rownames(out5)[out[,3]==1 | out[,2]==1 ]
nrow(out); colSums(out[,2:7])
out <- Common_markers(HCC6, HCC6_markers, Extern_table[Extern_table[,2]=="MSC",1], "Extern_MSC_in_HCC6.png")
nrow(out); colSums(out[,2:7])
out <- Common_markers(HCC6, HCC6_markers, Extern_table[Extern_table[,2]=="Chol",1], "Extern_Chol_in_HCC6.png")
nrow(out); colSums(out[,2:7])
out <- Common_markers(HCC6, HCC6_markers, Extern_table[Extern_table[,2]=="CLC",1], "Extern_CLC_in_HCC6.png")
nrow(out); colSums(out[,2:7])
out <- Common_markers(HCC6, HCC6_markers, Extern_table[Extern_table[,2]=="Mature",1], "Extern_C-Mature_in_HCC6.png")
nrow(out); colSums(out[,2:7])
out <- Common_markers(HCC6, HCC6_markers, Extern_table[Extern_table[,2]=="Chol_HPA",1], "Extern_Chol-HPA_in_HCC6.png")
nrow(out); colSums(out[,2:7])
out <- Common_markers(HCC6, HCC6_markers, Extern_table[Extern_table[,2]=="Hep_HPA",1], "Extern_Hep-HPA_in_HCC6.png")
nrow(out); colSums(out[,2:7])

out <- Common_markers(HCC10, HCC10_markers, Extern_table[Extern_table[,2]=="SC",1], "Extern_SC_in_HCC10.png")
HCC10_sc_agree <- rownames(out)[out[,2]==1 ]
nrow(out); colSums(out[,2:7])
out <- Common_markers(HCC10, HCC10_markers, Extern_table[Extern_table[,2]=="Hep",1], "Extern_Hep_in_HCC10.png")
HCC10_hep_agree <- rownames(out)[out[,6]==1 | out[,7]==1 ]
nrow(out); colSums(out[,2:7])
out <- Common_markers(HCC10, HCC10_markers, Extern_table[Extern_table[,2]=="MSC",1], "Extern_MSC_in_HCC10.png")
nrow(out); colSums(out[,2:7])
out <- Common_markers(HCC10, HCC10_markers, Extern_table[Extern_table[,2]=="Chol",1], "Extern_Chol_in_HCC10.png")
nrow(out); colSums(out[,2:7])
out <- Common_markers(HCC10, HCC10_markers, Extern_table[Extern_table[,2]=="CLC",1], "Extern_CLC_in_HCC10.png")
nrow(out); colSums(out[,2:7])
out <- Common_markers(HCC10, HCC10_markers, Extern_table[Extern_table[,2]=="Mature",1], "Extern_C-Mature_in_HCC10.png")
nrow(out); colSums(out[,2:7])
out <- Common_markers(HCC10, HCC10_markers, Extern_table[Extern_table[,2]=="Chol_HPA",1], "Extern_Chol-HPA_in_HCC10.png")
nrow(out); colSums(out[,2:7])
out <- Common_markers(HCC10, HCC10_markers, Extern_table[Extern_table[,2]=="Hep_HPA",1], "Extern_Hep-HPA_in_HCC10.png")
nrow(out); colSums(out[,2:7])


######## Trio Plot ##########
# Chol marker = ERO1A
# Hep marker 
trio_plot <- function(gene, title=FALSE) {
	vals <- c(exprs(CCA1)[fData(CCA1)$feature_symbol == gene,], exprs(HCC6)[fData(HCC6)$feature_symbol == gene,], exprs(HCC10)[fData(HCC10)$feature_symbol == gene,])
	splits <- quantile(vals, probs=seq(from=0, to=1, length=8))
	splits <- unique(splits);
        #splits <- seq(from=0, to=max(splits), length=5)
	bins1 <- cut(exprs(CCA1)[fData(CCA1)$feature_symbol == gene,], breaks=splits, include.lowest=TRUE);
	bins2 <- cut(exprs(HCC6)[fData(HCC6)$feature_symbol == gene,], breaks=splits, include.lowest=TRUE);
	bins3 <- cut(exprs(HCC10)[fData(HCC10)$feature_symbol == gene,], breaks=splits, include.lowest=TRUE);
	my_cols <- colorRampPalette(brewer.pal(6, "Blues"))(length(splits)-1);
	if (title) {
		title1 <- "CCA1"
		title2 <- "HCC6"
		title3 <- "HCC10"
	} else {
		title1 <- ""
		title2 <- ""
		title3 <- ""
	}

	plot(CCA1$dm1, -CCA1$dm2, col=my_cols[bins1], pch=16, main=title1, xaxt="n", yaxt="n", xlab="", ylab=gene, bty="n")
	plot(-HCC6$dm1, HCC6$dm2, col=my_cols[bins2], main=title2, pch=16, xaxt="n", yaxt="n", xlab="", ylab="", bty="n")
	plot(-HCC10$dm1, HCC10$dm2, col=my_cols[bins3], main=title3, pch=16, xaxt="n", yaxt="n", xlab="", ylab="", bty="n")
}

png("MarkerPlots_forMelbourne.png", width=8, height=6, units="in", res=300)
par(mfrow=c(2,3))
par(mar=c(1,4,3,1))
trio_plot("ERO1A", title=TRUE)
trio_plot("CYP27A1", title=FALSE)
dev.off()




# OLAPS
J <- function(A, B) {
	A_set <- unique(Extern_table[Extern_table[,2] == A,1])
	B_set <- unique(Extern_table[Extern_table[,2] == B,1])
	both = sum(A_set %in% B_set);
	j = both/ (length(A_set)+length(B_set)-both);
	return(both);
}
for(i in unique(Extern_table[,2])){
print(sapply(unique(Extern_table[,2]), J, B=i))
}




Fancy_SC <- HCC6_sc_agree[HCC6_sc_agree %in% CCA1_sc_agree & HCC6_sc_agree %in% HCC10_sc_agree]
# Best common markers = "ENSG00000164104" "ENSG00000173207"
Fancy_Hep <- rownames(out8)[1:2]
#Fancy_Chol <- 

fData(HCC6)$feature_symbol[rownames(HCC6) %in% Fancy_Hep]


intensity_plot(CCA1, "HMGB2")
intensity_plot(HCC6, "HMGB2")
intensity_plot(HCC10, "HMGB2")



##### Marker Breakdown ######
Camp <- readRDS("/lustre/scratch117/cellgen/team218/TA/scRNASeqDatasets/Liver/camp_ht_qced.rds") # QCed healthy human tissue
Camp_markers <- complex_markers(exprs(Camp), factor(Camp$cell_type1)) 
Camp_progenitor_markers <- rownames(Camp_markers)[Camp_markers$q.value < 0.05 & Camp_markers$AUC > 0.7 & Camp_markers[,2] == 0 & Camp_markers[,4] ==1 & rowSums(Camp_markers[,2:8]) == 1]
group_means <- my_row_mean_aggregate(exprs(Camp), factor(Camp$cell_type1))


CCA1_good_markers <- CCA1_markers[CCA1_markers$q.value < 0.05 & CCA1_markers$AUC > 0.7 & rowSums(CCA1_markers[,2:7]) < 3 & CCA1_markers[,2]+CCA1_markers[,5] > 0,]
HCC6_good_markers <- HCC6_markers[HCC6_markers$q.value < 0.05 & HCC6_markers$AUC > 0.7 & rowSums(HCC6_markers[,2:7]) < 3 & HCC6_markers[,5] > 0,]
HCC10_good_markers <- HCC10_markers[HCC10_markers$q.value < 0.05 & HCC10_markers$AUC > 0.7 & rowSums(HCC10_markers[,2:7]) < 3 & HCC10_markers[,2] > 0,]

Consistent_good_sc_markers <- rownames(HCC6_good_markers)[rownames(HCC6_good_markers) %in% rownames(CCA1_good_markers) & rownames(HCC6_good_markers) %in% rownames(HCC10_good_markers)]
Consistent_good_sc_markers_names <- unique(fData(CCA1)[rownames(CCA1) %in% Consistent_good_sc_markers,"feature_symbol"])

cellcycle <- read.table("~/Data/Whitfield_CC.txt")



intensity_plot(CCA1, "CYP2E1")
intensity_plot(HCC6, "CYP2E1")
intensity_plot(HCC10, "CYP2E1")

###### Marker Gene List Analysis #######
# (1) Gene-Gene correlations among list in each donor & across donors
# (2) cluster genes based on correlations -> core genes in the list
# 
#
# Or just glm results? +1,-1,0 foreach coefficient - across types and within types = 18-item vector + 3 item-vector -> recluster genes. 


# Slalom - requires R3.4
library("slalom")
genesets <- GSEABase::getGmt("~/Collaborations/LiverOrganoids/Gene-sets.gmt.csv")

my_convert <- function(SCE) {
	dat <- exprs(SCE)
	keep <- fData(SCE)$feature_symbol != ""
	dat <- dat[keep,]
	rownames(dat) <- fData(SCE)$feature_symbol[keep]
	SingleCellExperiment::SingleCellExperiment(assays = list(logcounts=dat))
}

cca1_sce <- my_convert(CCA1)
hcc6_sce <- my_convert(HCC6)
hcc10_sce <- my_convert(HCC10)
all_sce <- SingleCellExperiment::SingleCellExperiment(assays = list(logcounts=cbind(exprs(CCA1), exprs(HCC6), exprs(HCC10))))

set.seed(162)
cca1_model <- newSlalomModel(cca1_sce, genesets, n_hidden = 5, min_genes = 5)
cca1_model <- initSlalom(cca1_model)
cca1_model <- trainSlalom(cca1_model, nIterations = 5000, shuffle=TRUE, pretrain=TRUE, seed=666)
topTerms(cca1_model)
saveRDS(cca1_model, "Slalom_CCA1.rds")

set.seed(162)
hcc6_model <- newSlalomModel(hcc6_sce, genesets, n_hidden = 5, min_genes = 5)
hcc6_model <- initSlalom(hcc6_model)
hcc6_model <- trainSlalom(hcc6_model, nIterations = 5000, shuffle=TRUE, pretrain=TRUE, seed=666)
topTerms(hcc6_model)
saveRDS(hcc6_model, "Slalom_HCC6.rds")

set.seed(162)
hcc10_model <- newSlalomModel(hcc10_sce, genesets, n_hidden = 5, min_genes = 5)
hcc10_model <- initSlalom(hcc10_model)
hcc10_model <- trainSlalom(hcc10_model, nIterations = 5000, shuffle=TRUE, pretrain=TRUE, seed=666)
topTerms(hcc10_model)
saveRDS(hcc10_model, "Slalom_HCC10.rds")




