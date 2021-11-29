Hepa_Dimension <- readRDS("Camp_Data/Hepatocyte_Profile.rds")
Chol_Dimension <- readRDS("ExternalData/Sampaziotis_CholDimension.rds")

require("M3Drop")
require("scater")
require("matrixStats")
require("RColorBrewer")
require("CellTypeProfiles")
source("complex_markers.R")


set.seed(142)

CCA1 <- readRDS("CCA1_SC3_Prolif.rds")
HCC6 <- readRDS("HCC6_SC3_Prolif.rds")
HCC10 <- readRDS("HCC10_SC3_Prolif.rds")

source("/nfs/users/nfs_t/ta6/R-Scripts/Ensembl_Stuff.R")
fData(CCA1)$feature_symbol <- map_symbol_ensg(rownames(counts(CCA1)), is.org="Hsap", is.name="ensg")
fData(HCC6)$feature_symbol <- map_symbol_ensg(rownames(counts(HCC6)), is.org="Hsap", is.name="ensg")
fData(HCC10)$feature_symbol <- map_symbol_ensg(rownames(counts(HCC10)), is.org="Hsap", is.name="ensg")


my_custom_projection <- function(SCE) {
	require("RColorBrewer")

	dat <- exprs(SCE)
	scaled <- t(apply(dat, 1, function(x){(x-mean(x))/sd(x)}))
	keep = fData(SCE)$feature_symbol != "" & !is.na(rowSums(scaled))
	scaled <- scaled[keep,]
	g_sym <- fData(SCE)$feature_symbol[keep]	
	markers <- complex_markers(exprs(SCE), SCE$sc3_6_clusters)
	
	keep <- rownames(markers)[markers$q.value < 0.05]
	keep <- rownames(scaled) %in% keep
	scaled <- scaled[keep,]
	g_sym <- g_sym[keep]
	
	hep <- t(scaled[match(names(Hepa_Dimension$HepUp), g_sym),]) * rank(Hepa_Dimension$HepUp)
	hsc <- t(scaled[match(names(Hepa_Dimension$SCUp), g_sym),]) * rank(Hepa_Dimension$SCUp)
	chol <- t(scaled[match(as.character(Chol_Dimension$CholUp[,1]), g_sym),]) * rank(Chol_Dimension$CholUp[,2])
	csc <- t(scaled[match(as.character(Chol_Dimension$hPSCsUp[,1]), g_sym),]) * rank(Chol_Dimension$hPSCsUp[,2])

	labs <- factor(SCE$sc3_6_clusters)

	plot(rowMeans(hep, na.rm=TRUE)-rowMeans(chol, na.rm=TRUE), rowMeans(hsc, na.rm=TRUE)+rowMeans(csc, na.rm=TRUE),
		xlab="Cholangiocyte vs Hepatocyte", ylab="Stem-ness", pch=16, col=brewer.pal(6, "Set3")[labs])
	legend("topright", fill=brewer.pal(6, "Set3")[1:length(levels(labs))], levels(labs), bty="n")

}

png("Hep_Chol_Score_CCA1_Plot.png", width=6, height=6, units="in", res=300)
my_custom_projection(CCA1)
dev.off()
png("Hep_Chol_Score_HCC6_Plot.png", width=6, height=6, units="in", res=300)
my_custom_projection(HCC6)
dev.off()
png("Hep_Chol_Score_HCC10_Plot.png", width=6, height=6, units="in", res=300)
my_custom_projection(HCC10)
dev.off()
	#Hscore <- colMeans( scaled[g_sym %in% names(Hepa_Dimension$HepUp),])
	#HSscore <- colMeans(scaled[g_sym %in% names(Hepa_Dimension$SCUp),])
	#Cscore <- colMeans(scaled[g_sym %in% as.character(Chol_Dimension$CholUp[,1]),])
	#CSscore <- colMeans(scaled[g_sym %in% as.character(Chol_Dimension$hPSCsUp[,1]),])
pd <- new("AnnotatedDataFrame", data.frame(sc3_6_clusters=rep(c("CCA1","HCC6","HCC10"), times=c(ncol(CCA1), ncol(HCC6), ncol(HCC10)))))
allcounts <- cbind(counts(CCA1), counts(HCC6), counts(HCC10))
rownames(pd) <- colnames(allcounts)
Whole_data <- newSCESet(countData=allcounts, phenoData=pd)
fData(Whole_data)$feature_symbol <- fData(CCA1)$feature_symbol

png("Hep_Chol_Score_All_Plot.png", width=6, height=6, units="in", res=300)
my_custom_projection(Whole_data)
dev.off()





