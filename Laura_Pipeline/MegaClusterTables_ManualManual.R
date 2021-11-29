
####### Mega Marker Tables #######
source("~/Collaborations/LiverOrganoids/Laura_Pipeline/99_Functions.R")
require("CellTypeProfiles")
require("SingleCellExperiment")
require("scater")

CC_genes <- load_CC(set="cycling")
CC_genes <- c(as.character(CC_genes$Whitfield[,2]), as.character(CC_genes$Tirosh[,1]))

GO <- read.delim("/lustre/scratch117/cellgen/team218/TA/Mirrored_Annotations/GO/hsapiens_80_GO_Annoations_Emsembl.out", sep="\t", header=F)
GO_cc <- GO[GO[,3] == "cell cycle",1]

schwalie_StemCell <- c("MFAP2", "HSPB6", "GBP3", "MYC", "CCND2", "RAC2", "ZFP36L2", "VGLL4", "ADGRG1", "IDH2", "CDK6", "SH3BGRL", "IFITM3", "LCP1", "ETV6", "RPGRIP1", "ERI3", "SDCBP", "DAP", "JTB", "GSTO1", "DAPP1", "CTSZ")
GO_telo <- unique(GO[GO[,3] == "telomere maintenance",1])
LigRec <- read.delim("~/Data/LigandReceptorPairs.csv", sep=",", header=T)
GO_tf <- unique(GO[grep("transcription factor activity", GO[,3]),1])

#Lineage <- read.table("/lustre/scratch117/cellgen/team218/TA/LiverOrganoids/Analysis_Pipeline_Output/Cleaned_Lineage.txt")
Lineage <- read.table("/nfs/users/nfs_t/ta6/Collaborations/LiverOrganoids/Cleaned_PatLineage.txt")

Surface <- read.delim("~/Data/Combined_Surfaceome.csv", sep=",")

add_anno <- function(marker_tab) {
        marker_tab$is.cycle <- rownames(marker_tab) %in% GO_cc
        marker_tab$is.telo <- rownames(marker_tab) %in% GO_telo
        marker_tab$is.tf <- rownames(marker_tab) %in% GO_tf
        marker_tab$ReceptorOf <- LigRec[match(marker_tab$Symbol,LigRec$Receptor.ApprovedSymbol),"Ligand.ApprovedSymbol"]
        marker_tab$LigandTo <- LigRec[match(marker_tab$Symbol,LigRec$Ligand.ApprovedSymbol),"Receptor.ApprovedSymbol"]
        marker_tab$Lineage <- Lineage[match(marker_tab$Symbol,Lineage[,1]),2]
	tmp <- Surface[match(marker_tab$Symbol,rownames(Surface)),]
	tmp[is.na(tmp)] <- 0
        marker_tab <- cbind(marker_tab, tmp)
        return(marker_tab)
}

# Each line:
#	GeneID	GeneSymbol means, cell-type marker-ness, annotations, biotype
# Cluster table:
#	ClusterID ClusterName nCells, CCStage, nSigMarkers, most similar clusters, 
clustered_rds <- c("CCA1_manual_SC3.rds", 
		"CCA5_manual_SC3.rds", 
		"HCC6_manual_SC3.rds", 
		"HCC23_manual_SC3.rds", 
		"HCC10_manual_SC32.rds", 
		"HCC24_manual_SC3.rds", 
		"D3DM_manual_SC3.rds", 
		"D3EM_manual_SC3.rds", 
		"D9DM_manual_SC3.rds", 
		"D9EM_manual_SC3.rds");
rds_names <- c("CCA1", "CCA5", "HCC6", "HCC23", "HCC10", "HCC24", "D3DM", "D3EM", "D9DM", "D9EM")
expr_type <- "norm_exprs";

line_specific_groups <- list(CCA1=c("Progenitor", "Differentiated1",
                        "TICs", "Differentiated2"),
                        CCA5=c("Differentiated1", "Differentiated1", "TICs", "Quiescent", "Differentiated2"),
                        HCC10=c("TICs", "Progenitor", "Differentiated"),
                        HCC23=c("TICs", "Differentiated", "Progenitor"),
                        HCC6=c("Progenitor1", "Differentiated", "Progenitor2", "TICs"),
                        HCC24=c("Differentiated", "Progenitor1", "Progenitor2", "TICs")) # cluster names


Global <- readRDS("../Global_SCE.rds")

make_cluster_marker_tables <- function(line_name) {
        # Set up
	i <- which(rds_names == line_name);
	file = clustered_rds[i]

        SCE <- readRDS(paste("../", file, sep=""))
	# Cluster stats
	if (!identical(colnames(Global)[Global$Donor == line_name], colnames(SCE))) {stop("Global-local mismatch")}  # checking
	Cs_cc <- table(Global$CC_state[Global$Donor == line_name], SCE$Manual_Clusters)

        nCs <- factor_counts(SCE$Manual_Clusters)
        mark_cells <- SCE$Manual_Clusters %in% names(nCs)[nCs > 10]
	SCE <- SCE[,mark_cells]
	Cs_cc <- Cs_cc[,nCs > 10]

	Cs_names <- line_specific_groups[[line_name]]
	names(Cs_names) <- names(nCs)[nCs > 10]
	SCE$Manual_Manual_Clusters <- Cs_names[match(SCE$Manual_Clusters, names(Cs_names))]

	nCs <- factor_counts(SCE$Manual_Clusters)

	Cs_ids <- names(nCs);

        # Calculate Markers
	# change clusters
        mark_mat <- assays(SCE)[[expr_type]]
        marks <- complex_markers(mark_mat, factor(SCE$Manual_Manual_Clusters))
	marks$is.GoodMarker <- marks$AUC > 0.75 & marks$q.value < 0.05

	# Markers V2 (scfind)
	#index <- buildCellTypeIndex(SCE, dataset.name=line_name, cell.type.label="Manual_Clusters")
	#find_marks <- evaluateMarkers(index, rownames(SCE), unique(SCE$Manual_Clusters))

	# Mean Expression
	expr_table <- my_row_mean_aggregate(assays(SCE)[["logcounts"]], SCE$Manual_Manual_Clusters)

	# Assemble Marker Table
	MARKERS <- cbind(marks, expr_table, rowData(SCE)$feature_symbol, rowData(SCE)$biotype)
	colnames(MARKERS)[grep("symbol", colnames(MARKERS))] <- "Symbol"
	colnames(MARKERS)[grep("biotype", colnames(MARKERS))] <- "GeneType"
	MARKERS <- add_anno(MARKERS); 
	MARKERS <- MARKERS[MARKERS$AUC != -1 ,]

	write.table(MARKERS, file=paste(line_name, "ManualManualClustering_MarkerTable_Nov2018.csv", sep="_"), sep=",")

	Clust_marks <- MARKERS[MARKERS$is.GoodMarker, ]
	unique_marks <- Clust_marks[,2:( ncol(marks) -3)]
	unique_marks <- unique_marks[rowSums(unique_marks) == 1,]
	Cs_marks <- colSums(unique_marks)

	# Lineage markers by cluster
	C_assigned <- apply(unique_marks, 1, function(x){which(x==1)})
	#C_assigned <- factor(as.character(C_assigned), levels=colnames(Cs_cc))
	C_assigned <- colnames(unique_marks)[C_assigned]
	Cs_marks_lin <- table(MARKERS$Lineage[rowSums(MARKERS[,2:( ncol(marks) -3)]) == 1 & MARKERS$is.GoodMarker], C_assigned)

	# Cluster-cluster similarity
	c_expr_table <- expr_table[marks$is.GoodMarker,]
	Cs_cSim <- cor(c_expr_table)
	rownames(Cs_cSim) <- paste("Cor_with", rownames(Cs_cSim), sep="_")
		
	# Assemble Cluster Table
	#Cs_names <- Cs_names[match(names(nCs), names(Cs_marks))]
	nCs <- table(SCE$Manual_Manual_Clusters)
	nCs <- nCs[match(Cs_names, names(nCs))]
	Cs_marks_lin <- Cs_marks_lin[,match(Cs_names, colnames(Cs_marks_lin))]
	Cs_cSim <- Cs_cSim[,match(Cs_names, names(Cs_marks))]
	Cs_marks <- Cs_marks[match(Cs_names, names(Cs_marks))]

	CLUSTERS <- rbind(Cs_ids, Cs_names, nCs, Cs_cc, Cs_marks, Cs_marks_lin, round(Cs_cSim, digits=2))
	rownames(CLUSTERS)[rownames(CLUSTERS) == "nCs"] <- "n_cells"
	rownames(CLUSTERS)[rownames(CLUSTERS) == "Cs_ids"] <- "ID_Num"
	rownames(CLUSTERS)[rownames(CLUSTERS) == "Cs_names"] <- "Suggest_Name"
	rownames(CLUSTERS)[rownames(CLUSTERS) == "Cs_marks"] <- "Unique_Markers"
	CLUSTERS <- CLUSTERS[rownames(CLUSTERS) != "G0",]
	
	write.table(CLUSTERS, file=paste(line_name,"ManualManualClustering_ClusterTable.csv", sep="_"), sep=",")
}

tmp_names <- rds_names[2:10]
for (line in tmp_names) {
	print(line)
	make_cluster_marker_tables(line)
}

