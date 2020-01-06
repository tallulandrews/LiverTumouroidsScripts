# GO enrichments of markers
# External markers
# diff CC phase of same functional cell-type?
source("/nfs/users/nfs_t/ta6/Collaborations/LiverOrganoids/Laura_Pipeline/0_ColourScheme.R")
source("/nfs/users/nfs_t/ta6/R-Scripts/Blank_plot.R")

args <- commandArgs(trailingOnly=TRUE) # SCE RDSs
clusters_name <- "clusters_clean"
if (!grepl("rds$", args[1])) {
	clusters_name <- args[1]
	args <- args[-1]
}
nSCEs <- length(args)
expr_type <- "lognorm"

type_col <- type_col[1:nSCEs]

SCE_list <- list();
max_ngenes <- 0;
for (f in args) {
	require("scater")
	obj <- readRDS(f);
	if (class(obj)[1] == "SCESet") {
		obj <- toSingleCellExperiment(obj)
	}
	if (nrow(obj) > max_ngenes) {max_ngenes <- nrow(obj)}
	tmp <- unlist(strsplit(f, "\\."))
	SCE_list[[tmp[1]]] <- obj;
}

# Laura_LitMarkers
Chol_lineage <- read.table("/nfs/users/nfs_t/ta6/Collaborations/LiverOrganoids/Markers_130418_Chol.txt", header=TRUE)
Hep_lineage <- read.table("/nfs/users/nfs_t/ta6/Collaborations/LiverOrganoids/Markers_130418_Hep.txt", header=TRUE)
Chol_lineage[,2] <- as.character(Chol_lineage[,2])
Chol_lineage[Chol_lineage[,2] == "Prog",2] <- "Chol-Prog"
Chol_lineage[Chol_lineage[,2] == "Chol",2] <- "Chol-Mature"
Hep_lineage[,2] <- as.character(Hep_lineage[,2])
Hep_lineage[Hep_lineage[,2] == "Prog",2] <- "Hep-Prog"
Hep_lineage[Hep_lineage[,2] == "Hep",2] <- "Hep-Mature"
Lineage <- rbind(Chol_lineage, Hep_lineage)
Liver <- intersect( Lineage[ Lineage[,2] == "Chol-Mature", 1], Lineage[ Lineage[,2] == "Hep-Mature", 1] )
Prog <- intersect( Lineage[ Lineage[,2] == "Chol-Prog", 1], Lineage[ Lineage[,2] == "Hep-Prog", 1] )
Chol <- intersect( Lineage[ Lineage[,2] == "Chol-Prog", 1], Lineage[ Lineage[,2] == "Chol-Mature", 1] )
Hep <- intersect( Lineage[ Lineage[,2] == "Hep-Prog", 1], Lineage[ Lineage[,2] == "Hep-Mature", 1] )
exclude <- c(intersect( Lineage[ Lineage[,2] == "Hep-Prog", 1], Lineage[ Lineage[,2] == "Chol-Mature", 1] ),
	     intersect( Lineage[ Lineage[,2] == "Chol-Prog", 1], Lineage[ Lineage[,2] == "Hep-Mature", 1] ))
Lineage[ Lineage[,1] %in% Liver, 2] <- "Mature-Liver"
Lineage[ Lineage[,1] %in% Prog, 2] <- "Common-Prog"
Lineage[ Lineage[,1] %in% Chol, 2] <- "Common-Chol"
Lineage[ Lineage[,1] %in% Hep, 2] <- "Common-Hep"
Lineage <- Lineage[ !Lineage[,1] %in% exclude,]
Lineage[,2] <- factor(Lineage[,2])
Lineage <- Lineage[!duplicated(Lineage[,1]),]


# Cell-Type Markers
require("CellTypeProfiles")
source("/nfs/users/nfs_t/ta6/NetworkInferencePipeline/Dropouts/My_R_packages/CellTypeProfiles/R/Markers.R")
require("gProfileR")
require("proxy")
n_top_rich = 30;
for (i in 1:nSCEs) {
	obj <- SCE_list[[i]]
	obj_orig <- obj;

	background_genes_spec <- rownames(obj)[ rowData(obj)$pct_dropout <= 99];

	obj <- obj[rownames(obj) %in% background_genes_spec,]
	obj <- obj[, colData(obj)[,clusters_name] != "Outliers"]
	markers1 <- complex_markers(assays(obj)[[ expr_type ]], factor(colData(obj)[,clusters_name]), n_max=1)
	m <- markers1[match(rownames(obj_orig), rownames(markers1)),]
	colnames(m) <- paste("markers", clusters_name, colnames(m), sep="_")
	rowData(obj_orig) <- cbind(rowData(obj_orig), m)

	rowData(obj_orig)$Lineage <- Lineage[ match(rowData(obj_orig)$Symbol, Lineage[,1]) , 2]

	SCE_list[[i]] <- obj_orig;

	markers1 <- markers1[order(markers1$AUC, decreasing=TRUE),]
	markers1 <- markers1[ markers1$q.value < 0.05 & markers1$q.value >= 0 , ]
	
	# get gene list

	get_top_markers <- function(group) {
		m <- markers1[markers1$Group == group,]
		m <- m[1:10,]
		return(m)
	}
	get_richments <- function(group) {
		gene_list <- rownames(markers1)[markers1$Group == group]
		   enrichments <- gprofiler(gene_list, organism="hsapiens", 
		   ordered_query=T, significant=T, custom_bg=background_genes_spec, 
		   hier_filtering="moderate", max_set_size=10000,
		   src_filter=c("GO:BP", "KEGG", "REAC", "HPA"), 
		   correction_method="fdr", min_isect_size=3, min_set_size=10)
		enrichments<-enrichments[order(enrichments$p.value),]
		# enrichment filtering
		# remove hpa low
		exclude <- enrichments$domain == "hpa" & grepl("Low", enrichments$term.name)
		exclude <- exclude | enrichments$domain == "hpa" & grepl("Not detected", enrichments$term.name)
		exclude <- exclude | enrichments$domain == "hpa" & grepl("Uncertain", enrichments$term.name)
		enrichments <- enrichments[!exclude,]
		enrichments <- enrichments[1:min(nrow(enrichments), n_top_rich),]
		enrichments$GroupID = rep(group, times=nrow(enrichments))
		return(enrichments)
	}
	top_marker_table <- factor()
	Rich_table <- factor()
	for(g in levels(factor(markers1$Group))) {
		top_marker_table <- rbind(top_marker_table, get_top_markers(g))
		Rich_table <- rbind(Rich_table, get_richments(g))
	}	
	write.table(Rich_table, file=paste(names(SCE_list)[i], clusters_name, "marker_GOrich.txt", sep="_"))
	write.table(top_marker_table, file=paste(names(SCE_list)[i], clusters_name, "TopMarkerTable.txt", sep="_"))
}



for (i in 1:length(SCE_list)) {
	saveRDS(SCE_list[[i]], file=paste(names(SCE_list)[i], "Function.rds", sep="_"))
}
