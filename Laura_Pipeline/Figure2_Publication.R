source("~/Collaborations/LiverOrganoids/Laura_Pipeline/0_ColourScheme.R")



# Read in each manual clustering
expr_type <- "norm_exprs";
expr_mats <- list()
cell_colors <- list();
nn_graph <- list();
cell_coords <- list();

dimred_name <- "dm"

sce_objs <- list(CCA1="CCA1_manual_SC3.rds",
		  CCA5="CCA5_manual_SC3.rds",
		  HCC6="HCC6_manual_SC3.rds",
		  HCC10="HCC10_manual_SC3.rds",
		  HCC23="HCC23_manual_SC3.rds",
		  HCC24="HCC24_manual_SC3.rds",
		  D3DM="D3DM_manual_SC3.rds",
		  D3EM="D3EM_manual_SC3.rds",
		  D9DM="D9DM_manual_SC3.rds",
		  D9EM="D9EM_manual_SC3.rds")

dim_reduction <- list(CCA1="CCA1_1000_Visualizations_dims.rds",
		      CCA5="CCA5_1000_Visualizations_dims.rds",
		      HCC6="HCC6_1000_Visualizations_dims.rds",
		      HCC10="HCC10_1000_Visualizations_dims.rds",
		      HCC23="HCC23_1000_Visualizations_dims.rds",
		      HCC24="HCC24_1000_Visualizations_dims.rds",
		      D3DM="D3DM_1000_Visualizations_dims.rds",
		      D3EM="D3EM_1000_Visualizations_dims.rds",
		      D9DM="D9DM_1000_Visualizations_dims.rds",
		      D9EM="D9EM_1000_Visualizations_dims.rds")

line_specific_genes <- list(CCA1=c("ZWINT", "NDRG1", "EIF5A", "HERPUD1"),
			    CCA5=c("HMGB2", "LDHA", "CD81", "ANXA4"),
			    HCC6=c("MCM7", "GAL", "CST3", "IGFBP3"),
			    HCC23=c("HMGB2", "F10", "SERPINA1", "CBS"),
			    HCC10=c("HMGB2","CEBPA","ABCA1","AFP"),
			    HCC24=c("ZWINT","ALB","KLK1","TUSC3"),
			    D3DM=c("CDK2AP1","CD24","TESC","P4HA1"),
			    D3EM=c("YWHAB","ENO2","BTG1","MIF"),
			    D9DM=c("CEACAM6","SERPINA1","NDRG1","KCNT1"),
			    D9EM=c("EIF5A","HERPUD1","CA9","EFNA1")
			) # genes for dotplots


line_specific_groups <- list(CCA1=c("Prog", "Chol", "CSC", "Stress"),
			    CCA5=c("Chol", "Stress", "CSC", "Unk", "Hep"),
			    HCC6=c("Prog1", "Stress", "Prog2", "CSC"),
			    HCC23=c("CSC", "Clot-Hep", "Prog"),
			    HCC10=c("Prog1", "Hep", "CSC", "iHep", "Stress", "Prog2"),
			    HCC24=c("Clot-Hep", "Prog1", "Prog2", "CSC", "Hep1", "Hep2"),
			    D3DM=c("Chol", "Prog", "Stress", "Cycling"),
			    D3EM=c("Stress", "Chol1", "Chol2", "Prog"),
			    D9DM=c("Prog1", "Prog2", "Prog3", "Clot-Hep", "Chol1", "Chol2"),
			    D9EM=c("Chol1", "Cycling", "Chol2")
			) # cluster names
			

scmap_results <- list(CCA1="CCA1_scmap_output.rds",
		      CCA5="CCA5_scmap_output.rds",
		      HCC10="HCC10_scmap_output.rds",
		      HCC6="HCC6_scmap_output.rds",
		      HCC23="HCC23_scmap_output.rds",
		      HCC24="HCC24_scmap_output.rds",
		      D3DM="D3DM_scmap_output.rds",
		      D3EM="D3EM_scmap_output.rds",
		      D9DM="D9DM_scmap_output.rds",
		      D9EM="D9EM_scmap_output.rds"
		)


for (i in names(sce_objs)) {
	# save cluster IDs & markers
	# add lineage annotations to markers
	# do tSNE & sil_nn graph
	# save coords
	require("SingleCellExperiment")
	require("scater")
	require("Rtsne")
	require("M3Drop")
	require("CellTypeProfiles")
	require("Seurat")
	perp <- 30
	# Set up
	SCE <- readRDS(sce_objs[[i]])
	SCE <- SCE[!is.na(rowData(SCE)$biotype),]
	SCE <- SCE[rowData(SCE)$feature_symbol != "",]
	SCE <- SCE[!duplicated(rowData(SCE)$feature_symbol),]
	rownames(SCE) <- rowData(SCE)$feature_symbol 
	SCE <- SCE[rowData(SCE)$biotype == "protein_coding",]
	palette <- cluster_col(max(SCE$Manual_Clusters))
	cell_colours <- palette[SCE$Manual_Clusters]
	nCs <- factor_counts(SCE$Manual_Clusters)

	# save colours
	names(palette) <- 1:length(palette)
	SCE@metadata$palette <- palette
	SCE@metadata$C_keep <- nCs > 10
	SCE@metadata$C_names <- line_specific_groups[[i]]
	name_map <- SCE@metadata$C_names
	names(name_map) <- names(SCE@metadata$C_keep)[SCE@metadata$C_keep]


	coords <- readRDS(dim_reduction[[i]])
	
	coords <- coords[[dimred_name]]

	cell_keep <- SCE$Manual_Clusters %in% names(SCE@metadata$C_keep)[SCE@metadata$C_keep]
	SCE <- SCE[,cell_keep]
	SCE$named_clusters <- name_map[match(SCE$Manual_Clusters, names(name_map))]
	
	# ScatterPlot
	plot(coords$x[cell_keep], coords$y[cell_keep], col=cell_colours, cex=16, xlab="Dim 1", ylab="Dim 2")

	# DotPlot
	seurat <- as.Seurat(SCE, data=expr_type)
	DotPlot(seurat, features=line_specific_genes[[i]], group.by="named_clusters")

	# ScatterPlot + scmap results
	scmap_out <- readRDS(scmap_results[[i]])
	plot(coords$x, coords$y, col=scmap_out$scmap_cell_Cols);
	lgend <- unique(scmap_out$scmap_cluster_labs);
	lgend_col <- unique(scmap_out$scmap_cell_Cols);

	
}
