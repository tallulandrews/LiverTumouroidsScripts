source("~/Collaborations/LiverOrganoids/Laura_Pipeline/0_ColourScheme.R")

intensity_plot <- function(SCE1, SCE2, genes, name1="", name2="", ylab1="", ylab2="") {
        vals1 <- exprs(SCE1)[fData(SCE1)$Symbol %in% genes,];
	if (!is.null(dim(vals1))) {
		vals1 <- colMeans(vals1)
	}

        vals2 <- exprs(SCE2)[fData(SCE2)$Symbol %in% genes,];
	if (!is.null(dim(vals2))) {
		vals2 <- colMeans(vals2)
	}
        splits <- quantile(c(vals1, vals2), probs=seq(from=0, to=1, length=8))
        splits <- unique(splits);
#        splits <- seq(from=0, to=max(splits), length=7)

        bins <- cut(vals1, breaks=splits, include.lowest=TRUE);
        my_cols <- colorRampPalette(brewer.pal(8, "Blues"))(length(levels(bins)));
        plot(SCE1$coords1, SCE1$coords2, col=my_cols[bins], pch=16, main=name1, xaxt="n", yaxt="n", ylab=ylab1, xlab="", bty="n")

        bins <- cut(vals2, breaks=splits, include.lowest=TRUE);
        my_cols <- colorRampPalette(brewer.pal(8, "Blues"))(length(levels(bins)));
        plot(SCE2$coords1, SCE2$coords2, col=my_cols[bins], pch=16, main=name2, xaxt="n", yaxt="n", ylab=ylab2, xlab="", bty="n")
	return(splits)
}


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
		  HCC10="HCC10_manual_SC32.rds",
		  HCC23="HCC23_manual_SC3.rds",
		  HCC24="HCC24_manual_SC3.rds",
		  D3DM="D3DM_manual_SC3.rds",
		  D3EM="D3EM_manual_SC3.rds")
#		  D9DM="D9DM_manual_SC3.rds",
#		  D9EM="D9EM_manual_SC3.rds")

dim_reduction <- list(CCA1="CCA1_1000_Visualizations_dims.rds",
		      CCA5="CCA5_1000_Visualizations_dims.rds",
		      HCC6="HCC6_1000_Visualizations_dims.rds",
		      HCC10="HCC10_1500_Visualizations_dims.rds",
		      HCC23="HCC23_1000_Visualizations_dims.rds",
		      HCC24="HCC24_500_Visualizations_dims.rds",
		      D3DM="D3DM_1000_Visualizations_dims.rds",
		      D3EM="D3EM_1000_Visualizations_dims.rds")
#		      D9DM="D9DM_1000_Visualizations_dims.rds",
#		      D9EM="D9EM_1000_Visualizations_dims.rds")

line_specific_genes <- list(CCA1=c("ATP1B3", "DPAGT1", "CLDN2", "AQP5", "EZH2", "RECQL4", "TRAIP", "LMNB1", "CDCA7L", "IQGAP3", "DDIAS", "CA9", "NDRG1", "ATP2B4", "DAPK1", "HIST1H2AC"),
			    CCA5=c("HMGB2", "LDHA", "CD81", "ANXA4"),
			    HCC6=c("MCM7", "GAL", "CST3", "IGFBP3"),
			    HCC23=c("HMGB2", "F10", "SERPINA1", "CBS"),
			    HCC10=c("PLK4", "AURKB", "MKI67", "EZH2", "LMNB1", "HMGB2", "KPNA2", "LRR1", "TTK", "SLC1A5", "SLC19A1", "SLC1A4", "SLCO4A1", "GPM6A", "MXD4", "AFM", "ATP5L2", "GLI4", "SLC29A4", "HEXIM1", "PRRG3", "PBXIP1", "CAPNS1", "CLDN11", "HMGCS2", "MLXIPL", "ABAT", "ADH4"),
			    HCC24=c("NDC80", "MTFR2", "MAGEB1", "DIAPH3", "IQGAP3", "RAD51AP1", "TRAIP", "HMGB2", "DEPDC1B", "LRR1", "SLC19A3", "FRAS1", "MKI67", "DDIAS", "ITIH4", "IGFBP1", "ADH4", "GC", "AZGP1", "CD36", "CD44", "UPB1", "NR1H4", "HSD17B2", "TNFSF15", "TAT", "CYP2C19"),
			    D3DM=c("CDK2AP1","CD24","TESC","P4HA1"),
			    D3EM=c("YWHAB","ENO2","BTG1","MIF"),
			    D9DM=c("CEACAM6","SERPINA1","NDRG1","KCNT1"),
			    D9EM=c("EIF5A","HERPUD1","CA9","EFNA1")
			) # genes for dotplots


line_specific_groups <- list(CCA1=c("Progenitor", "Differentiated1", "TICs", "Differentiated2"),
			    CCA5=c("Chol", "Stress", "CSC", "Unk", "Hep"),
			    HCC6=c("Prog1", "Stress", "Prog2", "CSC"),
			    HCC23=c("CSC", "Clot-Hep", "Prog"),
			    #HCC10=c("Prog1", "Hep", "CSC", "iHep", "Stress", "Prog2"),
			    HCC10=c("TICs", "Progenitor", "Differentiated"),
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
	this_name <- paste(i, dimred_name, "Bespoke", sep="_")
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
	pdf(paste(this_name, "DimRedScatter.pdf", sep="_"), width=6, height=6)
	plot(coords$x[cell_keep], coords$y[cell_keep], col=cell_colours, pch=16, xlab="Dim 1", ylab="Dim 2")
	dev.off()

	# DotPlot
	require("ggplot2")
	pdf(paste(this_name, "MarkerDot.pdf", sep="_"))
	seurat <- as.Seurat(SCE, data=expr_type)
	a<-DotPlot(seurat, features=line_specific_genes[[i]], group.by="named_clusters")
	plot(a)
	dev.off()

	# ScatterPlot + scmap results
	pdf(paste(this_name, "ScmapScatter.pdf", sep="_"), width=6, height=6)
	scmap_out <- readRDS(scmap_results[[i]])
	plot(coords$x, coords$y, col="black", bg=scmap_out$scmap_cell_Cols, pch=21);
	lgend <- unique(scmap_out$scmap_cluster_labs);
	lgend_col <- unique(scmap_out$scmap_cell_Cols);
	dev.off()
}
