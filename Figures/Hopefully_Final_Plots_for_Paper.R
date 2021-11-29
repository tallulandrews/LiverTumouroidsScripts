source("/nfs/users/nfs_t/ta6/Collaborations/LiverOrganoids/Figures/0_ColourScheme.R")

colour_palette <- list(Chol1="#238b45", Chol2="#66c2a4", 
		Hep1="#cb181d", Hep2="#fb6a4a", 
		Prog1="#fe9929", Prog2="#cc4c02", 
		TICs="#7a0177")


# Read in each manual clustering
expr_type <- "norm_exprs";
expr_mats <- list()
cell_colors <- list();
nn_graph <- list();
cell_coords <- list();

dimred_name <- "dm"

sce_objs <- list(CCA1="CCA1_manual_SC3.rds", CCA5="CCA5_add_ProjObj_Metadata_manual_SC3.rds", 
	HCC10="HCC10_manual_SC32.rds", HCC24="HCC24_manual_SC32.rds", 
	HCC6="HCC6_manual_SC3.rds", HCC23="HCC23_manual_SC3.rds")

dim_reduction <- list(CCA1="CCA1_1000_Visualizations_dims.rds",
			HCC10="HCC10_1500_Visualizations_dims.rds",
			HCC24="HCC24_500_Visualizations_dims.rds",
			CCA5="CCA5_1000_Visualizations_dims.rds",
			HCC23="HCC23_1000_Visualizations_dims.rds",
			HCC6="HCC6_1000_Visualizations_dims.rds")

CC_genes <- c("DDIAS", "TTK", "AURKB", "RECQL4", "TRAIP", "LMNB1", "CDCA7L", 
		"ADCY3", "IQGAP3", "UBE2T", "EXOSC9", "DNAJC9", "CLDN2", "EZH2",
		"AQP5", "ATP1B3", "DPAGT1", "HMGB2", "DAPK1", "ATP2B4", "CA9", "NDRG1",
		"HIST1H2AC")

line_specific_genes <- list(CCA1=c("MKI67", "DDIAS", "TTK", "AURKB", "RECQL4", "TRAIP", "LMNB1", "ADCY3", 
	"IQGAP3", "UBE2T", "EXOSC9", "DNAJC9", "CLDN2", "EZH2", "AQP5", "ATP1B3", "HMGB2", "PLK4",
	"RAD51AP1", "UBE2C", "DAPK1", "ATP2B4", "CA9", "NDRG1", "HIST1H2AC", "ENO2", "VEGFA", 
	"FAM177B", "NDUFA4L2", "LGALS1", "EGLN3", "ARFGEF3"),
			CCA5=c("MKI67", "ASF1B", "PRC1", "TK1", "MCM6", "RNASEH2A", "SLC25A19", "EZH2",
	"LRR1", "HMGB2", "ADCY3", "NCAPH", "UBE2T", "EXOSC9", "DNAJC9", "COG3", "ARNT", "IRGQ",
	"ZBTB11", "MTRR", "SPTLC3", "ANKRD46", "NDRG1", "EGLN3", "ARFGEF3"),
			HCC10=c("MKI67", "PLK4", "AURKB", "NDC80", "TTK", "DIAPH3", "RAD51AP1", "DDIAS", 
	"NCAPH", "UBE2C", "IQGAP3", "SLC1A5", "SLC19A1", "LMNB1", "HMMR", "LRR1", 
	"EZH2","UBE2T", "HMGB2", "KPNA2", "SLC1A4", "SLCO4A1", "GPM6A", "AFM", "ATP5L2", 
	"GLI4", "SLC29A4", "HEXIM1", "PRRG3", "PBXIP1", "CAPNS1", "CLDN11", "HMGCS2", 
	"MLXIPL", "ABAT", "ADH4"),
			HCC24=c("MKI67", "NDC80", "MTFR2", "MAGEB1", "DIAPH3", "DDIAS", "NCAPH", 
	"IQGAP3", "RAD51AP1", "ASF1B", "TTK", "TRAIP", "AURKB", "UBE2C", "RECQL4", "DEPDC1B",
	"HMGB2", "PLK4", "LRR1", "SLC19A3", "FRAS1", "ITIH4", "IGFBP1", "NR1H4","ADH4", "GC", 
	"AZGP1", "HSD17B2", "CD36", "CD44", "TNFSF15", "TAT", "UPB1", "CYP2C19", "VEGFA", 
	 "EGLN3"),
			HCC6=c("MKI67", "NCAPH", "UBE2C", "TTK", "HMMR", "ASF1B", "RECQL4", "PRC1", 
	"IQGAP3", "UBE2T", "DNAJC9", "LRR1", "HMGB2", "EXOSC9", "NDRG1", "HIST1H2AC", "ENO2",
	"NDUFA4L2", "FAM177B", "LGALS1", "IGFBP1", "EGLN3", "ARFGEF3", "AZGP1"),
			HCC23=c("MKI67", "GPR4", "DDIAS", "TTK", "AURKB", "RECQL4", "TRAIP", "LMNB1", 
	"ADCY3", "IQGAP3", "UBE2T", "EXOSC9", "DNAJC9", "EZH2", "ATP1B3", "ASF1B", "PRC1",
	"TK1", "PLK4", "RAD51AP1", "NDC80", "UBE2C", "DAPK1", "ATP2B4", "CA9",
	"NDRG1", "HIST1H2AC", "VEGFA", "IGFBP1", "SLC29A4", "HMGCS2", "MLXIPL", "ABAT",
	"ITIH4", "HSD17B2", "CYP2C19"))

heatmap_genes <- list(CCA1=c("CALM1", "DEGS2", "FASN", "FUT2", "MAP1LC3B", "ROIK3", 
	"HERPUD1", "EIF5AL1", "EIF5A", "CCT3", "HSPE1", "GOT2", "C1QBP", "LDHB", "MAD2L1",
	"ZWINT", "ASF1B", "CDK1", "RRM2", "NCAPH", "FEN1", "TYMS", "ANLN", "HMGB2", "SCD",
	"SCD", "NDRG1","ERO1A", "NDUFA4L2", "P4HA1", "QSOX1", "BNIP3L", "FXYD3"),
		CCA5=c("CALM1", "DEGS2", "FASN", "FUT2", "MAP1LC3B", "ROIK3",
        "HERPUD1", "EIF5AL1", "EIF5A", "CCT3", "HSPE1", "GOT2", "C1QBP", "LDHB", "MAD2L1",
        "ZWINT", "ASF1B", "CDK1", "RRM2", "NCAPH", "FEN1", "TYMS", "ANLN", "HMGB2", "SCD",
        "SCD", "NDRG1","ERO1A", "NDUFA4L2", "P4HA1", "QSOX1", "BNIP3L", "FXYD3"),
                     HCC10=c("RAD51AP1", "PLK4", "AURKB", "MKI67", "EZH2", "LMNB1", "HMGB2", "KPNA2", "LRR1", "TTK", "SLC1A5", "SLC19A1", "SLC1A4", "SLCO4A1", "GPM6A", "MXD4", "AFM", "ATP5L2", "GLI4", "SLC29A4", "HEXIM1", "PRRG3", "PBXIP1", "CAPNS1", "CLDN11", "HMGCS2", "MLXIPL", "ABAT", "ADH4"),
		    HCC24=c("PLK4", "TTK", "NDC80", "MTFR2", "MAGEB1", "DIAPH3", "IQGAP3", "RAD51AP1", "TRAIP", "HMGB2", "DEPDC1B", "LRR1", "SLC19A3", "FRAS1", "MKI67", "DDIAS", "ITIH4", "IGFBP1", "ADH4", "GC", "AZGP1", "CD36", "CD44", "UPB1", "NR1H4", "HSD17B2", "TNFSF15", "TAT", "CYP2C19"),
		    HCC23=c("PLK4", "TTK", "NDC80", "MTFR2", "MAGEB1", "DIAPH3", "IQGAP3", "RAD51AP1", "TRAIP", "HMGB2", "DEPDC1B", "LRR1", "SLC19A3", "FRAS1", "MKI67", "DDIAS", "ITIH4", "IGFBP1", "ADH4", "GC", "AZGP1", "CD36", "CD44", "UPB1", "NR1H4", "HSD17B2", "TNFSF15", "TAT", "CYP2C19"),
		    HCC6=c("PLK4", "TTK", "NDC80", "MTFR2", "MAGEB1", "DIAPH3", "IQGAP3", "RAD51AP1", "TRAIP", "HMGB2", "DEPDC1B", "LRR1", "SLC19A3", "FRAS1", "MKI67", "DDIAS", "ITIH4", "IGFBP1", "ADH4", "GC", "AZGP1", "CD36", "CD44", "UPB1", "NR1H4", "HSD17B2", "TNFSF15", "TAT", "CYP2C19"))

line_specific_groups <- list(CCA1=c("Progenitor", "Differentiated1", 
			"TICs", "Differentiated2"),
			CCA5=c("Differentiated1", "Differentiated1", "TICs", "Quiescent", "Differentiated2"),
			HCC10=c("TICs", "Progenitor", "Differentiated"),
			HCC23=c("TICs", "Differentiated", "Progenitor"),
			HCC6=c("Progenitor1", "Differentiated", "Progenitor2", "TICs"),
			HCC24=c("Differentiated", "Progenitor1", "Progenitor2", "TICs")) # cluster names

line_specific_colours <- list(
	CCA1=c(colour_palette[["Prog1"]], colour_palette[["Chol1"]], 
		colour_palette[["TICs"]], colour_palette[["Chol2"]]),
	CCA5=c(colour_palette[["Chol1"]], colour_palette[["Chol1"]], 
		colour_palette[["TICs"]], "#f768a1", colour_palette[["Chol2"]]),
	HCC10=c(colour_palette[["TICs"]], colour_palette[["Prog1"]], 
		colour_palette[["Hep1"]]),
	HCC24=c(colour_palette[["Hep1"]], colour_palette[["Prog1"]], 
		colour_palette[["Prog2"]], colour_palette[["TICs"]]),
	HCC23=c(colour_palette[["TICs"]], colour_palette[["Hep1"]], 
		colour_palette[["Prog1"]]),
	HCC6=c(colour_palette[["Prog1"]], colour_palette[["Chol1"]], 
		colour_palette[["Prog2"]], colour_palette[["TICs"]])
)
	

scmap_results <- list(CCA1="CCA1_scmap_output.rds", 
			CCA5="CCA5_scmap_output.rds",
			HCC6="CCA5_scmap_output.rds",
			HCC23="HCC23_scmap_output.rds",
		      HCC10="HCC10_scmap_output.rds", 
		     HCC24="HCC24_scmap_output.rds")



for (i in names(sce_objs)) {
	this_name <- paste(i, dimred_name, "final_Bespoke", sep="_")
	print(this_name)
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
	#palette <- cluster_col(max(SCE$Manual_Clusters))
	palette <-  line_specific_colours[[i]]
	cell_colours <- palette[as.numeric(SCE$Manual_Clusters)]
	nCs <- factor_counts(SCE$Manual_Clusters)

	# save colours
	names(palette) <- 1:length(palette)
	SCE@metadata$palette <- palette
	SCE@metadata$C_keep <- nCs > 10
	SCE@metadata$C_names <- line_specific_groups[[i]]
	name_map <- SCE@metadata$C_names
	names(name_map) <- names(SCE@metadata$C_keep)[SCE@metadata$C_keep]


	if (i == "HCC10") {
		include <- SCE$Old_Manual_Clusters != 7
		SCE <- SCE[,include]
		cell_colours <- cell_colours[include]
	}

	coords <- readRDS(dim_reduction[[i]])
	
	coords <- coords[[dimred_name]]

	cell_keep <- SCE$Manual_Clusters %in% names(SCE@metadata$C_keep)[SCE@metadata$C_keep]
	SCE <- SCE[,cell_keep]
	SCE$named_clusters <- name_map[match(SCE$Manual_Clusters, names(name_map))]
	
	# ScatterPlot
	pdf(paste(this_name, "DimRedScatter.pdf", sep="_"), width=6, height=6)
	plot(coords$x[cell_keep], coords$y[cell_keep], col=cell_colours, 
		pch=16, xlab="DM_1", ylab="DM_2")
	dev.off()
	blank_plot <- function() {
		tmp <- par("mar")
		par(mar=c(0,0,0,0))
		plot(1,1, col="white", xlim=c(0,1), ylim=c(0,1), xaxt="n", yaxt="n", main="", xlab="", ylab="", bty="n")
		par(mar=tmp)
	}
	pdf(paste(this_name, "DimRedScatter_Legend.pdf", sep="_"), width=6, height=6)
	blank_plot();
	my_order <- order(name_map, decreasing=T)
	legend("center", name_map[my_order], col=palette[my_order], pch=16, cex=2, bty="n")
	dev.off()


	# Marker - DotPlot
	require("ggplot2")
	seurat <- as.Seurat(SCE, data=expr_type)
	thing <- cbind(coords$x, coords$y)
	thing <- apply(thing, 2, function(x){
		x=x-min(x); x<-x/max(x)*2; x<-x-1})
	rownames(thing) <- colnames(SCE)
	colnames(thing) <- c("DM1", "DM2")
	seurat[["dm"]] <- CreateDimReducObject(embeddings = thing, key = "DM", assay = DefaultAssay(seurat))
	

	pdf(paste(this_name, "MarkerDot.pdf", sep="_"), width=8, height=5)
	a<-DotPlot(seurat, 
		features=line_specific_genes[[i]], 
		group.by="named_clusters")+
		theme(axis.text.x = element_text(angle = 90, hjust = 1))
	plot(a+labs(x="Genes", y="Type"))
	dev.off()
	# Marker - Scatters
	#pdf(paste(this_name, "MarkerScatter.pdf", sep="_"), width=16, height=14)
	#FeaturePlot(seurat, reduction="dm", features = line_specific_genes[[i]])
	#dev.off()

	get_top_markers <- function(cluster, nmarks=5) {
		tmp <- FindMarkers(seurat, ident.1=cluster, group.by="named_clusters")
		tmp$detect_diff <- tmp$pct.1-tmp$pct.2
		tmp <- tmp[tmp$avg_logFC > 0 & tmp$detect_diff > 0,]
		return(rownames(tmp)[1:nmarks])
	}
	#heatmap_genes <- c();
	#for (n in unique(seurat@meta.data$named_clusters)) {
	#	heatmap_genes <- c(heatmap_genes, get_top_markers(n))
	#}

	#pdf(paste(this_name, "MarkerHeatmap.pdf", sep="_"), width=14, height=8)
	#seurat <- ScaleData(seurat)
	#my_order <- order(name_map, decreasing=F)
	#DoHeatmap(seurat, features = c(heatmap_genes,line_specific_genes[[i]]),
	#	 size = 3, group.by = "named_clusters", group.colors=palette[my_order])
	#dev.off()
	
	dat<- seurat@assays$RNA@data[rownames(seurat) %in% line_specific_genes[[i]],]
	dat <- data.frame(t(dat), seurat@meta.data$named_clusters)
	colnames(dat)[ncol(dat)] <- "type"
	
	# Marker Violin
	ggplot_palette <- palette
	names(ggplot_palette) <- name_map
	marker_violins <- list();
	for (gene in line_specific_genes[[i]]) {
		pdf(paste(this_name, gene, "MarkerScatter.pdf", sep="_"), width=8, height=8)
		print(FeaturePlot(seurat, reduction="dm", features = gene))
		dev.off()
	}
	#a <- ggplot(dat, aes_string(x="type", y=gene, fill="type"))+geom_violin()+scale_fill_manual(values=ggplot_palette)+ ggtitle(gene)+theme(plot.title = element_text(face="bold", size=30,hjust = 0.5 ), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
	#panel.background = element_blank(), axis.line = element_line(colour = "black"))
	#marker_violins[[gene]] <- a;
	#}
	#require(ggpubr)
	#theme_set(theme_pubr())
	#pdf(paste(this_name, "MarkerViolins.pdf", sep="_"), width=16, height=16)
	#ggarrange(plotlist=marker_violins)
	#dev.off()

	# ScatterPlot + scmap results
	scmap_out <- readRDS(scmap_results[[i]])
	pdf(paste(this_name, "ScmapScatter.pdf", sep="_"), width=6, height=6)
	plot(coords$x, coords$y, col="black", bg=scmap_out$scmap_cell_Cols, pch=21, xlab="DM_1", ylab="DM_2");
	lgend <- unique(scmap_out$scmap_cluster_labs);
	lgend_col <- unique(scmap_out$scmap_cell_Cols);
	dev.off()

	pdf(paste(this_name, "Scmap_Legend.pdf", sep="_"), width=6, height=6)
	blank_plot();
	legend("center",lgend, col=lgend_col, pch=16, cex=2, bty="n")
	dev.off()
}
