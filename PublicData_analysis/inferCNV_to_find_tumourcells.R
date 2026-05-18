
run_seurat_pipeline <- function(seur_obj) {
	#quickly cluster the object
	seur_obj <- Seurat::NormalizeData(seur_obj);
	seur_obj <- Seurat::ScaleData(seur_obj);
	seur_obj <- Seurat::FindVariableFeatures(seur_obj, selection.method = "vst", nfeatures = 2000)
	hvgs <- VariableFeatures(seur_obj); hvgs <- hvgs[!grepl("^MT-", hvgs)];
	VariableFeatures(seur_obj) <- hvgs;
	seur_obj <- Seurat::RunPCA(seur_obj, features = VariableFeatures(object = seur_obj))
	seur_obj <- Seurat::FindNeighbors(seur_obj, dims = 1:20)
	seur_obj <- Seurat::FindClusters(seur_obj)
	s.genes <- cc.genes$s.genes
	g2m.genes <- cc.genes$g2m.genes
	seur_obj <- Seurat::CellCycleScoring(seur_obj, s.features = s.genes, g2m.features=g2m.genes, set.ident=TRUE)
	return(seur_obj)
}

run_infercnv <- function(seur_obj, anno_file, out_tag) {
	gene_pos <- "/home/tandrew6/projects/def-tandrew6/tandrew6/External_Data/Ensembl/EnsemblBiomart_09Dec2020_HNGC_Gene_positions.txt"
	require(infercnv)
	infercnv_obj = CreateInfercnvObject(seur_obj@assays$RNA@counts,
                                    annotations_file=anno_file,
                                    delim="\t", 
				    gene_order_file=gene_pos,
				    ref_group_names=NULL,
				    chr_exclude=c("chrMT", "chrX", "chrY", "chrM"))
	# perform infercnv operations to reveal cnv signal

	infercnv_obj = infercnv::run(infercnv_obj,
				cutoff=0.1,  # use 1 for smart-seq, 0.1 for 10x-genomics
				out_dir=paste(out_tag, "infercnv_output", sep="_"),  # dir is auto-created for storing outputs
                               cluster_by_groups=T,   # cluster
                               denoise=T,
                               HMM=F
	                      )

	saveRDS(infercnv_obj, file=paste(out_tag, "infercnv_obj.rds", sep="_"))
}

create_anno_file <- function(seur_obj, out_tag) {
	anno_file <- paste(out_tag, "cell_annotation.txt", sep="_")
	tab <- cbind(rownames(seur_obj@meta.data), seur_obj@meta.data$inferCNV_groups)
	tmp <- table(tab[,2])
	exclude <- names(tmp)[tmp<10]
	if (length(exclude) > 0) {
		seur_obj <- seur_obj[,! (seur_obj@meta.data$inferCNV_groups %in% exclude)]
	}
	tab <- tab[tab[,1] %in% rownames(seur_obj@meta.data),]
	write.table(tab, anno_file, col.names=F, row.names=F, sep="\t", quote=F)
	return(list(obj=seur_obj, anno_file=anno_file))
}

##### MA Tumours #####

source("/home/tandrew6/projects/def-tandrew6/tandrew6/External_Data/Merixtell/MerixtellLoadData.R")

#seur_obj <- load_MaTumors(min.cells=20, min.features=500)
#out_tag <- "MaTumors"

#seur_obj <- run_seurat_pipeline(seur_obj);

#Cell Annotations
# Assuming that non-tumour cells and tumour cells are present in each sample, and non-tumour cells are likely to co-cluster between samples, then I want to intersect the clusters above with the sampleID and run inferCNV on those groups.
#seur_obj@meta.data$inferCNV_groups <- paste(seur_obj@meta.data$Sample, seur_obj@meta.data$seurat_clusters, sep="-")

#stuff <- create_anno_file(seur_obj, out_tag)
#seur_obj <- stuff$obj
#anno_file <- stuff$anno_file


#run_infercnv(seur_obj, anno_file, out_tag)


##### ZhangCC #####
ZhangMats <- load_ZhangCC_aslist()
out_tag <- "ZhangCC"

seur_obj <- convert_MatrixList2Seurat(ZhangMats, min.cells=20, min.features=500)
seur_obj <- run_seurat_pipeline(seur_obj)

anno_file <- paste(out_tag, "cell_annotation.txt", sep="_")
seur_obj@meta.data$inferCNV_groups <- paste(seur_obj@meta.data$origin, seur_obj@meta.data$seurat_clusters, sep="-")

stuff <- create_anno_file(seur_obj, out_tag)
seur_obj <- stuff$obj
anno_file <- stuff$anno_file
run_infercnv(seur_obj, anno_file, out_tag)


##### SuHCC #####
#seur_obj <- load_SuHCC(min.cells=20, min.features=500)
#out_tag <- "SuHCC"

#seur_obj <- run_seurat_pipeline(seur_obj);

#anno_file <- paste(out_tag, "cell_annotation.txt", sep="_")
#seur_obj@meta.data$inferCNV_groups <- paste(seur_obj@meta.data$patient, seur_obj@meta.data$seurat_clusters, sep="-")

#stuff <- create_anno_file(seur_obj, out_tag)
#seur_obj <- stuff$obj
#anno_file <- stuff$anno_file
#run_infercnv(seur_obj, anno_file, out_tag)
