Meri_data_dir <- "/home/tandrew6/projects/def-tandrew6/tandrew6/External_Data/Merixtell"

load_MaTumors <- function(min.cells=1, min.features=10) {
	require(Seurat)
	require(Matrix)
	obj1 <- Matrix::readMM(paste(Meri_data_dir,"GSE125449_Set1_matrix.mtx.gz", sep="/"))
	rownames(obj1) <- read.delim(paste(Meri_data_dir, "GSE125449_Set1_genes.tsv.gz", sep="/"), sep="\t", header=FALSE)[,2]
	colnames(obj1) <- read.delim(paste(Meri_data_dir, "GSE125449_Set1_barcodes.tsv.gz", sep="/"), sep="\t", header=FALSE)[,1]


	obj2 <- Matrix::readMM(paste(Meri_data_dir,"GSE125449_Set2_matrix.mtx.gz", sep="/"))
	rownames(obj2) <- read.delim(paste(Meri_data_dir, "GSE125449_Set2_genes.tsv.gz", sep="/"), sep="\t", header=FALSE)[,2]
	colnames(obj2) <- read.delim(paste(Meri_data_dir, "GSE125449_Set2_barcodes.tsv.gz", sep="/"), sep="\t", header=FALSE)[,1]

	metadata1 <- read.table(paste(Meri_data_dir, "GSE125449_Set1_samples.txt.gz", sep="/"), sep="\t", header=TRUE)
	metadata2 <- read.table(paste(Meri_data_dir, "GSE125449_Set2_samples.txt.gz", sep="/"), sep="\t", header=TRUE)
	origin <- read.delim(paste(Meri_data_dir, "GSE125449_manual_origin_annotation.txt", sep="/"), sep="\t", header=TRUE)
	MetaData <- rbind(metadata1, metadata2)
	MetaData$tumour_type <- unlist(origin[MetaData[,1]])


	common_genes <- sort(intersect(rownames(obj1), rownames(obj2)))

	obj1 <- obj1[match(common_genes, rownames(obj1)),]
	obj2 <- obj2[match(common_genes, rownames(obj2)),]
	count_mat <- cbind(obj1, obj2)
	if (!identical(MetaData[,2], colnames(count_mat))) {stop("Error: Metadata doesn't match count matrix!")}
	rownames(MetaData) <- colnames(count_mat)
	seur_obj <- CreateSeuratObject(count_mat, meta.data=MetaData, project="Ma_tumors", min.cells=min.cells, min.features=min.features)
	return(seur_obj)
}

load_ZhangCC_aslist <- function() {
	require(Matrix)
	files <- c("GSM4116580_ICC_18_Tumor_UMI.csv.gz", "GSM4116581_ICC_20_Tumor_UMI.csv.gz", "GSM4116583_ICC_23_Tumor_UMI.csv.gz", "GSM4116584_ICC_24_Tumor1_UMI.csv.gz", "GSM4116585_ICC_24_Tumor2_UMI.csv.gz", "GSM4116579_ICC_18_Adjacent_UMI.csv.gz", "GSM4116582_ICC_23_Adjacent_UMI.csv.gz", "GSM4116586_ICC_25_Adjacent_UMI.csv.gz")
	count_mats <- list();
	common_genes <- NULL;
	for (f in files) {
		count_mat <- read.table(paste(Meri_data_dir, f, sep="/"), sep=",", header=TRUE, stringsAsFactors=FALSE)
		rownames(count_mat) <- count_mat[,1];
		count_mat <- Matrix::Matrix(as.matrix(count_mat[,-1]))
		if (is.null(common_genes)) {
			common_genes <- rownames(count_mat)
		} else {
			common_genes <- sort(intersect(common_genes, rownames(count_mat)))
		}
		count_mats[[f]] <- count_mat
	}
	count_mats[["common_genes"]] <- common_genes
	return(count_mats)
}

load_SuHCC <- function(min.cells=1, min.features=10) {
	require(Matrix)
	require(Seurat)
	count_mat <- read.delim(paste(Meri_data_dir, "GSE146115_HCC1-2-5-9_count_with_ERCC.txt.gz", sep="/"), sep="\t", header=TRUE, stringsAsFactors=FALSE)
	count_mat <- count_mat[!grepl("^ERCC-", count_mat[,1]),]
	genes <- fix_Excel_gene_names(count_mat[,1])
	count_mat <- Matrix::Matrix(as.matrix(count_mat[,-1]))
	count_mat <- count_mat[!duplicated(genes),]
	

	rownames(count_mat) <- genes[!duplicated(genes)]
	metadata <- read.delim(paste(Meri_data_dir, "GSE146115_raw_to_prcessed_data_column.txt.gz", sep="/"), sep="\t")
	colnames(metadata) <- metadata[1,]
	metadata <- metadata[-1,]
	metadata <- metadata[grepl("^Sample", metadata[,1]),]
	sample_ids <- metadata[,3]
	sample_ids <- sub("L0", "L", sample_ids)
	sample_ids <- sub("W0", "W", sample_ids)
	rownames(metadata) <- sample_ids
	metadata <- metadata[match(colnames(count_mat), sample_ids),]
	colnames(metadata) <- c("Sample", "name", "cellID", "rawFile", "md5checksum")
	metadata$patient <- sapply(strsplit(metadata$name, ","), "[", 1)
	metadata <- metadata[,-c(5,6)]

	seur_obj <- CreateSeuratObject(count_mat, meta.data=metadata, project="Su_HCC", min.cells=min.cells, min.features=min.features)
	return(seur_obj)

}

fix_Excel_gene_names <- function(genes, species=c("Human")) {
	if (! species %in% c("Human")) { stop(paste("Error:", species, "currently not supported.")) }
	MARCH_nums <- sapply(strsplit(genes[grep("-Mar", genes)], "-"), "[", 1)
	MARCH_HGNC <- paste("MARCHF", MARCH_nums, sep="")
	genes[grep("-Mar", genes)] <- MARCH_HGNC

	SEPT_nums <- sapply(strsplit(genes[grep("-Sep", genes)], "-"), "[", 1)
	SEPT_HGNC <- paste("SEPTIN", SEPT_nums, sep="")
	genes[grep("-Sep", genes)] <- SEPT_HGNC

	DEC_nums <- sapply(strsplit(genes[grep("-Dec", genes)], "-"), "[", 1)
	DEC_HGNC <- paste("DEC", DEC_nums, sep="")
	genes[grep("-Dec", genes)] <- DEC_HGNC

	return(genes)
}

# Takes a list of count matricies and merges them together, keeping only genes that appear in all count matricies, and creates a seurat object from it.
convert_MatrixList2Seurat <- function(mat_list, common_genes=NULL, project="project", min.cells=1, min.features=10) {
	require(Seurat)
	if (is.null(common_genes)) {
		if (is.null(mat_list[["common_genes"]])) {
			common_genes<- rownames(mat_list[[1]]);
			for (m in mat_list) {
				common_genes <- intersect(common_genes, rownames(m))
			}
			common_genes <- sort(common_genes);
		} else {
			common_genes <- mat_list[["common_genes"]]
			mat_list[["common_genes"]] <- NULL
		}
	}
	if (length(common_genes) < 2) {stop("Error: Less than 2 common genes between matrices.")}
	CountMat <- NULL
	origin <- c()
	for (i in 1:length(mat_list)) {
		m <- mat_list[[i]]
		m <- m[match(common_genes, rownames(m)),]
		colnames(m) <- paste(colnames(m), i, sep="-")
		origin <- c(origin, rep(names(mat_list)[i], ncol(m)))
		if (is.null(CountMat)) {
			CountMat <- m
		} else {
			CountMat <- cbind(CountMat, m)
		}
	}
	MetaData <- data.frame(origin=origin); rownames(MetaData) <- colnames(CountMat);
	seur_obj <- CreateSeuratObject(CountMat, meta.data=MetaData, project=project, min.cells=min.cells, min.features=min.features)
	return(seur_obj)
}
