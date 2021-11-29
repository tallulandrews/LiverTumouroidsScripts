source("/nfs/users/nfs_t/ta6/My_R_packages/TreeOfCells/R/Fitting_ZINB.R")
require("Matrix")
require("methods")

# Cao et al. The dynamic transcriptional landscape of mammalian organogenesis at single cell resolution. Nature. 2019. 566(7745):496-502
# PMID:30787437
# GEO: GSE119945

require("Matrix");
c_data <- read.table("GSE119945_cell_annotate.csv.gz", sep=",", header=T)
g_data <- read.table("GSE119945_gene_annotate.csv.gz", sep=",", header=T)
counts <- readMM("GSE119945_gene_count.txt")

labs <- c_data$Main_Cluster
lab_names <- c("1"="Connective tissue progenitors",
		"2"="Chondrocytes/osteoblasts",
		"3"="Intermediate Mesoderm",
		"4"="Jaw/tooth progenitors",
		"5"="Excitatory neurons",
		"6"="Epithelial cells",
		"7"="Radial glia",
		"8"="Early mesenchyme",
		"9"="Neural progenitor cells",
		"10"="Postmitotic premature neurons",
		"11"="Oligodendrogyte progenitors",
		"12"="Isthmic organizer cells",
		"13"="Myocytes",
		"14"="Dorsal neural tube cells",
		"15"="Inhibitory neurons", 
		"16"="Stromal cells",
		"17"="Osteoblasts",
		"18"="Inhibitory neuron progenitors",
		"19"="Premature oligodendrocyte",
		"20"="Endothelial cells",
		"21"="Chondrocyte progenitors",
		"22"="Definitive erythroid lineage",
		"23"="Schwann cell precursor",
		"24"="Sensory neurons",
		"25"="Limb mesenchyme",
		"26"="Primitive erythroid lineage",
		"27"="Inhibitory interneurons",
		"28"="Granule neurons",
		"29"="Hepatocytes",
		"30"="Notochord and floor plate cells",
		"31"="White blood cells", 
		"32"="Ependymal cell",
		"33"="Cholinergic neurons",
		"34"="Cardiac muscle lineages",
		"35"="Megakaryocytes",
		"36"="Melanocytes",
		"37"="Lens",
		"38"="Neutrophils")
labs <- lab_names[labs]
c_data$cell_type1 <- labs;


#gene_names <- as.character(g_data[,1]);
#gene_names <- unlist(lapply(strsplit(gene_names, "\\."), function(x){return(x[[1]])}))
gene_names <- as.character(g_data$gene_short_name)
tidy <- !duplicated(gene_names);

rownames(counts) <- gene_names
colnames(counts) <- c_data[,1]
g_data <- g_data[tidy,]
rownames(g_data) <- gene_names[tidy];
tidy2 <- !is.na(labs);
rownames(c_data) <- c_data[,1]
c_data <- c_data[tidy2,]
counts <- coundt[tidy, tidy2];

require(SingleCellExperiment)
sce <- SingleCellExperiment(assays=list(counts=counts), colData=c_data, rowData=g_data);
saveRDS(sce, file="Cao_Mmus_Organogenesis.rds")

OUT_list <- list();
for (type in levels(factor(sce$cell_type1))[1:10]) {
	print(type)
	if (sum(sce$cell_type1==type) < 20) {next;}
	out <- fit_ZINB_to_matrix(as.matrix(assays(sce)$counts[,sce$cell_type1==type]));
	OUT_list[[type]] <- out;
}
saveRDS(OUT_list, file="Cao_Mmus_Organogenesis_Profiles_1_10.rds");
