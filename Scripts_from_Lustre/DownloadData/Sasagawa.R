dat1 <- read.table("Sasagawa_ExprMat1.txt", header=T)
ann1 <- read.table("Sasagawa_Ann1.txt", header=T, stringsAsFactors=FALSE)
keep <- ann1[4,] == "single cell" & ann1[7,] == "Quartz-Seq" & ann1[6,] != "All phase"
ann1 <- ann1[,keep]
ann1[6,ann1[6,]=="G2/M"] <- "G2M"
dat1 <- dat1[,keep]
ANN1 <- data.frame(Species=unlist(ann1[2,]), 
		cell_type2=unlist(ann1[4,]), 
		cell_type1=unlist(ann1[6,]))
rownames(ANN1) <- ann1[8,]

# Not single cell
#dat2 <- read.table("Sasagawa_ExprMat2.txt", header=T)
#ann2 <- read.table("Sasagawa_Ann2.txt", header=T, stringsAsFactors=FALSE)
#keep <- ann2[4,] == "single cell" & ann2[7,] == "Quartz-Seq" & ann2[6,] != "All phase"
#ann2 <- ann2[,keep]
#ann2[6,ann2[6,]=="G2/M"] <- "G2M"
#dat2 <- dat2[,keep]
#ANN2 <- data.frame(Species=unlist(ann2[2,]), 
#		cell_type2=unlist(ann2[4,]), 
#		cell_type1=unlist(ann2[6,]))
#rownames(ANN2) <- ann2[8,]

fpkms <- dat1
ann <-ANN1

source("~/R-Scripts/Ensembl_Stuff.R")
gene_symbol <- General_Map(rownames(fpkms), in.org="Mmus", in.name="ensg", out.org="Mmus", out.name="symbol")

require("scater")
SCE <- SingleCellExperiment(assays=list(fpkms=as.matrix(fpkms)), colData=ann)
rowData(SCE)$feature_symbol <- gene_symbol


saveRDS(SCE, file="Sasagawa.rds")
