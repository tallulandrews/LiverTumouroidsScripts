# Reference: Scialdone et al. (2016). Resolving early mesoderm diversitification through single-cell expression profiling

x1 = read.table("Gastrulation_counts.txt")
m1 = read.table("Gastrulation_metadata.txt", header=T, stringsAsFactor=FALSE)
rownames(m1) = m1[,1]
m1 <- m1[,2:5]
colnames(m1) <- c("location","individual","Sorting","cell_type1")
m1$Genotype <- rep("WT", times=length(m1[,1]))


x2 = read.table("Gastrulation_Tal1counts.txt")
m2 = read.table("Gastrulation_metadataTal1.txt", header=T, stringsAsFactor=FALSE)
rownames(m2) = m2[,1]
m2 <- m2[,2:5]
colnames(m2) <- c("location","individual","Sorting","cell_type1")
m2$Genotype <- rep("WT", times=length(m2[,1]))
m2[grep("Null", m2$Sorting),5] = "Tal1-KO"

M <- rbind(m1,m2)
M[M[,4] == "turquoise",4] <- "epiblast"
M[M[,4] == "pink",4] <- "extraembryonic ectoderm"
M[M[,4] == "magenta",4] <- "visceral endoderm"
M[M[,4] == "blue",4] <- "nascent mesoderm"
M[M[,4] == "darkorange",4] <- "posterior mesoderm"
M[M[,4] == "red",4] <- "endothelium"
M[M[,4] == "yellow",4] <- "blood progenitors"
M[M[,4] == "brown",4] <- "embryonic blood"
M[M[,4] == "black",4] <- "allantois"
M[M[,4] == "green",4] <- "pharyngeal mesoderm"
M[M[,4] == "grey",4] <- "unknown (mesoderm)"

Counts <- cbind(x1,x2)

require("scater")
pd <- new("AnnotatedDataFrame", data=M)
scialdone <- newSCESet(countData=as.matrix(Counts), phenoData=pd)
scialdone <- getBMFeatureAnnos(
    scialdone, filters="ensembl_gene_id",
    biomart="ensembl", dataset="mmusculus_gene_ensembl")

saveRDS(scialdone, "scialdone.rds")
