# Human blastocyst
# Blakeley et al. (2015) Defining the three cell lineages of the human blastocyst by single-cell RNA-seq. Development. 142(18):3151-65
# GEO: GSE66507
# PMID: 26293300


x <- read.table("GSE66507_human_blastocyst_rnaseq_counts.txt.gz", header=TRUE)

rownames(x) <- x[,1];
x <- x[,-1];
type <- colnames(x);
type[grep("TE", type)] <- "TE"
type[grep("PE", type)] <- "PE"
type[grep("EPI", type)] <- "EPI"



require("scater")
P <- data.frame(cell_type1=type); rownames(P) <- colnames(x);
pd <- new("AnnotatedDataFrame", data=P)
blakeley <- newSCESet(countData=as.matrix(x), phenoData=pd)
blakeley <- getBMFeatureAnnos(
    blakeley, attributes=c("ensembl_transcript_id", "ensembl_gene_id","chromosome_name", "transcript_biotype", "transcript_start", "transcript_end"), filters="ensembl_gene_id",
    biomart="ensembl", dataset="hsapiens_gene_ensembl")

saveRDS(blakeley, "blakeley.rds")

