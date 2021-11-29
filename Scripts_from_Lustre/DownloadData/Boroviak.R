# Boroviak et al. (2015) Lineage-Specific Profiling Delineates the Emergence and Progression of Naive Pluripotency in Mammalian Embryogenesis. Developmental Cell. 35(3):366-382
# PMID: 26555056
# ArrayExpress: E-MTAB-2958, E-MTAB-2959

data <- read.table("mmc2.csv", stringsAsFactors=FALSE, header=T, sep="\t")

rownames(data) <- data[,3]
data <- data[,-c(1:3)]
labs <- as.character(colnames(data))
time <- labs

time[grep("E2.5", time)] <- "E2.5"
time[grep("E3.5", time)] <- "E3.5"
time[grep("E4.5", time)] <- "E4.5"
time[grep("E5.5", time)] <- "E5.5"
time[grep("E5.5", time)] <- "E5.5"
time[grep("DIA", time)] <- "DIA"
time[grep("ESC", time)] <- "2iL"

labs[grep("MOR", labs)] <- "morula"
labs[grep("EPI", labs)] <- "EPI"
labs[grep("PrE", labs)] <- "PE"
labs[grep("ICM", labs)] <- "ICM"
labs[grep("ESC", labs)] <- "ESC"

require("scater")
P <- data.frame(cell_type1=labs, Time=time); rownames(P) <- colnames(data);
pd <- new("AnnotatedDataFrame", data=P)
BORO <- newSCESet(fpkmData=as.matrix(data), phenoData=pd, logExprsOffset=1)
fData(BORO)$feature_symbol <- rownames(exprs(BORO));
saveRDS(BORO, "Boroviak.rds")

