data <- read.table("E-MTAB-3929_counts.txt", header=TRUE)
lab <- strsplit(colnames(data), "\\.");

get_metadata <- function(x){ tmp <- unlist(x); c(x[1], x[2])}

labels <- matrix(unlist(lapply(lab, get_metadata)), ncol=2, byrow=T)

cell_type <- rep("blast", times=length(labels[,1]))
cell_type[labels[,1] == "E4"] = "16cell"
cell_type[labels[,1] == "E3"] = "8cell"


ANN <- data.frame(Species = rep("Homo sapiens", times=length(labels[,1])), cell_type1 = cell_type, Source = rep("preimplantation embryo", times=length(labels[,1])), age = labels[,1], batch=labels[,2])
rownames(ANN) <- colnames(data);

require("scater")
pd <- new("AnnotatedDataFrame", data=ANN)
petro <- newSCESet(countData=as.matrix(data), phenoData=pd)
fData(petro)$feature_symbol <- featureNames(petro)

saveRDS(petro, "petropoulos.rds")
