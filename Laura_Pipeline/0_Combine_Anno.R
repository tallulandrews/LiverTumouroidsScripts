args <- commandArgs(trailingOnly=TRUE);

CtW <- read.delim(args[1], sep="\t", header=F)
Plate <- read.delim(args[2], sep="\t", header=F)


CtW <- CtW[CtW[,4] != 0, ]
keep <- CtW[ CtW[,2] %in% Plate[,1] , 2]
keep <- sort(keep)
CtW <- CtW[match(keep, CtW[,2]),] 
Plate <- Plate[match(keep, Plate[,1]),] 

stuff <- matrix(unlist(strsplit(as.character(CtW[,2]), "_")), ncol=2, byrow=TRUE)
Ann <- cbind(as.character(CtW[,1]), stuff[,1], stuff[,2], as.character(CtW[,4]), as.character(Plate[,2]))

colnames(Ann) <- c("CellID", "Plate", "Well", "File", "Type")

write.table(Ann, file=args[3], quote=TRUE, row.names=FALSE, col.names=TRUE, sep="\t")
