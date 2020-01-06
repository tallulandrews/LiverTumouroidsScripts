# Lineage markers
Chol_lineage <- read.table("/nfs/users/nfs_t/ta6/Collaborations/LiverOrganoids/Markers_130418_Chol.txt", header=TRUE)
Hep_lineage <- read.table("/nfs/users/nfs_t/ta6/Collaborations/LiverOrganoids/Markers_130418_Hep.txt", header=TRUE)

Hep_both <- Hep_lineage[ grepl("Prog", Hep_lineage[,2]) & grepl("Hep", Hep_lineage[,2]), 1]
Chol_both <- Chol_lineage[ grepl("Prog", Chol_lineage[,2]) & grepl("Chol", Chol_lineage[,2]), 1]
Prog_both <- Hep_lineage[ grepl("Prog", Hep_lineage[,2]) & Hep_lineage[,1] %in% Chol_lineage[Chol_lineage[,2] == "Prog",1], 1]

Conflict1 <- Hep_lineage[ grepl("Hep", Hep_lineage[,2]) & Hep_lineage[,1] %in% Chol_lineage[Chol_lineage[,2] == "Chol",1], 1]
Conflict2 <- Hep_lineage[ grepl("Prog", Hep_lineage[,2]) & Hep_lineage[,1] %in% Chol_lineage[Chol_lineage[,2] == "Chol",1], 1]
Conflict3 <- Hep_lineage[ grepl("Hep", Hep_lineage[,2]) & Hep_lineage[,1] %in% Chol_lineage[Chol_lineage[,2] == "Prog",1], 1]
Conflict4 <- Chol_lineage[ grepl("Prog", Chol_lineage[,2]) & Chol_lineage[,1] %in% Hep_lineage[Hep_lineage[,2] == "Hep",1], 1]

Conflicts <- c(as.character(Conflict1), as.character(Conflict2), as.character(Conflict3), as.character(Conflicts4))

Chol_lineage <- Chol_lineage[Chol_lineage[,1] %in% marker_genes & Chol_lineage[,1] %in% keep_genes,]
Hep_lineage <- Hep_lineage[Hep_lineage[,1] %in% marker_genes & Hep_lineage[,1] %in% keep_genes,]

Chol_lineage[,2] <- as.character(Chol_lineage[,2])
Chol_lineage[Chol_lineage[,2] == "Prog",2] <- "Chol-Prog"
Chol_lineage[Chol_lineage[,2] == "Chol",2] <- "Chol-Mature"
Chol_lineage[Chol_lineage[,1] %in% Chol_both,2] <- "Chol-Both"
Chol_lineage <- Chol_lineage[!(Chol_lineage[,1] %in% Conflicts),]

Hep_lineage[,2] <- as.character(Hep_lineage[,2])
Hep_lineage[Hep_lineage[,2] == "Prog",2] <- "Hep-Prog"
Hep_lineage[Hep_lineage[,2] == "Hep",2] <- "Hep-Mature"
Hep_lineage[Hep_lineage[,1] %in% Hep_both,2] <- "Hep-Both"
Hep_lineage <- Hep_lineage[!(Hep_lineage[,1] %in% Conflicts),]

Lineage <- rbind(Chol_lineage,Hep_lineage)
Lineage[Lineage[,1] %in% Prog_both,2] <- "Common-Prog"
Lineage <- Lineage[!duplicated(Lineage[,1]),]
Lineage[,1] <- as.character(Lineage[,1])
Lineage[Lineage[,1] == "05-Mar",1] <- "MARCH5"
Lineage <- unique(Lineage)

write.table(Lineage, file="Cleaned_Lineage.txt", row.names=F, col.names=F)

