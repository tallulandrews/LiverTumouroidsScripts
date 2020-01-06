# Lineage markers
Lineage <- read.table("/nfs/users/nfs_t/ta6/Collaborations/LiverOrganoids/PatsGeneSets.txt", header=FALSE)

hep_mature <- Lineage[Lineage[,2]=="Hep-Mature",1]
Conflict1 <- hep_mature[hep_mature %in% Lineage[grep("Chol",Lineage[,2]),1]]
chol_mature <- Lineage[Lineage[,2]=="Chol-Mature",1]
Conflict2 <- chol_mature[chol_mature %in% Lineage[grep("Hep",Lineage[,2]),1]]

Lineage <- Lineage[! Lineage[,1] %in% c(as.character(Conflict1), as.character(Conflict2)), ]

write.table(Lineage, file="Cleaned_PatLineage.txt", row.names=F, col.names=F)

