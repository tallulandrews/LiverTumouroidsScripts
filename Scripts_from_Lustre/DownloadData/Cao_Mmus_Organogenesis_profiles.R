#source("/nfs/users/nfs_t/ta6/My_R_packages/TreeOfCells/R/Fitting_ZINB.R")
require("TreeOfCells")
require("Matrix")
require("methods")

require(SingleCellExperiment)

infile <- commandArgs(trailingOnly=T)

#sce <- readRDS("Cao_Mmus_Organogenesis_p1.rds");
sce <- readRDS(infile);

require(stringr)
outfile <- str_replace(infile, ".rds", "")

set <- as.numeric(commandArgs(trailingOnly=TRUE));

OUT_list <- list();
for (type in levels(factor(sce$cell_type1))[set]) {
        print(type)
        if (sum(sce$cell_type1==type) < 20) {next;}
        out <- fit_NB_to_matrix(as.matrix(assays(sce)$counts[,sce$cell_type1==type]));
        OUT_list[[type]] <- out;
}
saveRDS(OUT_list, file=paste(outfile,"_NB_",set,".rds", sep=""));

