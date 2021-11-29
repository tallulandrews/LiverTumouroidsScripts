source("/nfs/users/nfs_t/ta6/My_R_packages/TreeOfCells/R/Fitting_ZINB.R")
require("Matrix")
require("scater")
require("methods")

args <- commandArgs(trailingOnly=TRUE)
f <- args[1]
#for (f in files) {
#print(f)

sce <- readRDS(f)

outfile <- sub("^.*/", "", f)
outfile <- sub(".rds", "_profiles.rds", outfile);

type_labs <- sce$cell_type1;
if (length(type_labs) != ncol(sce)) {
	type_col <- grep("celltype", colnames(colData(sce)), ignore.case=TRUE);
	type_labs <- colData(sce)[,type_col]
}

OUT_list <- list()
i = 0
for (type in levels(factor(type_labs))) {
i <- i+1
if (i <= 21) {next;}
	out <- fit_ZINB_to_matrix(sce@assays[["counts"]][,type_labs == type & !is.na(type_labs)]);
	OUT_list[[paste(sce$tissue[1], type)]] <- out;
if ( i %% 3 == 0) {
saveRDS(OUT_list, file=sub(".rds", paste(i,".rds", sep=""), outfile))


}

}
saveRDS(OUT_list, file=outfile)

#}
