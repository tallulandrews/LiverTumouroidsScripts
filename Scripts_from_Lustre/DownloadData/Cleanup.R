files <- Sys.glob("*10X.rds");
for ( f in files ) {
	require("scater")
	require("SingleCellExperiment")
	dat <- readRDS(f);
	exclude <- is.na(dat$cell_type1) | dat$cell_type1 == "unknown"
	if (sum(exclude) > 0 & sum(exclude) < length(exclude)) {
		dat <- dat[,!exclude]
		saveRDS(dat, file=f);
	}
	rm(dat)
}
