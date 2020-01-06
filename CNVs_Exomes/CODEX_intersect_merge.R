
intersect_merge <- function(reg_e, reg_l, threshold=1, size_t=500000){
	both <- intersect(reg_e, reg_l)
	both <- both[wdith(both) > size_t];
	out <- as.data.frame(both)
	if (nrow(out)==0) {return(rep(NA, 5))}

	region_width <- width(both);
	gap <- gaps(both); 
	gap <- gap[start(gap) > 1]; 
	gap <- gap[end(gap) < max(end(both))]
	gap_width <- width(gap); 
	lens <- apply(cbind(region_width[1:(length(both)-1)], region_width[2:length(both)]), 1, min)

	out[,2] <- as.numeric(out[,2])
	out[,3] <- as.numeric(out[,3])
	out[,4] <- as.numeric(out[,4])
	out <- out[out[,4] > size_t,]
	if (nrow(out)==0) {return(rep(NA, 5))}
	lens <- cbind(out[1:(nrow(out)-1),4], out[2:nrow(out),4])
	lens <- apply(lens, 1, min);
	score <- (out[2:nrow(out),2]-out[1:(nrow(out)-1),3])/lens
	gaps <- which(score > threshold)
	gaps <- c(1, gaps+1)

	if (length(gaps) == 1) {
		out[1,3] <- out[nrow(out),3]
		out[1,4] <- out[1,3]-out[1,2]
		out <- out[1,]
	} else {
		for (i in 2:length(gaps)) {
			out[gaps[i-1], 3] <- out[max(1,gaps[i]-1),3]
			out[gaps[i-1], 4] <- out[gaps[i-1], 3] - out[gaps[i-1], 2]
			if (gaps[i] > gaps[i-1]+1) {
				out[(gaps[i-1]+1):(gaps[i]-1), 4] <- NA;
			}
		}
		out[(gaps[length(gaps)]+1):nrow(out), 4] <- NA
		out <- out[!is.na(out[,4]),];
	}
	return(out)
}




args <- commandArgs(trailingOnly=TRUE);
OUT <- vector();
for (f in args) {
print(f);

line <- unlist(strsplit(f, "_"))[2]
a <- readRDS(f);

require("GenomicRanges")
tmp1 <- matrix(a[grepl("eOrgan",a[,1]) & a[,3] == "dup", 1:5], ncol=5)
colnames(tmp1) <- colnames(a)[1:5]
tmp2 <- matrix(a[grepl("lOrgan",a[,1]) & a[,3] == "dup", 1:5], ncol=5)
colnames(tmp2) <- colnames(a)[1:5]
if (nrow(tmp1) > 0 & nrow(tmp2) > 0) {
	dups_e <- makeGRangesFromDataFrame(as.data.frame(tmp1, stringsAsFactors=FALSE), 
		start.field="st_bp", end.field="ed_bp")
	dups_l <- makeGRangesFromDataFrame(as.data.frame(tmp2, stringsAsFactors=FALSE), 
		start.field="st_bp", end.field="ed_bp")

	dups <-intersect_merge(dups_e, dups_l)
	if (is.null(dim(dups))) {
		dups <- data.frame(matrix(dups, ncol=5))
		colnames(dups) <- colnames(a)[1:5]
	}
	dups$type <- rep("dup", max(1,nrow(dups), na.rm=T))
	dups$source <- rep(line, max(1,nrow(dups), na.rm=T))
} else {dups <- rep(NA, 7)}

tmp1 <- matrix(a[grepl("eOrgan",a[,1]) & a[,3] == "del", 1:5], ncol=5)
colnames(tmp1) <- colnames(a)[1:5]
tmp2 <- matrix(a[grepl("lOrgan",a[,1]) & a[,3] == "del", 1:5], ncol=5)
colnames(tmp2) <- colnames(a)[1:5]
if (nrow(tmp1) > 0 & nrow(tmp2) > 0) {
	dels_e <- makeGRangesFromDataFrame(as.data.frame(tmp1, stringsAsFactors=FALSE), 
		start.field="st_bp", end.field="ed_bp")
	dels_l <- makeGRangesFromDataFrame(as.data.frame(tmp2, stringsAsFactors=FALSE), 
		start.field="st_bp", end.field="ed_bp")

	dels <-intersect_merge(dels_e, dels_l);
	if (is.null(dim(dels))) {
		dels <- data.frame(matrix(dels, ncol=5))
		colnames(dels) <- colnames(a)[1:5]
	}

	dels$type <- rep("del", max(1,nrow(dels), na.rm=T))
	dels$source <- rep(line, max(1,nrow(dels), na.rm=T))
} else {dels <- rep(NA, 7)}

if (!is.na(dups[1])) {
	OUT <- rbind(OUT,dups);
} 
if (!is.na(dels[1])) {
	OUT <- rbind(OUT,dels);
}

}
saveRDS(OUT, "Cleaned_Codex_CNVs.rds");
