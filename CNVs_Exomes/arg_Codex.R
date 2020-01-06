args <- commandArgs(trailingOnly=TRUE)

library(CODEX)
library(methods)
dirPath <- "/lustre/scratch117/cellgen/team218/TA/LiverOrganoids/ExomeCNVs/DNA_BAMS"
bamFile <- list.files(dirPath, pattern = '*.bam$')
bamdir <- file.path(dirPath, bamFile)
samples <- read.table("/lustre/scratch117/cellgen/team218/TA/LiverOrganoids/ExomeCNVs/DNA_BAMS/sample_names.txt", header=T)
sampname <- as.matrix(paste(samples[,2], samples[,3], sep="-"), ncol=1);
ctrl <- which(samples[,2] %in% c("Donor3", "Donor4"))

bedFile <- "/lustre/scratch117/cellgen/team218/TA/genomebuilding/Homo_sapiens.GRCh38.79.fixed.bed"
for (chr in paste("chr",args, sep="")) {
	out <- vector();
	out_intersected_dups <- vector()
	out_intersected_dels <- vector()
for (s in unique(samples$Donor)) {
	out_file=paste(chr,s, "Codex_vsCtrl_allCalls.rds", sep="_")
	if(file.exists(out_file)) {next;}
	samp_ids <- c(ctrl, which(samples$Donor == s))
	bamdir_s <- bamdir[samp_ids];
	samp_names <- as.matrix(sampname[samp_ids,], ncol=1)
	ctrl_s <- 1:length(ctrl);
	

	#chr <- "chr22"
	# Set Up data
	bambedObj <- getbambed(bamdir = bamdir_s, bedFile = bedFile,
        	 sampname = samp_names, projectname = paste("CODEX",s, "vs_ctrl", sep="_"), chr)
	bamdir_s <- bambedObj$bamdir; samp_names <- bambedObj$sampname
	ref <- bambedObj$ref; chr <- bambedObj$chr

	# get data for confounders
	coverageObj <- getcoverage(bambedObj, mapqthres = 20)
	Y <- coverageObj$Y; readlength <- coverageObj$readlength
	gc <- getgc(chr, ref)
	gc[is.na(gc)] <- gc[min(which(is.na(gc)))-1]
	mapp <- getmapp(chr, ref)
	mapp[is.na(mapp)] <- mapp[min(which(is.na(mapp)))-1]
	# QC
	qcObj <- qc(Y[,ctrl_s], samp_names, chr, ref, mapp, gc, cov_thresh = c(20, 10000),
	         length_thresh = c(20, 10000), mapp_thresh = 0.9, gc_thresh = c(20, 80))
	kept <- as.logical(qcObj$qcmat[,4]);

	Y_qc <- Y[kept,]; sampname_qc <- samp_names; gc_qc <- gc[kept]
	mapp_qc <- mapp[kept]; ref_qc <- ref[kept,]; qcmat <- qcObj$qcmat
	# This bit actually fitting thing.
	normObj <- normalize2(Y_qc, gc_qc, K = 1:(length(ctrl_s)-1), normal_index=ctrl_s)
	Yhat <- normObj$Yhat; AIC <- normObj$AIC; BIC <- normObj$BIC
	RSS <- normObj$RSS; K <- normObj$K
	optK = K[which.max(BIC)]
	# Call CNVs
	finalcall <- segment(Y_qc, Yhat, optK = optK, K = optK, sampname_qc,
	         ref_qc, chr, lmax = 5000, mode = "fraction")
	call.qc <- finalcall[as.numeric(finalcall[,6]) > 500 & abs(as.numeric(finalcall[,11])-2) > 0.1,]
	saveRDS(call.qc, file=out_file)
	# Merged
#	dups <- call.qc[call.qc[,3] == "dup",]
#	dels <- call.qc[call.qc[,3] == "del",]
#	if (is.null(nrow(dups))) {
#		dups <- rbind(dups,dups)
#	}
#	if (is.null(nrow(dels))) {
#		dels <- rbind(dels,dels)
#	}
#	require("bedr")
#	regions.dup <- paste(dups[,2], ":", 
#			     dups[,4], "-", 
#			     dups[,5], sep="")
#	regions.del <- paste(dels[,2], ":", 
#			     dels[,4], "-", 
#			     dels[,5], sep="")
#	regions.dup <- bedr.merge.region(bedr.sort.region(regions.dup, method = "natural"));
#	regions.del <- bedr.merge.region(bedr.sort.region(regions.del, method = "natural"));
#	out_intersected_dups <- rbind(out_intersected_dups, regions.dup)
#	out_intersected_dels <- rbind(out_intersected_dels, regions.del)

}
#	saveRDS(list(dup=out_intersected_dups, del=out_intersected_dels), file=paste(chr, "Codex_vsCtrl_mergeCalls.rds", sep="_"))

}
