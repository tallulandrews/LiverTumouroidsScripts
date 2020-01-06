library(CODEX)
library(methods)
dirPath <- "/lustre/scratch117/cellgen/team218/TA/LiverOrganoids/ExomeCNVs/DNA_BAMS"
bamFile <- list.files(dirPath, pattern = '*.bam$')
bamdir <- file.path(dirPath, bamFile)
samples <- read.table("/lustre/scratch117/cellgen/team218/TA/LiverOrganoids/ExomeCNVs/DNA_BAMS/sample_names.txt", header=T)
sampname <- as.matrix(paste(samples[,2], samples[,3], sep="-"), ncol=1);
ctrl <- which(samples[,2] %in% c("Donor3", "Donor4"))

bedFile <- "/lustre/scratch117/cellgen/team218/TA/genomebuilding/Homo_sapiens.GRCh38.79.fixed.bed"
for (chr in paste("chr",18:22, sep="")) {
	#chr <- "chr22"
	# Set Up data
	bambedObj <- getbambed(bamdir = bamdir, bedFile = bedFile,
        	 sampname = sampname, projectname = "CODEX_demo", chr)
	bamdir <- bambedObj$bamdir; sampname <- bambedObj$sampname
	ref <- bambedObj$ref; projectname <- bambedObj$projectname; chr <- bambedObj$chr
	coverageObj <- getcoverage(bambedObj, mapqthres = 20)
	Y <- coverageObj$Y; readlength <- coverageObj$readlength
	gc <- getgc(chr, ref)
	gc[is.na(gc)] <- gc[min(which(is.na(gc)))-1]
	mapp <- getmapp(chr, ref)
	mapp[is.na(mapp)] <- mapp[min(which(is.na(mapp)))-1]
	# QC
	qcObj <- qc(Y[,ctrl], sampname, chr, ref, mapp, gc, cov_thresh = c(20, 10000),
	         length_thresh = c(20, 10000), mapp_thresh = 0.9, gc_thresh = c(20, 80))
	kept <- as.logical(qcObj$qcmat[,4]);

	Y_qc <- Y[kept,]; sampname_qc <- sampname; gc_qc <- gc[kept]
	mapp_qc <- mapp[kept]; ref_qc <- ref[kept,]; qcmat <- qcObj$qcmat
	# This bit actually fitting thing.
	normObj <- normalize2(Y_qc, gc_qc, K = 1:(length(ctrl)-1), normal_index=ctrl)
	Yhat <- normObj$Yhat; AIC <- normObj$AIC; BIC <- normObj$BIC
	RSS <- normObj$RSS; K <- normObj$K
	optK = K[which.max(BIC)]
	# Call CNVs
	finalcall <- segment(Y_qc, Yhat, optK = optK, K = K, sampname_qc,
	         ref_qc, chr, lmax = 2000, mode = "fraction")
	saveRDS(finalcall, file=paste(chr, "Codex_allCalls.rds", sep="_"))
	# The below must be updated for multiple samples.
	# Filter CNVs
	call.qc <- finalcall[as.numeric(finalcall[,6])>500,]
	list_dups <- list();
	list_dels <- list();
	for (line in samples[,2]) {
		if (line %in% c("Donor3", "Donor4")) {next;}
		for (sample in sampname[grep(line, sampname[,1]),]) {
		call.qc.sample.dup <- call.qc[call.qc[,1] == line & call.qc[,3] == "dup",]
		call.qc.sample.del <- call.qc[call.qc[,1] == line & call.qc[,3] == "del",]
		if (is.null(nrow(call.qc.sample1.dup))) {
			call.qc.sample.dup <- rbind(call.qc.sample.dup,call.qc.sample.dup)
		}
		if (is.null(nrow(call.qc.sample1.del))) {
			call.qc.sample.del <- rbind(call.qc.sample.del,call.qc.sample.del)
		}
		# Intersect
		require("bedr")
		regions.dup <- paste(call.qc.sample.dup[,2], ":", 
				     call.qc.sample.dup[,4], "-", 
				     call.qc.sample.dup[,5], sep="")
		regions.del <- paste(call.qc.sample.del[,2], ":", 
				     call.qc.sample.del[,4], "-", 
				     call.qc.sample1.del[,5], sep="")
		regions.dup <- bedr.merge.region(bedr.sort.region(regions.dup, method = "natural"));
		regions.del <- bedr.merge.region(bedr.sort.region(regions.del, method = "natural"));
	}
	intesect.dup.regions <- bedr.join.multiple.region(x=list(regions1.dup, regions2.dup))
	intesect.del.regions <- bedr.join.multiple.region(x=list(regions1.del, regions2.del))
	out <- list(int.dup=intesect.dup.regions, int.del=intesect.del.regions)
	saveRDS(out, file=paste(chr, "Codex_intCalls.rds", sep="_"))
}
