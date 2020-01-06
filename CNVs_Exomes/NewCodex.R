library(CODEX)
library(methods)
dirPath <- "/lustre/scratch117/cellgen/team218/TA/LiverOrganoids/ExomeCNVs/DNA_BAMS"
bamFile <- list.files(dirPath, pattern = '*.bam$')
bamdir <- file.path(dirPath, bamFile)
sampname <- as.matrix(c("HCC10early_organoid", "HCC10late_organoid", "HCC10_TISSUE"), ncol=1)
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
	qcObj <- qc(Y, sampname, chr, ref, mapp, gc, cov_thresh = c(20, 4000),
	         length_thresh = c(20, 2000), mapp_thresh = 0.9, gc_thresh = c(20, 80))
	Y_qc <- qcObj$Y_qc; sampname_qc <- qcObj$sampname_qc; gc_qc <- qcObj$gc_qc
	mapp_qc <- qcObj$mapp_qc; ref_qc <- qcObj$ref_qc; qcmat <- qcObj$qcmat
	# This bit actually fitting thing.
	Y_qc <- cbind(Y_qc, Y_qc[,3])
	normObj <- normalize2(Y_qc, gc_qc, K = 1, normal_index=c(3,4))
	Yhat <- normObj$Yhat; AIC <- normObj$AIC; BIC <- normObj$BIC
	RSS <- normObj$RSS; K <- normObj$K
	optK = K[which.max(BIC)]
	# Call CNVs
	finalcall <- segment(Y_qc, Yhat, optK = optK, K = K, sampname_qc,
	         ref_qc, chr, lmax = 2000, mode = "fraction")
	saveRDS(finalcall, file=paste(chr, "Codex_allCalls.rds", sep="_"))
	# Filter CNVs
	call.qc <- finalcall[as.numeric(finalcall[,6])>500,]
	call.qc.sample1.dup <- call.qc[call.qc[,1] == sampname[1] & call.qc[,3] == "dup",]
	call.qc.sample1.del <- call.qc[call.qc[,1] == sampname[1] & call.qc[,3] == "del",]
	call.qc.sample2.dup <- call.qc[call.qc[,1] == sampname[2] & call.qc[,3] == "dup",]
	call.qc.sample2.del <- call.qc[call.qc[,1] == sampname[2] & call.qc[,3] == "del",]
	# Intersect
	if (is.null(nrow(call.qc.sample1.dup))) {
		call.qc.sample1.dup <- rbind(call.qc.sample1.dup,call.qc.sample1.dup)
	}
	if (is.null(nrow(call.qc.sample1.del))) {
		call.qc.sample1.del <- rbind(call.qc.sample1.del,call.qc.sample1.del)
	}
	if (is.null(nrow(call.qc.sample2.dup))) {
		call.qc.sample2.dup <- rbind(call.qc.sample2.dup,call.qc.sample2.dup)
	}
	if (is.null(nrow(call.qc.sample2.del))) {
		call.qc.sample2.del <- rbind(call.qc.sample2.del,call.qc.sample2.del)
	}
	require("bedr")
	regions1.dup <- paste(call.qc.sample1.dup[,2], ":", call.qc.sample1.dup[,4], "-", call.qc.sample1.dup[,5], sep="")
	regions1.del <- paste(call.qc.sample1.del[,2], ":", call.qc.sample1.del[,4], "-", call.qc.sample1.del[,5], sep="")
	regions2.dup <- paste(call.qc.sample2.dup[,2], ":", call.qc.sample2.dup[,4], "-", call.qc.sample2.dup[,5], sep="")
	regions2.del <- paste(call.qc.sample2.del[,2], ":", call.qc.sample2.del[,4], "-", call.qc.sample2.del[,5], sep="")
	regions1.dup <- bedr.merge.region(bedr.sort.region(regions1.dup, method = "natural"));
	regions1.del <- bedr.merge.region(bedr.sort.region(regions1.del, method = "natural"));
	regions2.dup <- bedr.merge.region(bedr.sort.region(regions2.dup, method = "natural"));
	regions2.del <- bedr.merge.region(bedr.sort.region(regions2.del, method = "natural"));
	intesect.dup.regions <- bedr.join.multiple.region(x=list(regions1.dup, regions2.dup))
	intesect.del.regions <- bedr.join.multiple.region(x=list(regions1.del, regions2.del))
	out <- list(int.dup=intesect.dup.regions, int.del=intesect.del.regions)
	saveRDS(out, file=paste(chr, "Codex_intCalls.rds", sep="_"))
}
