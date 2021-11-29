# Human placenta
# 
# GEO: GSE87726
# PMID: 28174237


x1 = read.table("GSE87720_Placenta_singlecell_Batch1_tpm.txt.gz", header=TRUE)
x2 = read.table("GSE87720_Placenta_singlecell_Batch2_tpm.txt.gz", header=TRUE)
y = read.delim("GSE87708_humanplacenta_cufflinks.txt.gz", header=TRUE)
z = read.table("GSE87725_Syncytio.txt.gz", header=TRUE)

fix_Excel_names <- function(genes) {
	genes <- as.character(genes);
	#Mar, Sep, Dec
	fix <- grep("-Mar", genes)
	if (length(fix) > 0) {
		fix_with <- matrix(unlist(strsplit(genes[fix],"-")), ncol=2, byrow=T)
		genes[fix] <- paste("MARCH", fix_with[,1], sep="");
	}

	fix <- grep("-Sep", genes)
	if (length(fix) > 0) {
		fix_with <- matrix(unlist(strsplit(genes[fix],"-")), ncol=2, byrow=T)
		genes[fix] <- paste("SEPT", fix_with[,1], sep="");
	}

	fix <- grep("-Dec", genes)
	if (length(fix) > 0) {
		fix_with <- matrix(unlist(strsplit(genes[fix],"-")), ncol=2, byrow=T)
		genes[fix] <- paste("DEC", fix_with[,1], sep="");
	}
	return(genes);
}

fix_mat <- function(x) {
	fixed_names1 <- fix_Excel_names(x[,1])
	x <- x[,-1]
	remove <- duplicated(fixed_names1);
	x <- x[!remove,]
	fixed_names1 <- fixed_names1[!remove];
	rownames(x) <- fixed_names1
	return(x)
}

x1 <- fix_mat(x1)
x2 <- fix_mat(x2)
y <- fix_mat(y[,-c(1,4,6)])
z <- fix_mat(z)

rownames(y) <- sub(" ","", rownames(y)) # yes do this twice
rownames(y) <- sub(" ","", rownames(y))

# Bulk
rownames(y) <- sub(" ","", rownames(y))
saveRDS(list(placenta=y, syncytio=z), file="pavlicev_bulk.rds")


# Single-cell
batch <- rep(c(1,2), times=c(ncol(x1), ncol(x2)))
Full <- cbind(x1, x2)
Full <- Full[-c(30775:nrow(Full)),] # remove weird stuff

A <- data.frame(Source=rep("placenta", times=ncol(Full)), Species=rep("Human", times=ncol(Full)), Batch=batch);
rownames(A) <- colnames(Full)


require("scater")
pd <- new("AnnotatedDataFrame", data=as.data.frame(A))
pavlicev <- newSCESet(fpkmData=as.matrix(Full), phenoData=pd, logExprsOffset=1)

saveRDS(pavlicev, "pavlicev.rds")

