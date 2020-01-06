library(karyoploteR)
library(regioneR)

a <- readRDS("Cleaned_Codex_CNVs.rds")

donor <- unique(a$source)

for (d in donor) {
	gain <- toGRanges(a[a$source==d & a$type == "dup",1:3])
	loss <- toGRanges(a[a$source==d & a$type == "del",1:3])
	kp <- plotKaryotype(genome="hg19")
	kpPlotRegions(kp, gain, col="red")
	kpPlotRegions(kp, loss, col="blue")
}

	
