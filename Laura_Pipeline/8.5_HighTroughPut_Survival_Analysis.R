downloadTCGA <- function() {
	require("TCGAbiolinks")

	# Download
	hep_meta_query <- GDCquery(project="TCGA-LIHC", data.category="Clinical")
	chol_meta_query <- GDCquery(project="TCGA-CHOL", data.category="Clinical")
	hep_expr_query <- GDCquery(project="TCGA-LIHC", data.category="Transcriptome Profiling", data.type="Gene Expression Quantification", experimental.strategy="RNA-Seq", workflow.type="HTSeq - FPKM", sample.type="Primary solid Tumor")
	chol_expr_query <- GDCquery(project="TCGA-CHOL", data.category="Transcriptome Profiling", data.type="Gene Expression Quantification", experimental.strategy="RNA-Seq", workflow.type="HTSeq - FPKM", sample.type="Primary solid Tumor")

#	GDCdownload(hep_expr_query, method="api", directory="/lustre/scratch117/cellgen/team218/TA/LiverOrganoids/ExternalData/TCGA/Hepato", files.per.chunk=6)
#	GDCdownload(hep_meta_query, method="api", directory="/lustre/scratch117/cellgen/team218/TA/LiverOrganoids/ExternalData/TCGA/Hepato", files.per.chunk=6)
#	GDCdownload(chol_expr_query, method="api", directory="/lustre/scratch117/cellgen/team218/TA/LiverOrganoids/ExternalData/TCGA/Chol", files.per.chunk=6)
#	GDCdownload(chol_meta_query, method="api", directory="/lustre/scratch117/cellgen/team218/TA/LiverOrganoids/ExternalData/TCGA/Chol", files.per.chunk=6)

	# Process
	GDCprepare(hep_expr_query, save=TRUE, save.filename="HepExpr.RData", directory="/lustre/scratch117/cellgen/team218/TA/LiverOrganoids/ExternalData/TCGA/Hepato", summarizedExperiment=TRUE)
	GDCprepare(hep_meta_query, save=TRUE, save.filename="HepMeta.RData", directory="/lustre/scratch117/cellgen/team218/TA/LiverOrganoids/ExternalData/TCGA/Hepato", summarizedExperiment=TRUE)
}



loadTCGA <- function(file=c("/warehouse/team218_wh01/tallulah/TCGA/All_TCGA_Transcriptome.rds",
			    "/warehouse/team218_wh01/tallulah/TCGA/TCGA_Chol_Transcriptome.rds",
			    "/warehouse/team218_wh01/tallulah/TCGA/TCGA_Transcriptome.rds ")) {

	tcga_obj <- readRDS(file[1])
	my_stage <- as.character(tcga_obj$ann[,3])
	my_stage[my_stage == "stage:stage i, "] = 1
	my_stage[my_stage == "stage:stage ii, "] = 2
	my_stage[my_stage == "stage:stage iii, "] = 3
	my_stage[my_stage == "stage:stage iiia, "] = 3
	my_stage[my_stage == "stage:stage iiib, "] = 3
	my_stage[my_stage == "stage:stage iiic, "] = 3
	my_stage[my_stage == "stage:stage iv, "] = 4
	my_stage[my_stage == "stage:stage iva, "] = 4
	my_stage[my_stage == "stage:stage ivb, "] = 4
	my_stage[my_stage == "stage:not reported, "] = NA

	my_death_time <- as.character(tcga_obj$ann[,5])
	my_death_time<-matrix(unlist(strsplit(my_death_time, ":")), ncol=2, byrow=T)
	my_death_time<-sub(", ","", my_death_time[,2])
	my_death_time[my_death_time =="null"] = NA
	my_death_time<-as.numeric(my_death_time)

	my_followup_time <- as.character(tcga_obj$ann[,7])
	my_followup_time<-matrix(unlist(strsplit(my_followup_time, ":")), ncol=2, byrow=T)
	my_followup_time<-sub(", ","", my_followup_time[,2])
	my_followup_time[my_followup_time =="null"] = NA
	my_followup_time<-as.numeric(my_followup_time)


	sf<- colSums(tcga_obj$data)
	norm <- t(t(tcga_obj$data)/sf*1000000)
	norm <- log(norm+1)/log(2)
	return(list(counts=tcga_obj$data, norm=norm, death=my_death_time, followup=my_followup_time, stage=my_stage))
}


calcScores <- function(tcga, genes) {
	source("~/R-Scripts/Ensembl_Stuff.R")
	if (sum(genes %in% rownames(tcga$counts)) == 0){
		new_gene_names <- General_Map(rownames(tcga$counts), in.org="Hsap", out.org="Hsap", in.name="ensg", out.name="symbol")
		genes <- genes[genes != ""];
		rows <- new_gene_names %in% genes
	} else {
		rows <- rownames(tcga$counts) %in% genes
	}
	
	# scaling?

	# scores
	if (sum(rows) > 1) {
		scores <- colMeans(tcga$norm[rows,])
	} else if (sum(rows) > 0) {
		scores <- tcga$norm[rows,]
	} else {
		print ("No genes mapped.")
		scores <- rep(1, ncol(tcga$norm))
	}
	return(scores)
}

score_vs_stage <- function(tcga, score) {
	boxes <- boxplot(score~tcga$stage, notch=TRUE, xlab="Stage", ylabl="Score", col=c("#4575b4","#fee090","#d73027","#a50026"))
	p = wilcox.test(score[tcga$stage=="1"], score[tcga$stage=="3"])
	yes <- max(boxes$stats[4,1:3])*1.05
	xes <- c(1-0.15,3+0.15)
	lines(xes, c(yes, yes), lwd=2)
	text(2, yes, paste("p =", signif(p$p.value, digits=2)), font=2, pos=3)
	return(p)
}

score_vs_death <- function(tcga, score, top_quantile=0.5) {
	hi <- score > quantile(score, probs=top_quantile)
	low <- score <= quantile(score, probs=1-top_quantile)
	boxes <- boxplot(tcga$death[hi], tcga$death[low], notch=TRUE, names=c("High", "Low"), xlab="Score", ylab="Time to Death (days)", col=c("#d73027", "#4575b4"))
	p <- wilcox.test(tcga$death[hi], tcga$death[low])
	yes <- max(boxes$stats[4,1:2])*1.1
	xes <- c(1-0.15,2+0.15)
	lines(xes, c(yes, yes), lwd=2)
	text(1.5, yes, paste("p =", signif(p$p.value, digits=2)), font=2, pos=3)
}
	
do_survival_test <- function(tcga, score, top_quantile=0.5) {
	require("survival")
	dead = !is.na(tcga$death);
	times <- tcga$followup;
	times[dead] <- tcga$death[dead]

	# top X vs bottom X
	predictor <- as.numeric(score > quantile(score, probs=top_quantile))
	predictor[score < quantile(score, probs=top_quantile) & score > quantile(score, probs=1-top_quantile)] <- NA

	keep <- !is.na(times) & !is.na(predictor) & times > 0;
	predictor<-predictor[keep];
	times<-times[keep]
	dead <- dead[keep]

	surv_tcga_obj <- Surv(times, dead)
	surv_test <- survdiff(surv_tcga_obj~predictor) # Discrete predictor
	curves <- survfit(surv_tcga_obj ~ predictor)
	plot(curves, col=c("blue", "red"), mark.time=TRUE, conf.int=TRUE, xlab="Days", ylab="Survival")
	legend("bottomleft", c("High", "Low"), col=c("red","blue"), title="Score", lty=1, bty="n")
	legend("topright", bty="n", paste("p =", signif(1-pchisq(surv_test$chisq, df=1), digits=2)))
} 

# For Patricia
tcga <- loadTCGA()
gene_files <- paste("/nfs/users/nfs_t/ta6/Collaborations/LiverOrganoids/", c("Patricia_CC1_CSC_Markers.txt", "Patricia_CC1_Chol_Markers.txt", "Patricia_HCC10_CSC_Markers.txt", "Patricia_HCC10_Hep_Markers.txt"), sep="")
for (f in gene_files) {
	genes <- read.table(f, stringsAsFactors=FALSE);
	for (g in genes[,1]) {
		score <- calcScores(tcga, g)
		png(paste(g,"_Survival.png", sep=""), width=5*2.75, height=5, units="in", res=300)
		par(mfrow=c(1,3))
		do_survival_test(tcga, score, 0.5)
		score_vs_death(tcga, score)
		score_vs_stage(tcga, score)
		dev.off();
	}
}

# For Montreal
tcga <- loadTCGA()
genes <- read.table("/lustre/scratch117/cellgen/team218/TA/LiverOrganoids/Analysis_Pipeline_Output/Consistent_Prolif_notCC_Alt.txt")
score <- calcScores(tcga, rownames(genes))
png("Consistent_Prolif_notCC_Survival_Alt.png", width=5*1.75, height=5, units="in", res=300)
par(mfrow=c(1,2))
do_survival_test(tcga, score, 0.5)
score_vs_death(tcga, score)
dev.off()

png("Consistent_Prolif_notCC_Stage_Alt.png", width=5, height=5, units="in", res=300)
score_vs_stage(tcga, score)
dev.off()

# For Patricia

