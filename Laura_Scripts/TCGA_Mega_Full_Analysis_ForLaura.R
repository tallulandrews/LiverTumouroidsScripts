obj <- readRDS("/warehouse/team218_wh01/tallulah/TCGA/All_TCGA_Transcriptome.rds")
my_stage <- as.character(obj$ann[,3])
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

my_death_time <- as.character(obj$ann[,5])
my_death_time<-matrix(unlist(strsplit(my_death_time, ":")), ncol=2, byrow=T)
my_death_time<-sub(", ","", my_death_time[,2])
my_death_time[my_death_time =="null"] = NA
my_death_time<-as.numeric(my_death_time)

my_followup_time <- as.character(obj$ann[,7])
my_followup_time<-matrix(unlist(strsplit(my_followup_time, ":")), ncol=2, byrow=T)
my_followup_time<-sub(", ","", my_followup_time[,2])
my_followup_time[my_followup_time =="null"] = NA
my_followup_time<-as.numeric(my_followup_time)


norm<- colSums(obj$data)
normed <- t(t(obj$data)/norm*1000000)
normed <- log(normed+1)/log(2)

profiles <- read.delim("/lustre/scratch117/cellgen/team218/TA/LiverOrganoids/WeightedMean_Profiles.csv", sep=",", header=T)

profile_score <- (profiles$CCA1_4+profiles$CCA1_1+profiles$HCC6_4+profiles$HCC10_1)-rowSums(profiles[,!(colnames(profiles) %in% c("CCA1_4","CCA1_1","HCC6_4","HCC10_1"))])

# Laura-esque
normal <- apply(profiles[,!(colnames(profiles) %in% c("CCA1_4","CCA1_1","HCC6_4","HCC10_1"))], 1, mean)
stemy <- apply(profiles[,(colnames(profiles) %in% c("CCA1_4","CCA1_1","HCC6_4","HCC10_1"))], 1, min)
#profile_score2 <- (stemy-normal)/abs(normal) # no work b/c of vals close to 0.
profile_score2<- stemy-normal

profile_score = profile_score2

profiles <- profiles[order(-1*profile_score),]

# Consistent gene set #
normed <- normed[match(rownames(profiles),rownames(normed)),]
profiles<-profiles[rownames(profiles) %in% rownames(normed),]
normed<-normed[rownames(normed) %in% rownames(profiles),]

# General score #
stem_score_weights <- (profiles$CCA1_4+profiles$CCA1_1+profiles$HCC6_4+profiles$HCC10_1)/4
pat_stem_score <- t(normed) %*% stem_score_weights
png("Stage_vs_Score.png", width=6, height=6, units="in", res=300)
boxplot(pat_stem_score~my_stage, notch=TRUE, xlab="Stage", ylab="Stem-ness Score", col=c("#4575b4","#fee090","#d73027","#a50026"))
p = wilcox.test(pat_stem_score[my_stage=="1"], pat_stem_score[my_stage=="3"])
lines(c(1-0.15,3+0.15), c(11000,11000), lwd=2)
text(2, 11000, paste("p =", round(p$p.value, digits=2)), font=2, pos=3)
dev.off()

# Exteme score #
pat_extreme_score <- colSums(normed[1:100,]) - colSums(normed[length(normed[,1])-1:100,]) # tried 20, 50, 100


hi <- pat_stem_score > median(pat_stem_score)
png("Score_vs_Death.png", width=6, height=6, units="in", res=300)
boxplot(my_death_time~hi, notch=TRUE, names=c("Low", "High"), xlab="Stem-ness Score", ylab="Time to Death (days)", col=c("#4575b4","#d73027"))
p <- wilcox.test(my_death_time[hi], my_death_time[!hi])$p.value
lines(c(1-0.15,2+0.15), c(1600,1600), lwd=2)
text(1.5, 1600, paste("p =",round(p, digits=4)), font=2, pos=3)
dev.off()

hi <- pat_extreme_score > median(pat_extreme_score)
png("Extreme_Score_vs_Death.png", width=6, height=6, units="in", res=300)
boxplot(my_death_time~hi, notch=TRUE, names=c("Low", "High"), xlab="Stem-ness Score", ylab="Time to Death (days)", col=c("#4575b4","#d73027"))
p <- wilcox.test(my_death_time[hi], my_death_time[!hi])$p.value
lines(c(1-0.15,2+0.15), c(1600,1600), lwd=2)
text(1.5, 1600, paste("p =",round(p, digits=4)), font=2, pos=3)
dev.off()


# Survival Test
require("survival")
dead = !is.na(my_death_time);
times <- my_followup_time;
times[dead] <- my_death_time[dead]

predictor <- as.numeric(pat_stem_score > median(pat_stem_score))
#predictor <- as.numeric(pat_extreme_score > median(pat_extreme_score))
predictor <- factor(predictor)

keep <- (dead | !is.na(my_followup_time)) & !is.na(predictor);
predictor<-predictor[keep];
times<-times[keep]
dead <- dead[keep]

surv_obj <- Surv(times, dead)
surv_test <- survdiff(surv_obj~predictor) # Discrete predictor
curves <- survfit(surv_obj ~ predictor)
png("Stem-ness_Survival_Curve.png", width=6, height=6, units="in", res=300)
plot(curves, col=c("blue", "red"), mark.time=TRUE, conf.int=TRUE, xlab="Days", ylab="Survival")
legend("bottomleft", c("High", "Low"), col=c("red","blue"), title="Stem-ness Score", lty=1, bty="n")
legend("topright", bty="n", paste("p =", signif(1-pchisq(surv_test$chisq, df=1), digits=1)))
dev.off()


predictor<-pat_stem_score
predictor<-predictor[keep];
surv_test1 <- coxph(surv_obj~predictor) # Continuous predictor?
surv_test2 <- survreg(surv_obj~predictor, dist="logistic") # Continuous predictor

sng_gene_test <- function(x) {
	require("survival")
	dead = !is.na(my_death_time);
	times <- my_followup_time;
	times[dead] <- my_death_time[dead]
	predictor<-x
	keep <- (dead | !is.na(my_followup_time)) & !is.na(predictor);
	predictor<-predictor[keep];
	times<-times[keep]
	dead <- dead[keep]

	surv_obj <- Surv(times, dead)
	surv_test1 <- coxph(surv_obj~predictor) # Cox proportional hazards model 
#	surv_test2 <- survreg(surv_obj~predictor, dist="logistic") # Continuous predictor
	p = pchisq(surv_test1$wald.test, df=1, lower.tail=FALSE)
	eff <- surv_test1$coefficients[1]
	return(c(eff,p))
}

out <- t(apply(normed, 1, sng_gene_test))
q.vals <- p.adjust(out[,2], method="fdr")
pat_association <- data.frame(eff_size = out[,1], p.value=out[,2], q.value=q.vals)

require("CellTypeProfiles")
require("scater")
get_log_norm <- function(SCE) {
        dat <- counts(SCE);
        fac <- colSums(dat);
        norm <- t(t(dat)/fac*median(fac))
        return(log(norm +1)/log(2));
}

CCA1 <- readRDS("/lustre/scratch117/cellgen/team218/TA/LiverOrganoids/CCA1_SC3.rds")
HCC6 <- readRDS("/lustre/scratch117/cellgen/team218/TA/LiverOrganoids/HCC6_SC3.rds")
HCC10 <- readRDS("/lustre/scratch117/cellgen/team218/TA/LiverOrganoids/HCC10_SC3.rds")
exprs(CCA1) <- get_log_norm(CCA1)
exprs(HCC6) <- get_log_norm(HCC6)
exprs(HCC10) <- get_log_norm(HCC10)
all_markers <- readRDS("/lustre/scratch117/cellgen/team218/TA/LiverOrganoids/Laura_Complex_Markers.rds")

CCA1_Means_mat <- my_row_mean_aggregate(exprs(CCA1), pData(CCA1)$sc3_6_clusters)
HCC6_Means_mat <- my_row_mean_aggregate(exprs(HCC6), pData(HCC6)$sc3_6_clusters)
HCC10_Means_mat <- my_row_mean_aggregate(exprs(HCC10), pData(HCC10)$sc3_6_clusters)
colnames(CCA1_Means_mat) <- paste("CCA1",colnames(CCA1_Means_mat), sep="_")
colnames(HCC6_Means_mat) <- paste("HCC6",colnames(HCC6_Means_mat), sep="_")
colnames(HCC10_Means_mat) <- paste("HCC10",colnames(HCC10_Means_mat), sep="_")

CCA1_strict_markers <- complex_markers(exprs(CCA1), pData(CCA1)$sc3_6_clusters, n_max=5, strict_only=TRUE)
HCC6_strict_markers <- complex_markers(exprs(HCC6), pData(HCC6)$sc3_6_clusters, n_max=5, strict_only=TRUE)
HCC10_strict_markers <- complex_markers(exprs(HCC10), pData(HCC10)$sc3_6_clusters, n_max=5, strict_only=TRUE)

source("~/R-Scripts/Ensembl_Stuff.R")
Symbol <- map_symbol_ensg(rownames(exprs(CCA1)), is.org="Hsap", is.name="ensg")

# Full Table
FULL <- data.frame(Symbol=Symbol, cbind(CCA1_Means_mat,HCC6_Means_mat,HCC10_Means_mat), CCA1_AUC=all_markers$CCA1[,1], HCC6_AUC=all_markers$HCC6[,1], HCC10_AUC=all_markers$HCC10[,1])

annotate_markers <- function(i) {
	FULL[i,2:22] <- signif(FULL[i,2:22], digits=2)
#	if (sum(CCA1_strict_markers[i,2:7]) > 0) {
#		col_indices <- which(CCA1_strict_markers[i,2:7]==1)+1
#		FULL[i,col_indices] <- paste("*", FULL[i,col_indices], sep="")
#	}
#	if (sum(HCC6_strict_markers[i,2:7]) > 0) {
#		col_indices <- which(HCC6_strict_markers[i,2:7]==1)+1+6
#		FULL[i,col_indices] <- paste("*", FULL[i,col_indices], sep="")
#	}
#	if (sum(HCC10_strict_markers[i,2:7]) > 0) {
#		col_indices <- which(HCC10_strict_markers[i,2:7]==1)+1+6+6
#		FULL[i,col_indices] <- paste("*", FULL[i,col_indices], sep="")
#	}
	return(FULL[i,])
}
annotated <- sapply(1:length(FULL[,1]), annotate_markers)
annotated<- as.data.frame(t(annotated))
rownames(annotated) = rownames(FULL)
annotated$Symbol <- FULL$Symbol
#here
CCA1_combo_labs <- get_combo_names(CCA1_strict_markers)
HCC6_combo_labs <- get_combo_names(HCC6_strict_markers)
HCC10_combo_labs <- get_combo_names(HCC10_strict_markers)

reorder <- c(1,7,12,16,10,19,18,17,4,9,8,2,11,5,14,20,21,22,6,13,15,3)
annotated <- annotated[,reorder]
FULL <- FULL[,reorder]

annotated$CCA1_Good_Mark <- CCA1_combo_labs
annotated$HCC6_Good_Mark <- HCC6_combo_labs
annotated$HCC10_Good_Mark <- HCC10_combo_labs
write.table(as.matrix(annotated), file="ForLaura_Mega_FullTable.csv", sep=",")
saveRDS(annotated, file="Mega_Analysis.rds")

death_association <- pat_association[match(rownames(FULL), rownames(pat_association)),]

# Cholandrocyte Table
denom <- apply(FULL[,5:15], 1, max)
denom[denom == 0] = 1/1000;
FC <- apply(FULL[,2:3], 1, min)-denom
C_tab <- cbind(annotated, death_association)
C_tab <- C_tab[order(-FC),]
C_tab<-C_tab[1:sum(FC>=0),]
write.table(as.matrix(C_tab), file="Mega_Chol_Table.csv", sep=",")

# Hepatocyte Table
denom <- apply(FULL[,c(2:4,12:15)], 1, max)
denom[denom == 0] = 1/1000;
FC <- apply(FULL[,5:11], 1, min)-denom
H_tab <- cbind(annotated, death_association)
H_tab <- H_tab[order(-FC),]
H_tab<-H_tab[1:sum(FC>=0),]
write.table(as.matrix(H_tab), file="Mega_Hepato_Table.csv", sep=",")

# Stem-cell Table
denom <- apply(FULL[,2:11], 1, max)
denom[denom == 0] = 1/1000;
FC <- apply(FULL[,12:15], 1, min)-denom
S_tab <- cbind(annotated, death_association)
S_tab <- S_tab[order(-FC),]
S_tab<-S_tab[1:sum(FC>=0),]
write.table(as.matrix(S_tab), file="Mega_Stem_Table.csv", sep=",")





#source("~/R-Scripts/Ensembl_Stuff.R")
#ForLaura <- profiles[,c(6,17,9,5,18,8,2,15,12,11,9,3,14,13,1,16,4,7)]
#ForLaura <- data.frame(ForLaura)
#ForLaura$profile_score <- profile_score
#ForLaura$effect_size =out[,1] 
#ForLaura$p.value=out[,2]
#ForLaura$q.value=q.vals
#ForLaura$Symbol <- map_symbol_ensg(rownames(ForLaura), is.org="Hsap", is.name="ensg")
#
#write.table(ForLaura, file="WeightedMean_Profile_PlusTCGA.csv", sep=",")




