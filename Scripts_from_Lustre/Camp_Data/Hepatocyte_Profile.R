# array data for hiPSC-derived cholangiocyte-like cells 
# Sampaziotis et al. (2015). Cholangiocytes derived from human induced pluripotent stem cells for disease modeling and drug validation. Nature Biotechnology. 33:845-852
# PMID: 26167629
# ArrayExpress: E-MTAB-2965
#require("lumi") # I think this is already QCed & quantile normalized

require("scater")
Camp_QC = readRDS("Camp_QC.rds")
Camp_QC<-Camp_QC[,pData(Camp_QC)$cell_type1 != "Kupffer"]
# Shorten names so combination names aren't ridiculously long
my_shortened_celltype <- as.character(pData(Camp_QC)$cell_type1)
my_shortened_celltype[my_shortened_celltype == "definitive endoderm"] <- "de"
my_shortened_celltype[my_shortened_celltype == "hepatic endoderm"] <- "he"
my_shortened_celltype[my_shortened_celltype == "mesenchymal stem cell"] <- "msc"
my_shortened_celltype[my_shortened_celltype == "immature hepatoblast"] <- "ih"
my_shortened_celltype[my_shortened_celltype == "mature hepatocyte"] <- "Mhep"
my_shortened_celltype[my_shortened_celltype == "erythroblasts"] <- "B_erythro"
my_shortened_celltype[my_shortened_celltype == "lymphoblasts"] <- "B_lympho"
my_shortened_celltype[my_shortened_celltype == "fetal hepatocytes"] <- "Fhep"
my_shortened_celltype[my_shortened_celltype == "stellate"] <- "stellate"
pData(Camp_QC)$cell_type2 <- factor(my_shortened_celltype);
my_shortened_source <- as.character(pData(Camp_QC)$Source)
my_shortened_source[my_shortened_source=="iPSC line TkDA3-4"] <- "iPSC line"
my_shortened_source[my_shortened_source=="Mesenchymal stem cell"] <- "MSC"
my_shortened_source[my_shortened_source=="liver bud"] <- "Liver bud"
pData(Camp_QC)$Source2 <- factor(my_shortened_source);



require("limma")
require("edgeR")
labels = factor(my_shortened_celltype)
design = model.matrix(~labels, data=labels)
limma_obj <- DGEList(counts=counts(Camp_QC))
limma_obj <- calcNormFactors(limma_obj, method="TMM")
voomed <- voom(limma_obj, design, plot=FALSE)
fit <- lmFit(voomed, design)
contrast.matrix <- makeContrasts(labelsFhep-labelsMhep, labelsipsc-labelshepatic, labelsipsc-labelsMhep, labelshe-labelshepatic, labelshe-labelsMhep, labelsipsc+labelshe-labelsMhep-labelshepatic, labelsde-labelsMhep, levels=design)
fit_contrasts <- contrasts.fit(fit, contrast.matrix)
fit_contrasts <- eBayes(fit_contrasts)
results <- decideTests(fit_contrasts)
require("CellTypeProfiles")
averages <- my_row_mean_aggregate(exprs(Camp_QC), labels)
#fit <- eBayes(fit)

SC_min <- apply(averages[,c("de","ipsc","he")], 1, min)
SC_max <- apply(averages[,c("de","ipsc","he")], 1, max)
Hep_max <- apply(averages[,c("Mhep","hepatic")], 1, max)
Hep_min <- apply(averages[,c("Mhep","hepatic")], 1, min)
SC_min[SC_min == 0] = 1/1000000
SC_max[SC_max == 0] = 1/1000000
Hep_max[Hep_max == 0] = 1/1000000
Hep_min[Hep_min == 0] = 1/1000000

Hep_score <- Hep_min/SC_max
SC_score <- SC_min/Hep_max

Hep_up <- Hep_min/SC_max > 2 & results[,6] < 0
SC_up <- SC_min/Hep_max > 2 & results[,6] > 0

OUT <- list(HepUp=Hep_score[Hep_up], SCUp=SC_score[SC_up])

saveRDS(OUT, "Hepatocyte_Profile.rds")

