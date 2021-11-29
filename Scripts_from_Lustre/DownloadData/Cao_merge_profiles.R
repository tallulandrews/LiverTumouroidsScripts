source("~/My_R_packages/TreeOfCells/R/Merge_Profiles.R")
profiles <- Sys.glob("Cao_Mmus_Organogenesis_Profiles*.rds")
all <- readRDS(profiles[1])
for (i in 2:length(profiles)) {
	b <- readRDS(profiles[i])
	all <- merge_lists(all, b, name1="cao1", name2="cao2")
}
saveRDS(all, "Cao_Mmus_Organogenesis_All_Profiles.rds")
merged <- merge_profiles(all)
saveRDS(merged, "Cao_Mmus_Organogenesis_Merged_Profiles.rds")
