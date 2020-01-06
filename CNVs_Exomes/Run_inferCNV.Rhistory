require(inferCNV)
require(InferCNV)
require(Infercnv)
require(infercnv)
infercnv_obj <- CreateInfercnvObject(raw_counts_matrix="inferCNV_input_counts.tsv", annotations_file="inferCNV_input_ann.txt", delim="\t", gene_order_file="inferCNV_input_geneorder.txt", ref_group_names=c("D3DM", "D3EM", "D9DM", "D9EM"))
infercnv_obj = infercnv::run(infercnv_obj, cutoff=1, out_dir="inferCNVRes/", cluster_by_groups=T, include.spike=F)
infercnv_obj = infercnv::run(infercnv_obj, cutoff=1, out_dir="inferCNVRes/", cluster_by_groups=T, include.spike=T)
infercnv_obj
savehistory("Run_inferCNV.Rhistory")
