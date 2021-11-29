# Human embryos and ESCs
# Yan et al. (2013) Single-cell RNA-Seq profiling of human preimplantation embryos and embryonic stem cells. Nature Structural Molecular Biology. 20(9):1131-9
# GEO: GSE36552
# PMID: 23934149

wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE36nnn/GSE36552/matrix/GSE36552_series_matrix.txt.gz
gunzip GSE36552_series_matrix.txt.gz
perl /nfs/users/nfs_t/ta6/Data_Processing_Scripts/parse_series_matrix.pl GSE36552_series_matrix.txt 4 > Anno.txt

