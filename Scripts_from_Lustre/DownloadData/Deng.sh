# Deng et al. Single-cell RNA-seq reveals dynamic, ranomd monoallelic gene expression in mammalian cells. Science (2014). 343(6167):193-6
# PMID: 24408435
# GEO: GSE45719

wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE45nnn/GSE45719/matrix/GSE45719_series_matrix.txt.gz
gunzip GSE45719_series_matrix.txt.gz
perl ~/Data_Processing_Scripts/parse_series_matrix.pl GSE45719_series_matrix.txt 2 > Deng_ann.txt
mv ExprMat.txt Deng_counts.txt

