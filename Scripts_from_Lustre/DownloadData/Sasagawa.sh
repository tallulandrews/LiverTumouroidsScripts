# Mouse ESCs cell-cycle staged
# Sasagawa et al. (2013) Quartz-Seq: a highly reproducible and sensitive single-cell RNA sequencing method, reveals non-genetic gene-expression heterogeneity. Genome Biology. 14:3097
# doi: 10.1186/gb-2013-14-4-r31

wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE42nnn/GSE42268/matrix/GSE42268-GPL13112_series_matrix.txt.gz -O Sasagawa_matrix1.txt.gz
gunzip Sasagawa_matrix1.txt.gz
perl ~/Data_Processing_Scripts/parse_series_matrix.pl Sasagawa_matrix1.txt > Sasagawa_Ann1.txt 2
mv ExprMat.txt Sasagawa_ExprMat1.txt
rm Sasagawa_matrix1.txt
wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE42nnn/GSE42268/matrix/GSE42268-GPL15103_series_matrix.txt.gz -O Sasagawa_matrix2.txt.gz
gunzip Sasagawa_matrix2.txt.gz
perl ~/Data_Processing_Scripts/parse_series_matrix.pl Sasagawa_matrix2.txt > Sasagawa_Ann2.txt 2
mv ExprMat.txt Sasagawa_ExprMat2.txt
rm Sasagawa_matrix2.txt
