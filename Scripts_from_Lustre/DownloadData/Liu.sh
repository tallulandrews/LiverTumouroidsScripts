#Liu et al. (2016) Identification of key factors conquering developmental arrest of somatic cell cloned embryos by combining embryo biopsy and single-cell sequencing. Cell Discovery. 2 : 16010
# The Gene Expression Omnibus (GEO) accession numbers for the single-cell RNA-seq, ULI-NChIP-seq datasets are GSE70605, GSE70606, and GEO number for the single-cell RRBS is GSE70607.
# PMID: 27462457

wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE70nnn/GSE70605/matrix/GSE70605_series_matrix.txt.gz
gunzip GSE70605_series_matrix.txt.gz
perl ~/Data_Processing_Scripts/parse_series_matrix.pl  GSE70605_series_matrix.txt 8 3 > Liu_anno.txt

