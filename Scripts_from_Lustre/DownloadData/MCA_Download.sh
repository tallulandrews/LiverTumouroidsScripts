wget  https://ndownloader.figshare.com/files/10233567?private_link=865e694ad06d5857db4b -O MCA_raw.rar
#wget https://ndownloader.figshare.com/files/10351110?private_link=865e694ad06d5857db4b
wget https://ndownloader.figshare.com/files/10514584?private_link=865e694ad06d5857db4b -O MCA_All_batch_corrected.zip
unzip MCA_All_batch_corrected.zip
unrar x MCA_raw.rar

unzip MCA_DGE/2XSequenced_mES_CJ7.dge.zip
gzip MCA_DGE/mES_CJ7-3_600cells.dge.txt

#wget https://ndownloader.figshare.com/files/10351110 -O Annotation.tar.gz

#wget http://www.cell.com/cms/attachment/2119021449/2088127046/mmc4.xlsx -O Figure2_Markers_Anno_88.xlsx
#wget http://www.cell.com/cms/attachment/2119021449/2088127049/mmc5.xlsx -O Figure2_Markers_Anno_Each.xlsx

# Annotations from Angela
wget https://raw.githubusercontent.com/czi-hca-comp-tools/easy-data/MCA/datasets/mouse-cell-atlas/MCA-cell-ann.csv
