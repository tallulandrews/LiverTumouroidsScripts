wget http://gastrulation.stemcells.cam.ac.uk/data/counts.gz
wget http://gastrulation.stemcells.cam.ac.uk/data/metadata.txt
wget http://gastrulation.stemcells.cam.ac.uk/data/tal1ExpCounts.gz
wget http://gastrulation.stemcells.cam.ac.uk/data/metadataTal1.txt

mv metadata.txt Gastrulation_metadata.txt
mv metadataTal1.txt Gastrulation_metadataTal1.txt
gunzip tal1ExpCounts.gz 
mv tal1ExpCounts Gastrulation_Tal1counts.txt
gunzip counts.gz 
mv counts Gastrulation_counts.txt
