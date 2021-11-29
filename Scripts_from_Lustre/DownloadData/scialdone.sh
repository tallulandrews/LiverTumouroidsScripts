wget http://gastrulation.stemcells.cam.ac.uk/data/counts.gz
wget http://gastrulation.stemcells.cam.ac.uk/data/metadata.txt
wget http://gastrulation.stemcells.cam.ac.uk/data/tal1ExpCounts.gz
wget http://gastrulation.stemcells.cam.ac.uk/data/metadataTal1.txt

mv metadata.txt Gastrulation_metadata.txt
mv metadataTal1.txt Gastrulation_metadataTal1.txt
tar -xvzf tal1ExpCounts.gz 
mv tal1ExpCounts.txt Gastrulation_Tal1counts.txt
tar -xvzf counts.gz 
mv counts.txt Gastrulation_counts.txt
