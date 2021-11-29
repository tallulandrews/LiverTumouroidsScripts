ls ../BAMS/*.bam > all_bams.txt
#sort AllKeptCells.txt > kept_cells.txt
comm -23 all_bams.txt kept_cells.txt > bad_bam.txt
while read -r filename; do
	rm "$filename"
done <bad_bam.txt

rm bad_bam.txt
