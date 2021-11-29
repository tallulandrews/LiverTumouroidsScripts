#!/bin/bash
TYPE=$1 # e.g. Donor-9-DM
ANNFILE=$2 # e.g. ../Annotation_table_Exp2.2.out

echo $TYPE

#Setup
export REF_CACHE=/lustre/scratch117/cellgen/team218/TA/TemporaryFileDir
SAMTOOLS=/nfs/users/nfs_t/ta6/RNASeqPipeline/software/CRAM/samtools-1.3.1/samtools
OUTDIR=/lustre/scratch117/cellgen/team218/TA/LiverOrganoids/BAMS

fileIDs=($(grep -P $TYPE $ANNFILE | awk -F'"' '{print $2}'))
 echo "grep -P $TYPE $ANNFILE"

echo "${fileIDs[1]}"


for file in "${fileIDs[@]}"; do
	# Find the file amongst all the raw data
	FILE=$(find /warehouse/team218_wh01/tallulah/Meritxell_Gurdon/Laura -name $file.cram)
	FILE1=$(basename $FILE)

	# Copy & Convert
	cp $FILE $OUTDIR/$FILE1
	$SAMTOOLS view -b -h $OUTDIR/$FILE1 -o $OUTDIR/$FILE1.bam

	# remove temporary file
	rm $OUTDIR/$FILE1

done
