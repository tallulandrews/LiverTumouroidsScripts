#!/bin/bash

export REF_CACHE=/lustre/scratch117/cellgen/team218/TA/TemporaryFileDir
SAMTOOLS=/nfs/users/nfs_t/ta6/RNASeqPipeline/software/CRAM/samtools-1.3.1/samtools

INPUTDIR=$1
OUTDIR=$2

FILEStoMAP=($INPUTDIR/*.bam)
ARRAYINDEX=$(($LSB_JOBINDEX-1))
INPUTBAM=${FILEStoMAP[$ARRAYINDEX]}
OUTFILE=$(basename $INPUTBAM).meta

$SAMTOOLS view -H $INPUTBAM > $OUTDIR/$OUTFILE
