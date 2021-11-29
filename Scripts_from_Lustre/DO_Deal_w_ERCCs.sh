#!/bin/bash

INPUTDIR=/lustre/scratch117/cellgen/team218/TA/LiverOrganoids/BAMS
INPUTFILES=($INPUTDIR/*.bam)
NUMFILES=${#INPUTFILES[@]}
MAXJOBS=$NUMFILES
OUTDIR=/lustre/scratch117/cellgen/team218/TA/LiverOrganoids/FastQs
mkdir -p $OUTDIR

bsub -J"featurecountsjobarray[1-$MAXJOBS]%50" -R"select[mem>1000] rusage[mem=1000]" -M1000 -q normal -o extract_unmapped.%J.%I  /lustre/scratch117/cellgen/team218/TA/LiverOrganoids/Deal_w_ERCCs.sh $INPUTDIR $OUTDIR
