#!/bin/bash

INPUTDIR=/lustre/scratch117/cellgen/team218/TA/LiverOrganoids/BAMS
INPUTFILES=($INPUTDIR/*.bam)
NUMFILES=${#INPUTFILES[@]}
MAXJOBS=$NUMFILES
OUTDIR=/lustre/scratch117/cellgen/team218/TA/LiverOrganoids/Tmp
mkdir -p $OUTDIR

bsub -J"featurecountsjobarray[1-$MAXJOBS]%50" -R"select[mem>1000] rusage[mem=1000]" -M1000 -q normal -o extract_meta.%J.%I  /lustre/scratch117/cellgen/team218/TA/LiverOrganoids/Extract_Metadata.sh $INPUTDIR $OUTDIR
