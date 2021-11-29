#!/bin/bash

MAXJOBS=384
bsub -J"mappingwithstararrayjob[1-$MAXJOBS]%10" -R"select[mem>2000] rusage[mem=2000]" -M2000 -q normal -o Convert.%J.%I /lustre/scratch117/cellgen/team218/TA/LiverOrganoids/Convert_CRAM_to_BAM.sh
