#!/usr/local/bin/bash
#bsub -R"select[mem>5000] rusage[mem=5000] span[hosts=1]" -M5000 -n 1 -q normal -o clean.%J.out -e clean.%J.err

source ~/.bashrc
conda activate gatk
export PATH=/lustre/scratch117/cellgen/team218/TA/LiverOrganoids/VariantCallingSoftware/gatk-4.1.2.0/:$PATH

BAMpath=/lustre/scratch117/cellgen/team218/TA/LiverOrganoids/ExomeCNVs/DNA_BAMS

REFfa=/lustre/scratch117/cellgen/team218/TA/genomebuilding/GATK/Homo_sapiens_assembly38.fasta


mutect_out=$1


stats_file=$mutect_out\.stats
tbi_file=$mutect_out\.tbi

mv $stats_file ${stats_file/.gz/}
stats_file=${stats_file/.gz/}

mv $tbi_file ${tbi_file/.gz/}
tbi_file=${tbi_file/.gz/}

gunzip $mutect_out
mutect_out=${mutect_out/.gz/}


gatk --java-options "-Xmx5G" FilterMutectCalls \
   -R $REFfa \
   -V $mutect_out \
   -O filtered.$mutect_out

gzip $mutect_out;
gzip filtered.$mutect_out;

# Intersect calls from each sample.
# read vcfs into R and combined columns then intersect.


# Step3: GetPileUpSummaries: next key step.
