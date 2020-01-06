source ~/.bashrc
conda activate gatk
export PATH=/lustre/scratch117/cellgen/team218/TA/LiverOrganoids/VariantCallingSoftware/gatk-4.1.2.0/:$PATH

BAMpath=/lustre/scratch117/cellgen/team218/TA/LiverOrganoids/ExomeCNVs/DNA_BAMS

REFfa=/lustre/scratch117/cellgen/team218/TA/genomebuilding/GATK/Homo_sapiens_assembly38.fasta

#hcc10_exome1=$BAMpath/SLX-13156.A007.HHMJHBBXX.Cut.QCed.Merged.bam.PD.bam.RG.bam.REASSIGNED.bam.GATKrealign.bam
#vcf=hcc10_exome1.vcf.gz
bam_file=$1
vcf=$2

bsub -R"select[mem>5000] rusage[mem=5000] span[hosts=1]" -M5000 -n 1 -q long -o vcf.%J.out -e vcf.%J.out \
gatk --java-options "-Xmx5G" Mutect2 \
   -R $REFfa \
   -I $bam_file \
   -O $vcf # Works, not fast...



# Step2: gatk FilterMutectCalls -R ref.fasta -V unfiltered.vcf -O filtered.vcf
# Step3: GetPileUpSummaries: next key step.
