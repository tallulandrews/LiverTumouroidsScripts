featureCOUNT=/nfs/users/nfs_t/ta6/RNASeqPipeline/software/subread-1.4.6-p2-Linux-x86_64/bin/featureCounts
ANNOTATIONgtf=/lustre/scratch117/cellgen/team218/TA/genomebuilding/Homo_sapiens.GRCh38.79.gtf
NUMTHREADS=1
INDIR=/lustre/scratch117/cellgen/team218/TA/LiverOrganoids/ExomeCNVs/DNA_BAMS
OUTDIR=./

#yes multimap, single end, no quality threshold
#$featureCOUNT -O -M -T $NUMTHREADS -a $ANNOTATIONgtf -o $OUTDIR/SLX-13156.A016.HHMJHBBXX.Cut.QCed.Merged.bam.PD.bam.RG.bam.REASSIGNED.bam.GATKrealign.fragmentcounts $INDIR/SLX-13156.A016.HHMJHBBXX.Cut.QCed.Merged.bam.PD.bam.RG.bam.REASSIGNED.bam.GATKrealign.bam
$featureCOUNT -O -M -T $NUMTHREADS -a $ANNOTATIONgtf -g exon_id -o $OUTDIR/SLX-13156.A016.HHMJHBBXX.Cut.QCed.Merged.bam.PD.bam.RG.bam.REASSIGNED.bam.GATKrealign.exoncounts $INDIR/SLX-13156.A016.HHMJHBBXX.Cut.QCed.Merged.bam.PD.bam.RG.bam.REASSIGNED.bam.GATKrealign.bam
$featureCOUNT -O -M -T $NUMTHREADS -a $ANNOTATIONgtf -g exon_id -o $OUTDIR/SLX-13156.A012.HHMJHBBXX.Cut.QCed.Merged.bam.PD.bam.RG.bam.REASSIGNED.bam.GATKrealign.exoncounts $INDIR/SLX-13156.A012.HHMJHBBXX.Cut.QCed.Merged.bam.PD.bam.RG.bam.REASSIGNED.bam.GATKrealign.bam
$featureCOUNT -O -M -T $NUMTHREADS -a $ANNOTATIONgtf -g exon_id -o $OUTDIR/SLX-13156.A007.HHMJHBBXX.Cut.QCed.Merged.bam.PD.bam.RG.bam.REASSIGNED.bam.GATKrealign.exoncounts $INDIR/SLX-13156.A007.HHMJHBBXX.Cut.QCed.Merged.bam.PD.bam.RG.bam.REASSIGNED.bam.GATKrealign.bam
