#!/bin/bash

INPUTDIR=$1
TMPDIR=$2
#STARPARAMFILE=$3
#ANNOTATIONgtf=$4

export REF_CACHE=/lustre/scratch117/cellgen/team218/TA/TemporaryFileDir
FILEStoMAP=($INPUTDIR/*.bam)
ARRAYINDEX=$(($LSB_JOBINDEX-1))
FILE=${FILEStoMAP[$ARRAYINDEX]}
SAMTOOLS=/nfs/users/nfs_t/ta6/RNASeqPipeline/software/CRAM/samtools-1.3.1/samtools
FASTQ1=$( basename $FILE )_1.fastq
FASTQ2=$( basename $FILE )_2.fastq


#Get upmapped reads and write to fastq
#$SAMTOOLS bam2fq -f 4 -1 $TMPDIR/$FASTQ1 -2 $TMPDIR/$FASTQ2 -n $FILE
$SAMTOOLS bam2fq -1 $TMPDIR/$FASTQ1 -2 $TMPDIR/$FASTQ2 -n $FILE

gzip $TMPDIR/$FASTQ1
gzip $TMPDIR/$FASTQ2

#INFQ1=$TMPDIR/$FASTQ1
#INFQ2=$TMPDIR/$FASTQ2

#GENOME=/lustre/scratch117/cellgen/team218/TA/STRIPED_GENOMES/ERCCs
#WORKINGDIR=/lustre/scratch117/cellgen/team218/TA/STRIPED_GENOMES/TemporaryFileDir/$LSB_JOBINDEX

#STAR=/nfs/users/nfs_t/ta6/RNASeqPipeline/software/STAR-STAR_2.4.0j/bin/Linux_x86_64_static/STAR
#if [[ $INFQ1 =~ \.gz$ ]] ; then
#    FILEnopath=`basename ${INFQ1%.fastq.gz}`
#    $STAR --runThreadN $NUMTHREADS --runMode alignReads --genomeDir $GENOME --readFilesIn $INFQ1 $INFQ2 --readFilesCommand zcat --parametersFiles $STARPARAMFILE --outFileNamePrefix $TMPDIR/$FILEnopath --outTmpDir $WORKINGDIR
#else
#    FILEnopath=`basename ${INFQ1%.fastq}`
#    $STAR --runThreadN $NUMTHREADS --runMode alignReads --genomeDir $GENOME --readFilesIn $INFQ1 $INFQ2 --parametersFiles $STARPARAMFILE --outFileNamePrefix $TMPDIR/$FILEnopath --outTmpDir $WORKINGDIR
#fi

#BAMFILE=""
#rm $INFQ1
#rm $INFQ2

##### Feature Counts #####
#featureCOUNT=/nfs/users/nfs_t/ta6/RNASeqPipeline/software/subread-1.4.6-p2-Linux-x86_64/bin/featureCounts
#OUTPUTFILE=$(basename "${INPUTBAM%.bam}.fragmentcounts")

#$featureCOUNT -p -T $NUMTHREADS -a $ANNOTATIONgtf -o $OUTDIR/$OUTPUTFILE $INPUTBAM # paired-end, no multimap
