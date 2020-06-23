#!/bin/bash
set -o pipefail
set -e 

workdir=$1
pipeline=$2
sampleID=$3
laneName=$4

ref_fasta=$pipeline/ref/hg19.fasta

Workdir=$workdir/$sampleID
mkdir -p $Workdir/bwa/$laneName
bwa \
    mem -K 100000000 -t 8 -M -Y \
    -R "@RG\tID:$sampleID\tSM:$sampleID\tLB:$laneName\tPL:COMPLETE" \
    ${ref_fasta} \
    $Workdir/filter/$laneName/$sampleID.$laneName.filter_1.fq.gz \
    $Workdir/filter/$laneName/$sampleID.$laneName.filter_2.fq.gz \
    | samtools view -S -b -o $Workdir/bwa/$laneName/$sampleID.$laneName.bam -
