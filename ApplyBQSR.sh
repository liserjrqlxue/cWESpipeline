#!/bin/bash
set -e

workdir=$1
pipeline=$2
sampleID=$3

ref_fasta=$pipeline/ref/hg19.fasta

Workdir=$workdir/$sampleID
output_bam_prefix=$Workdir/bwa/$sampleID

gatk \
  ApplyBQSR \
  --tmp-dir=$Workdir/javatmp \
  -R ${ref_fasta} \
  -I ${output_bam_prefix}.dup.bam \
  -O ${output_bam_prefix}.bqsr.bam \
  -bqsr ${output_bam_prefix}.recal_data.grp \
  --create-output-bam-md5 \
  --create-output-bam-index \
  --static-quantized-quals 10 \
  --static-quantized-quals 20 \
  --static-quantized-quals 30 \
  --showHidden
