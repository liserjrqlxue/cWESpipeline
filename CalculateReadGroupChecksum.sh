#!/bin/bash
set -e

workdir=$1
pipeline=$2
sampleID=$3

Workdir=$workdir/$sampleID
output_bam_prefix=$Workdir/bwa/$sampleID

gatk \
  CalculateReadGroupChecksum \
  -I ${output_bam_prefix}.bqsr.bam \
  -O ${output_bam_prefix}.bqsr.bam.read_group_md5
