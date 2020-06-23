#!/bin/bash
set -e

workdir=$1
pipeline=$2
sampleID=$3

Workdir=$workdir/$sampleID
output_bam_prefix=$Workdir/bwa/$sampleID

gatk \
  MarkDuplicates \
  -I ${output_bam_prefix}.merge.bam \
  -O ${output_bam_prefix}.dup.bam \
  -M ${output_bam_prefix}.dup.metrics \
  -ASO coordinate \
  -MAX_FILE_HANDLES 8000 \
  --MAX_OPTICAL_DUPLICATE_SET_SIZE 300000 \
  --OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500 \
  --CREATE_INDEX \
  --CLEAR_DT false \
  --showHidden
