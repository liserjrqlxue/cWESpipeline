#!/bin/bash
set -e

workdir=$1
pipeline=$2
sampleID=$3

Workdir=$workdir/$sampleID
output_bam_prefix=$Workdir/bwa/$sampleID
INPUT=""
for i in $Workdir/bwa/*/$sampleID.*.bam;do
    INPUT="$INPUT -I $i"
done

gatk \
  MergeSamFiles \
  --TMP_DIR $Workdir/javatmp \
  $INPUT \
  -O ${output_bam_prefix}.merge.bam \
  --USE_THREADING \
  -SO coordinate \
  --showHidden=true 
