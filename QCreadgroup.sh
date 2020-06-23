#!/bin/bash
set -e

workdir=$1
pipeline=$2
sampleID=$3

Workdir=$workdir/$sampleID
output_bam_prefix=$Workdir/bwa/$sampleID

gatk \
  CollectMultipleMetrics \
  -I ${output_bam_prefix}.merge.bam \
  -O ${output_bam_prefix}.readgroup \
  -AS true \
  -PROGRAM null \
  -PROGRAM CollectBaseDistributionByCycle \
  -PROGRAM CollectInsertSizeMetrics \
  -PROGRAM MeanQualityByCycle \
  -PROGRAM QualityScoreDistribution \
  -LEVEL null \
  -LEVEL ALL_READS \
  --showHidden=true 
