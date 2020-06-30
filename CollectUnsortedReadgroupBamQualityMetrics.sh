#!/bin/bash
set -e

workdir=$1
pipeline=$2
sampleID=$3
laneName=$4

Workdir=$workdir/$sampleID
gatk \
  CollectMultipleMetrics \
  -I $Workdir/bwa/$laneName/$sampleID.$laneName.bam \
  -O $Workdir/bwa/$laneName/$sampleID.$laneName.readgroup \
  -AS true \
  -PROGRAM null \
  -PROGRAM CollectBaseDistributionByCycle \
  -PROGRAM CollectInsertSizeMetrics \
  -PROGRAM MeanQualityByCycle \
  -PROGRAM QualityScoreDistribution \
  -METRIC_ACCUMULATION_LEVEL null \
  -METRIC_ACCUMULATION_LEVEL ALL_READS
