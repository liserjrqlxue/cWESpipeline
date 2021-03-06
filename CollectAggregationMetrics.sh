#!/bin/bash
set -e

workdir=$1
pipeline=$2
sampleID=$3

ref_fasta=$pipeline/ref/hg19.fasta

Workdir=$workdir/$sampleID
output_bam_prefix=$Workdir/bwa/$sampleID

gatk \
  CollectMultipleMetrics \
  -I ${output_bam_prefix}.bqsr.bam \
  -R ${ref_fasta} \
  -O ${output_bam_prefix} \
  -AS true \
  --PROGRAM null \
  --PROGRAM CollectAlignmentSummaryMetrics \
  --PROGRAM CollectInsertSizeMetrics \
  --PROGRAM CollectSequencingArtifactMetrics \
  --PROGRAM QualityScoreDistribution \
  --PROGRAM CollectGcBiasMetrics \
  -LEVEL null \
  -LEVEL SAMPLE \
  -LEVEL LIBRARY \
  --showHidden
