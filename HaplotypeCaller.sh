#!/bin/bash
set -e

workdir=$1
pipeline=$2
sampleID=$3

ref_fasta=$pipeline/ref/hg19.fasta
interval_list=$pipeline/cfg/target_region.bed

Workdir=$workdir/$sampleID
output_bam_prefix=$Workdir/bwa/$sampleID
vcf_basename=$Workdir/gatk/$sampleID
mkdir -p $Workdir/gatk

gatk \
  --java-options "-Xms6000m -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10" \
  HaplotypeCaller \
  --tmp-dir=$Workdir/javatmp \
  -R ${ref_fasta} \
  -O ${vcf_basename}.g.vcf.gz \
  -I ${output_bam_prefix}.bqsr.bam \
  -L ${interval_list} \
  -ip 500 \
  -ERC GVCF \
  -contamination 0.0 \
  -G StandardAnnotation -G StandardHCAnnotation \
  -new-qual \
  -GQB 10 -GQB 20 -GQB 30 -GQB 40 -GQB 50 -GQB 60 -GQB 70 -GQB 80 -GQB 90 \
  --base-quality-score-threshold 6 \
  -DF MappingQualityReadFilter \
  --max-reads-per-alignment-start 0 \
  -mbq 10 \
  -stand-call-conf 30.0 \
  -OVM true \
  --showHidden
