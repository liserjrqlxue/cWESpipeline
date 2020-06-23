#!/bin/bash
set -e

workdir=$1
pipeline=$2
sampleID=$3

ref_fasta=$pipeline/ref/hg19.fasta
interval_list=$pipeline/cfg/target_region.bed
dbSNP=$pieline/bundle/dbsnp_138.hg19.vcf.gz
mills_resource=$pipeline/bundle/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf.gz

Workdir=$workdir/$sampleID
output_bam_prefix=$Workdir/bwa/$sampleID


gatk \
  BaseRecalibrator \
  --tmp-dir=$Workdir/javatmp \
  -I ${output_bam_prefix}.dup.bam \
  -O ${output_bam_prefix}.recal_data.grp \
  -OQ \
  --known-sites ${dbSNP} \
  --known-sites ${mills_resource} \
  -R ${ref_fasta} \
  -L ${interval_list} \
  -ip 500 \
  --showHidden
