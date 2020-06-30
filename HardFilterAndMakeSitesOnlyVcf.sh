#!/bin/bash
set -e

workdir=$1
pipeline=$2
sampleID=$3

excess_het_threshold = 54.69

Workdir=$workdir/$sampleID
vcf_basename=$Workdir/gatk/$sampleID

gatk \
  VariantFiltration \
  --filter-expression "ExcessHet > $excess_het_threshold" \
  --filter-name ExcessHet \
  -O ${vcf_basename}.variant_filtered.vcf.gz \
  -V ${vcf_basename}.vcf.gz \
  --showHidden
  
  gatk \
  MakeSitesOnlyVcf \
  -I $Workdir/gatk/$sampleID.variant_filtered.vcf.gz \
  -O $Workdir/gatk/$sampleID.sites_only.variant_filtered.vcf.gz \
  --showHidden
