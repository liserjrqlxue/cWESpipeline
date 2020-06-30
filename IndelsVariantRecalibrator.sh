#!/bin/bash
set -e

workdir=$1
pipeline=$2
sampleID=$3

Workdir=$workdir/$sampleID
vcf_basename=$Workdir/gatk/$sampleID

gatk \
  VariantRecalibrator \
  -V ${vcf_basename}.sites_only.variant_filtered.vcf.gz \
  -O ${vcf_basename}.indels.recal \
  --tranches-file ${vcf_basename}.indels.tranches \
  --trust-all-polymorphic \
  -tranche 100 \
  -tranche 99.95 \
  -tranche 99.9 \
  -tranche 99.5 \
  -tranche 99.0 \
  -tranche 97.0 \
  -tranche 96.0 \
  -tranche 95.0 \
  -tranche 94.0 \
  -tranche 93.5 \
  -tranche 93.0 \
  -tranche 92.0 \
  -tranche 91.0 \
  -tranche 90.0 \
  -an FS \
  -an ReadPosRankSum \
  -an MQRankSum \
  -an QD \
  -an SOR \
  -an DP \  
  -mode INDEL \
  --max-gaussians 4 \
  -resource mills,known=false,training=true,truth=true,prior=12:${mills_resource_vcf} \
  -resource axiomPoly,known=false,training=true,truth=false,prior=10:${axiomPoly_resource_vcf} \
  -resource dbsnp,known=true,training=false,truth=false,prior=2:${dbsnp_resource_vcf}
