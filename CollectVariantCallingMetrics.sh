#!/bin/bash
set -e

workdir=$1
pipeline=$2
sampleID=$3

Workdir=$workdir/$sampleID
vcf_basename=$Workdir/gatk/$sampleID

gatk \
  CollectVariantCallingMetrics \
  -I ${vcf_basename}.g.vcf.gz \
  -O ${vcf_basename}.g \
  --GVCF_INPUT true \
  --DBSNP ${dbsnp_resource_vcf} \
  -SD ${ref_dict} \
  --THREAD_COUNT 8 \
  -TI ${interval_list}
