#!/bin/bash
set -e

workdir=$1
pipeline=$2
sampleID=$3

Workdir=$workdir/$sampleID
vcf_basename=$Workdir/gatk/$sampleID

gatk \
  GenotypeGVCFs \
  --tmp-dir=$Workdir/javatmp \
  -R $ref \
  -O ${vcf_basename}.vcf.gz \
  -V ${vcf_basename}.g.vcf.gz \
  -L $interval_list -ip 500 \
  -G StandardAnnotation \
  -stand-call-conf 0.0 \
  --showHidden
