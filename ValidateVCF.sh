
#!/bin/bash
set -e

workdir=$1
pipeline=$2
sampleID=$3

Workdir=$workdir/$sampleID
vcf_basename=$Workdir/gatk/$sampleID

gatk \
  ValidateVariants \
  -V ${vcf_basename}.g.vcf.gz \
  -R ${ref_fasta} \
  -L ${interval_list} -ip 500 \
  -gvcf \
  --validation-type-to-exclude ALLELES \
  --dbsnp ${dbsnp_resource_vcf}
