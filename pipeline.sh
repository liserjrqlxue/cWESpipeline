set -o pipefail
set -e

ref_fasta=$pipeline/hg19/hg19.fasta
ref_dict=$pipeline/hg19/hg19.dict

interval_list=$pipeline/config/cns_region_hg19_bychr/for500_region.bed
dbsnp_resource_vcf=$pipeline/hg19/dbsnp_138.hg19.vcf
mills_resource_vcf=$pipeline/hg19/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf
known_indels_vcf=$pipeline/hg19/Homo_sapiens_assembly19.known_indels_20120518.vcf
axiomPoly_resource_vcf=$pipeline/hg19/Axiom_Exome_Plus.genotypes.all_populations.poly.vcf.gz
hapmap_resource_vcf=$pipeline/hg19/hapmap_3.3.hg19.vcf.gz
omni_resource_vcf=$pipeline/hg19/1000G_omni2.5.hg19.vcf.gz
one_thousand_genomes_resource_vcf=$pipeline/hg19/1000G_phase1.snps.high_confidence.hg19.vcf.gz
excess_het_threshold = 54.69
indel_filter_level=99.7
snp_filter_level=99.7
tensor_type=reference
info_key=CNN_1D

Workdir=$workdir/$sampleID
vcf_basename=$Workdir/gatk/$sampleID
output_bam_prefix=$Workdir/bwa/$sampleID
export PATH=$pipeline/tools:$PATH

mkdir -p $Workdir
mkdir -p $Workdir/bwa
mkdir -p $Workdir/gatk
mkdir -p $Workdir/coverage

# Filter 每个lane一次，laneName不能重复（如可能重复需要把其他下机标签加入保持唯一）
mkdir -p $Workdir/filter/$laneName

SOAPnuke filter -o $Workdir/filter/$laneName \
    --fq1 $fq1 \
    --fq2 $fq2 \
    --adapter1 AAGTCGGAGGCCAAGCGGTCTTAGGAAGACAA \
    --adapter2 AAGTCGGATCGTAGCCATGTCGTTCTGTGAGCCAAGGAGTTG \
    --cleanFq1 $sampleID.$laneName.filter_1.fq.gz \
    --cleanFq2 $sampleID.$laneName.filter_2.fq.gz \
    --nRate 0.05 --lowQual 10  --seqType 0 --qualRate 0.5 -Q 2 -G 2

# Map reads to reference
# Alignment 每个lane一次，Filter后续步骤
mkdir -p $Workdir/bwa/$laneName

bwa \
    mem -K 100000000 -t 8 -M -Y \
    -R "@RG\tID:$sampleID\tSM:$sampleID\tLB:$laneName\tPL:COMPLETE" \
    ${ref_fasta} \
    $Workdir/filter/$laneName/$sampleID.$laneName.filter_1.fq.gz \
    $Workdir/filter/$laneName/$sampleID.$laneName.filter_2.fq.gz \
    | samtools view -S -b \
    -o $Workdir/bwa/$laneName/$sampleID.$laneName.bam  \
    -

# CollectUnsortedReadgroupBamQualityMetrics
# QC the aligned but unsorted readgroup BAM，每个lane一次，Alignment后续步骤，可与MergeSamFiles并行
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
  
# MergeSamFiles sam files and sort ，每个样品所有数据前序步骤alignment输出的所有bam作为输入-I
INPUT=""
for i in $Workdir/bwa/*/$sampleID.*.bam;do
    INPUT="$INPUT -I $i"
done
gatk \
  MergeSamFiles \
  --TMP_DIR $Workdir/javatmp \
  $INPUT \
  -O ${output_bam_prefix}.merge.bam \
  --USE_THREADING \
  -SO coordinate \
  --showHidden=true 

# QC the aligned but unsorted readgroup BAM ，额外QC步骤可以与MarkDuplicates并行
# Collect base quality and insert size metrics
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


# Mark duplicate reads to avoid counting non-independent observations
# MarkDuplicates
gatk \
  MarkDuplicates \
  -I ${output_bam_prefix}.merge.bam \
  -O ${output_bam_prefix}.dup.bam \
  -M ${output_bam_prefix}.dup.metrics \
  -ASO coordinate \
  -MAX_FILE_HANDLES 8000 \
  --MAX_OPTICAL_DUPLICATE_SET_SIZE 300000 \
  --OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500 \
  --CREATE_INDEX \
  --CLEAR_DT false \
  --showHidden

# Perform Base Quality Score Recalibration (BQSR) on the sorted BAM
# Generate Base Quality Score Recalibration (BQSR) model
# BQSR
gatk \
  BaseRecalibrator \
  --tmp-dir=$Workdir/javatmp \
  -I ${output_bam_prefix}.dup.bam \
  -O ${output_bam_prefix}.recal_data.grp \
  -OQ true \
  --known-sites $dbsnp_resource_vcf \
  --known-sites $mills_resource_vcf \
  --known-sites $known_indels_vcf \
  -R ${ref_fasta} \
  -L ${interval_list} -ip 500 \
  --showHidden

gatk \
  ApplyBQSR \
  --tmp-dir=$Workdir/javatmp \
  -R ${ref_fasta} \
  -I ${output_bam_prefix}.dup.bam \
  -O ${output_bam_prefix}.bqsr.bam \
  -bqsr ${output_bam_prefix}.recal_data.grp \
  -OQ true \
  --create-output-bam-md5 \
  --create-output-bam-index \
  --add-output-sam-program-record \
  --static-quantized-quals 10 \
  --static-quantized-quals 20 \
  --static-quantized-quals 30 \
  --showHidden

# bamdst QC ，bam QC可与HaplotypeCaller步骤并行
bamdst -p ${interval_list} --uncover 5 -o $Workdir/coverage ${output_bam_prefix}.bqsr.bam --cutoffdepth 20 

# AggregatedBamQC 

# QC the final BAM (consolidated after scattered BQSR) ，bam QC可与HaplotypeCaller步骤并行
# Collect alignment summary and GC bias quality metrics
# CollectReadgroupBamQualityMetrics
gatk \
  CollectMultipleMetrics \
  -I ${output_bam_prefix}.bqsr.bam \
  -R ${ref_fasta} \
  -O ${output_bam_prefix}.readgroup \
  -AS true \
  --PROGRAM null \
  --PROGRAM CollectAlignmentSummaryMetrics \
  --PROGRAM CollectGcBiasMetrics \
  -LEVEL null \
  -LEVEL READ_GROUP \
  --showHidden

# QC the final BAM some more (no such thing as too much QC) ，bam QC可与HaplotypeCaller步骤并行
# Collect quality metrics from the aggregated bam
# CollectAggregationMetrics
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

# Generate a checksum per readgroup in the final BAM ，bam QC可与HaplotypeCaller步骤并行
# Generate a checksum per readgroup
# CalculateReadGroupChecksum
gatk \
  CalculateReadGroupChecksum \
  -I ${output_bam_prefix}.bqsr.bam \
  -O ${output_bam_prefix}.bqsr.bam.read_group_md5

# HaplotypeCaller GATK4 
# --java-options 这个没特别考虑过，可能其他gatk步骤也要加，或者这个可以去掉，待深入研究确定

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
  -G StandardAnnotation -G StandardHCAnnotation -G AS_StandardAnnotation \
  -new-qual \
  -GQB 10 -GQB 20 -GQB 30 -GQB 40 -GQB 50 -GQB 60 -GQB 70 -GQB 80 -GQB 90 \
  --base-quality-score-threshold 6 \
  -DF MappingQualityReadFilter \
  --max-reads-per-alignment-start 0 \
  -mbq 10 \
  -stand-call-conf 0.0 \
  -OVM true \
  --showHidden
  
# Validate a gVCF with -gvcf specific validation
gatk \
  ValidateVariants \
  -V ${vcf_basename}.g.vcf.gz \
  -R ${ref_fasta} \
  -L ${interval_list} -ip 500 \
  -gvcf \
  --validation-type-to-exclude ALLELES \
  --dbsnp ${dbsnp_resource_vcf}

# Collect variant calling metrics from GVCF output
gatk \
  CollectVariantCallingMetrics \
  -I ${vcf_basename}.g.vcf.gz \
  -O ${vcf_basename}.g \
  --GVCF_INPUT true \
  --DBSNP ${dbsnp_resource_vcf} \
  -SD ${ref_dict} \
  --THREAD_COUNT 8 \
  -TI ${interval_list}

# GenotypeGVCFs 
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

# HardFilterAndMakeSitesOnlyVcf
gatk \
  VariantFiltration \
  --filter-expression "ExcessHet > $excess_het_threshold" \
  --filter-name ExcessHet \
  -O ${vcf_basename}.variant_filtered.vcf.gz \
  -V ${vcf_basename}.vcf.gz \
  --showHidden

gatk \
  MakeSitesOnlyVcf \
  -I ${vcf_basename}.variant_filtered.vcf.gz \
  -O ${vcf_basename}.sites_only.variant_filtered.vcf.gz \
  --showHidden

# IndelsVariantRecalibrator 
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

# VariantRecalibrator
gatk \
  VariantRecalibrator \
  -V ${vcf_basename}.sites_only.variant_filtered.vcf.gz \
  -O ${vcf_basename}.snps.recal \
  --tranches-file ${vcf_basename}.snps.tranches \
  --trust-all-polymorphic \
  -tranche 100 \
  -tranche 99.95 \
  -tranche 99.9 \
  -tranche 99.8 \
  -tranche 99.6 \
  -tranche 99.5 \
  -tranche 99.4 \
  -tranche 99.3 \
  -tranche 99.0 \
  -tranche 98.0 \
  -tranche 97.0 \
  -tranche 90.0 \
  -an QD \
  -an MQRankSum \
  -an ReadPosRankSum \
  -an FS \
  -an MQ \
  -an SOR \
  -an DP \  
  -mode SNP \
  --max-gaussians 6 \
  -resource hapmap,known=false,training=true,truth=true,prior=15:${hapmap_resource_vcf} \
  -resource omni,known=false,training=true,truth=true,prior=12:${omni_resource_vcf} \
  -resource 1000G,known=false,training=true,truth=false,prior=10:${one_thousand_genomes_resource_vcf} \
  -resource dbsnp,known=true,training=false,truth=false,prior=7:${dbsnp_resource_vcf}

# ApplyRecalibration 
gatk \
  ApplyVQSR \
  -O ${vcf_basename}.indel.recalibrated.vcf \
  -V ${vcf_basename}.variant_filtered.vcf.gz \
  --recal-file ${vcf_basename}.indels.recal \
  --tranches-file ${vcf_basename}.indels.tranches \
  --truth-sensitivity-filter-level ${indel_filter_level} \
  --create-output-variant-index true \
  -mode INDEL
  
gatk \
  ApplyVQSR \
  -O ${vcf_basename}.filtered.vcf.gz \
  -V ${vcf_basename}.indel.recalibrated.vcf \
  --recal-file ${vcf_basename}.snps.recal \
  --tranches-file ${vcf_basename}.snps.tranches \
  --truth-sensitivity-filter-level ${snp_filter_level} \
  --create-output-variant-index true \
  -mode SNP

# CollectVariantCallingMetrics 
gatk \
  CollectVariantCallingMetrics \
  -I ${vcf_basename}.filtered.vcf.gz \
  -O ${vcf_basename} \
  --DBSNP ${dbsnp_resource_vcf} \
  -SD ${ref_dict} \
  --THREAD_COUNT 8 \
  -TI ${interval_list}

# 另一种过滤
# HardFilterVcf
gatk \
  VariantFiltration \
  -V ${vcf_basename}.vcf.gz\
  -L ${interval_list} \
  --filter-expression "QD < 2.0 || FS > 30.0 || SOR > 3.0 || MQ < 40.0 || MQRankSum < -3.0 || ReadPosRankSum < -3.0" \
  --filter-name "HardFiltered" \
  -O ${vcf_basename}.SS.filtered.vcf.gz 


# CNNScoreVariants
gatk \
  CNNScoreVariants \
  -V ${vcf_basename}.SS.filtered.vcf.gz \
  -R ${ref_fasta} \
  -O ${vcf_basename}.SS.scored.vcf.gz \
  -tensor-type ${tensor_type}


# FilterVariantTranches
gatk \
  FilterVariantTranches \
  -V ${vcf_basename}.SS.scored.vcf.gz \
  -O ${vcf_basename}.SS.CNNS.filtered.vcf.gz \
  --snp-tranche 100 \
  --snp-tranche 99.95 \
  --snp-tranche 99.9 \
  --snp-tranche 99.8 \
  --snp-tranche 99.6 \
  --snp-tranche 99.5 \
  --snp-tranche 99.4 \
  --snp-tranche 99.3 \
  --snp-tranche 99.0 \
  --snp-tranche 98.0 \
  --snp-tranche 97.0 \
  --snp-tranche 90.0 \
  --indel-tranche 100 \
  --indel-tranche 99.95 \
  --indel-tranche 99.9 \
  --indel-tranche 99.5 \
  --indel-tranche 99.0 \
  --indel-tranche 97.0 \
  --indel-tranche 96.0 \
  --indel-tranche 95.0 \
  --indel-tranche 94.0 \
  --indel-tranche 93.5 \
  --indel-tranche 93.0 \
  --indel-tranche 92.0 \
  --indel-tranche 91.0 \
  --indel-tranche 90.0 \
  --resource ${hapmap_resource_vcf} \
  --resource ${omni_resource_vcf} \
  --resource ${one_thousand_genomes_resource_vcf} \
  --resource ${dbsnp_resource_vcf} \
  --info-key ${info_key} \
  --create-output-variant-index true
