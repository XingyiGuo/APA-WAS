#!/bin/sh

gtexDir=
wgsDir=${gtexDir}/RootStudyConsentSet_phs000424.GTEx.v8.p2.c1.GRU/GenotypeFiles/phg001219.v1.GTEx_v8_WGS.genotype-calls-vcf.c1
plinkDir=

### Filter by PASS
bcftools view -f PASS ${wgsDir}/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_866Indiv.vcf.gz | bgzip -c > ${gtexDir}/VCF/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_866Indiv.PASS.vcf.gz

### extract all European samples with Whole blood and WGS data (n=660)
#EUR_FEMALE_BLOOD_WGS_IDs=${gtexDir}/GTEx_v8_EUR_FEMALE_WHOLE_BLOOD_WGS_IDs.txt

EUR_BLOOD_WGS_IDs=${gtexDir}/GTEx_v8_EUR_WHOLE_BLOOD_WGS_IDs.txt

bcftools view -Oz --threads 48 -S ${EUR_BLOOD_WGS_IDs} ${gtexDir}/VCF/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_866Indiv.PASS.vcf.gz > ${gtexDir}/VCF/GTEx_v8_EUR_Blood_WGS.vcf.gz

### filter SNPs

ml load GCC/5.4.0-2.26 VCFtools/0.1.14

vcftools --gzvcf ${gtexDir}/VCF/GTEx_v8_EUR_Blood_WGS.vcf.gz \
  --maf 0.05 \
  --hwe 0.0001 \
  --max-missing 0.05 \
  --remove-indels \
  --min-alleles 2 --max-alleles 2 \
  --plink-tped \
  --out ${gtexDir}/VCF/GTEx_v8_EUR_Blood_WGS.maf_0.05_missing_0.05_hwe_1e4_noINDEL

ml load PLINK/1.9b_5.2

plink --tfile ${gtexDir}/VCF/GTEx_v8_EUR_Blood_WGS.maf_0.05_missing_0.05_hwe_1e4_noINDEL \
  --out ${plinkDir}/GTEx_v8_EUR_Blood_WGS.maf_0.05_missing_0.05_hwe_1e4_noINDEL \
  --make-bed

plink --tfile ${gtexDir}/VCF/GTEx_v8_EUR_Blood_WGS.maf_0.05_missing_0.05_noINDEL \
  --out ${plinkDir}/GTEx_v8_EUR_Blood_WGS.maf_0.05_missing_0.05_noINDEL \
  --make-bed

## extract EUR female Samples

plink --bfile ${plinkDir}/GTEx_v8_EUR_Blood_WGS.maf_0.05_missing_0.05_hwe_1e4_noINDEL \
      --keep ${gtexDir}/GTEx_v8_EUR_FEMALE_WHOLE_BLOOD_WGS_IDs_plink.txt \
      --make-bed \
      --out ${plinkDir}/GTEx_v8_EUR_Female_Blood_WGS.maf_0.05_missing_0.05_hwe_1e4_noINDEL

plink --bfile ${plinkDir}/GTEx_v8_EUR_Female_Blood_WGS.maf_0.05_missing_0.05_hwe_1e4_noINDEL \
      --freq \
      --out ${plinkDir}/GTEx_v8_EUR_Female_Blood_WGS.allsnpsfromEUR.frq

## exclude MAF<0.05 in EUR female
awk '{if ($5 < 0.05) print $2}' ${plinkDir}/GTEx_v8_EUR_Female_Blood_WGS.allsnpsfromEUR.frq.frq > ${plinkDir}/GTEx_v8_EUR_Female_Blood_WGS.MAF0.05.snps

plink --bfile ${plinkDir}/GTEx_v8_EUR_Female_Blood_WGS.maf_0.05_missing_0.05_hwe_1e4_noINDEL \
      --exclude ${plinkDir}/GTEx_v8_EUR_Female_Blood_WGS.MAF0.05.snps \
      --make-bed \
      --out ${plinkDir}/GTEx_v8_EUR_Female_Blood_WGS_final

### Extract EUR male samples

plink --bfile ${plinkDir}/GTEx_v8_EUR_Blood_WGS.maf_0.05_missing_0.05_hwe_1e4_noINDEL \
      --remove ${gtexDir}/GTEx_v8_EUR_FEMALE_WHOLE_BLOOD_WGS_IDs_plink.txt \
      --make-bed \
      --out ${plinkDir}/GTEx_v8_EUR_male_Blood_WGS.maf_0.05_missing_0.05_hwe_1e4_noINDEL

plink --bfile ${plinkDir}/GTEx_v8_EUR_male_Blood_WGS.maf_0.05_missing_0.05_hwe_1e4_noINDEL \
      --freq \
      --out ${plinkDir}/GTEx_v8_EUR_male_Blood_WGS.allsnpsfromEUR.frq

## exclude MAF<0.05 in EUR male
awk '{if ($5 < 0.05) print $2}' ${plinkDir}/GTEx_v8_EUR_male_Blood_WGS.allsnpsfromEUR.frq.frq > ${plinkDir}/GTEx_v8_EUR_male_Blood_WGS.MAF0.05.snps

plink --bfile ${plinkDir}/GTEx_v8_EUR_Female_Blood_WGS.maf_0.05_missing_0.05_hwe_1e4 \
      --exclude ${plinkDir}/GTEx_v8_EUR_Female_Blood_WGS.MAF0.05.snps \
      --make-bed \
      --out ${plinkDir}/GTEx_v8_EUR_Female_Blood_WGS_final
