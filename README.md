# Large-scale alternative polyadenylation (APA)-wide association studies to identify putative susceptibility genes for cancers of the breast, ovary, prostate, colorectal, lung, and pancreas
---
* [Introduction](#Introduction)
* [Data Resource](#Resource)
* [APA-WAS Pipeline](#Pipeline)

<a name="Introduction"/>

# Introduction

In this study, we used RNA-sequencing (RNA-seq) data generated in multiple normal tissues, along with the matched whole genome sequencing (WGS) data generated in blood samples from the Genotype-Tissue Expression (GTEx), and large-scale GWAS data for cancers of breast, ovary, prostate, colorectum, lung, and pancreas to conduct APA-WAS to search for susceptibility genes and loci in these common cancers. 

Xingyi Guo1,2 †, *, Jie Ping1†, Yaohua Yang1†, Xiao-ou Shu1, Wanqing Wen1, Zhishan Chen1, Ran Tao3, Guochong Jia1, Jingni He4, Qiuyin Cai1, Quan Long4,5, Graham G Giles6,  Rachel Pearlman7,  Gad Rennert8, Pavel Vodicka9, Amanda Phipps10,11, Stephen B Gruber12, Graham Casey13, Ulrike Peters10,11, Jirong Long1, Wei Zheng1, *


# Data Resources

### R1. Whole genome sequencing (WGS) in blood samples and RNA sequencing (RNA-seq) data generated in normal tissues from the GTEx project (version 8).  

### R2. Summary statistics of GWAS data of European descendants for breast, ovary, prostate, colorectal, lung, and pancreatic cancer.



<a name="Pipeline"/>

# APA-WAS Pipeline 
We provided example codes and expected results for our reported candiate susceptiblity gene, KDSR, in breast cancer.

---

## Prerequisites
1. `Python (version >= 2.7.14)`
2. `bedtools v2.25.0-119-ga0dc5db`
3. `samtools v1.9`
4. `R v4.0.5`
5. `PEER https://github.com/PMBio/peer`
6. `USCS tools for extracting data v1.04.00`

## Step 1: Download BAM files of RNA-seq data in GTEX

To download RNA sequencing data in normal tissues from the GTEx project via Anvil Gen3 

```bash
check - gen3-client 
curl https://api.github.com/repos/uc-cdis/cdis-data-client/releases/latest | grep browser_download_url.*linux |  cut -d '"' -f 4 | wget -qi -
./gen3-client configure --profile=AnVIL --cred=credentials.json --apiendpoint=https://gen3.theanvil.io
./gen3-client download-multiple --profile=AnVIL --manifest=file-manifest.json --download-path=. --protocol=s3
```

## Step 2: Covert BAM files to BedGraph


bedtools v2.25.0-119-ga0dc5db

USCS tools for extracting data v1.04.00 (http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/)

```bash
genomeCoverageBed -bg -ibam ${SAMPLE}.bam -g ${CHR_SIZE_FILE} -split > ${SAMPLE}.bedgraph

### Extract KDSR gene (chr18:63327726-63367327)
grep chr18 ${SAMPLE}.bedgraph > ${SAMPLE}.chr18.bedgraph
bedSort ${SAMPLE}.chr18.bedgraph ${SAMPLE}.chr18.sorted.bedgraph
bedToBigBed ${SAMPLE}.chr18.sorted.bedgraph ${chromsize} ${SAMPLE}.chr18.bb
bigBedToBed -chrom=chr18 -start=63327726 -end=63367327 ${SAMPLE}.chr18.bb example/KDSR/${SAMPLE}.KDSR.bedgraph
rm ${SAMPLE}.chr18.bedgraph ${SAMPLE}.chr18.sorted.bedgraph ${SAMPLE}.chr18.bb
```

CHR_SIZE_FILE is downloaded from UCSC: http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/chromInfo.txt.gz

## Step 3: Calculate read depth for each sample

samtools v1.9

```bash
wigLoc=KDSR/${SAMPLE}.KDSR.bedgraph
echo "$wigLoc" >> ${SAMPLE}.depth
samtools view -c -F 260 ${SAMPLE}.bam >> ${SAMPLE}.depth
```

## Step 4: Generate BAM location with depth information

```bash
for f in $(ls *.depth)
do
    cat ${f} | awk '{printf $0"\t"}' | awk '{print $1"\t"$2}' >> example/breast.KDSR.mapping_bam_location_with_depth.txt
done
```

## Step 5: Generate DaPars2 configure file

```bash
#!/bin/sh

Annotated_3UTR=example/KDSR.3UTR.bed
outDir=example/
Coverage_threshold=10
Num_Threads=16

sequencing_depth_file=example/breast.KDSR.mapping_bam_location_with_depth.txt

wigfilelist=$(awk '{print $1}' ${sequencing_depth_file} | tr "\n" "," | sed 's/.$//' )

echo "Annotated_3UTR=${Annotated_3UTR}" > ${outDir}/KDSR.configure

echo "Aligned_Wig_files=${wigfilelist}" >> ${outDir}/KDSR.configure

echo "Output_directory=PDUI/" >> ${outDir}/KDSR.configure
echo "Output_result_file=PDUI" >> ${outDir}/KDSR.configure
echo "Coverage_threshold=${Coverage_threshold}" >> ${outDir}/KDSR.configure
echo "Num_Threads=${Num_Threads}" >> ${outDir}/KDSR.configure
echo "sequencing_depth_file=${sequencing_depth_file}" >> ${outDir}/KDSR.configure

```

## Step 6: Run DaPars2 to generate APA events and levels

Python v2.7.14

DaPars2 http://bioinfo.szbl.ac.cn/DaPars2/DaPars2.html

Running Time: less than 1 minute

```bash
ml load GCC/6.4.0-2.28 OpenMPI/2.1.1 Python/2.7.14 numpy/1.13.1-Python-2.7.14 scipy/0.19.1-Python-2.7.14 R/3.4.3
dapar2py=DaPars2/src/DaPars2_Multi_Sample_Multi_Chr.py

cd example

python2 ${dapar2py} KDSR.configure

cat example/PDUI_chr18/PDUI_result_temp.chr18.txt > example/KDSR.PDUI.txt


```

## Step 7: Perform Normalization of APA levels

R v4.0.5

```r
library(preprocessCore)

setwd("example")
PDUI <- read.table("KDSR.PDUI.txt", sep = "\t", header = T, row.names = 1, as.is = T)

PDUI <- PDUI[, -c(1:3)]
colnames(PDUI) <- sapply(colnames(PDUI), function(x) strsplit(strsplit(x, "KDSR\\.")[[1]][2], ".KDSR")[[1]][1])
colnames(PDUI) <- gsub("\\.", "-", colnames(PDUI))

na.no <- apply(PDUI, 1, function(x) length(which(is.na(x))))
PDUI_new <- PDUI[which(na.no < ncol(PDUI) * 0.5), ]

PDUI_QN <- normalize.quantiles(as.matrix(PDUI_new))
colnames(PDUI_QN) <- colnames(PDUI_new)
rownames(PDUI_QN) <- rownames(PDUI_new)	

PDUI_QN <- t(PDUI_QN)
PDUI_QN <- as.data.frame(cbind(rownames(PDUI_QN), PDUI_QN))
colnames(PDUI_QN)[1] <- "SUBJID"

write.table(PDUI_QN, file = "KDSR.PDUI.QN.txt", sep = "\t", quote = F, row.names = F, col.names = T, eol = "\n")

```

## Step 8: Calculate PEER factors for APA 

Please see detailed doucuments of PEER from here: https://github.com/PMBio/peer

```r
set.seed(1024)
library(peer)

setwd("/gpfs52/nobackup/sbcs/pingj2/APA_TWAS/example/")

PDUI_QN <- read.table("KDSR.PDUI.QN.txt", sep = "\t", header = T, row.names = 1, as.is = T)

PDUI_QN_complete <- PDUI_QN[, which(complete.cases(t(PDUI_QN)))]

model <- PEER()
PEER_setPhenoMean(model, as.matrix(PDUI_QN_complete))

peerN <- 15

if ( nrow(PEER_getPhenoMean(model)) < 150) {
    peerN <- 15
} else if ( nrow(PEER_getPhenoMean(model)) < 250 & nrow(PEER_getPhenoMean(model)) >= 150 ) {
    peerN <- 30
} else if ( nrow(PEER_getPhenoMean(model)) < 350 & nrow(PEER_getPhenoMean(model)) >= 250 ) {
    peerN <- 45
} else if ( nrow(PEER_getPhenoMean(model)) >= 350 ) {
    peerN <- 60
}

PEER_setNk(model, peerN)
PEER_update(model)

factors <- PEER_getX(model)

rownames(factors) <- rownames(PDUI_QN_complete)
colnames(factors) <- paste0("PEER", 1:ncol(factors))

save(factors, file = "KDSR.complete_data_PEER_factors.rda")
write.table(factors, file = "KDSR.complete_data_PEER_factors.txt", sep = "\t", row.names = T, quote = F)

pdf("Breast.complete_data_diagnostics_peer.pdf")
PEER_plotModel(model)
dev.off()
```

PEER factors calculated from the entire gene expression data.

## Step 9: Prepare covariates using genotype data from GTEx

Implemented in scrpits/5.Prepare_Covariates.R
   ```
    Rscript scrpits/Prepare_Covariates.R
  ```  

## Step 10: Run assocation analysis based on TWAS framework

### 10.1. Build model for the APA mark "NM_002035.4|KDSR|chr18|-" using buildmodel.r

#### 10.1.1 Input files
    1.1.1 qn_pdui.rds: quantile-normalized PDUI values of "NM_002035.4|KDSR|chr18|-".
    1.1.2 Breast_complete_data_PEER_factors.txt: PEER factors of PDUI values of all APA marks within breast tissue samples.
    1.1.3 Breast_covs.txt: Covariates in addition to PEER factors.
    1.1.4 dosage.rds: Alleic dosage of SNPs within 500Kb flanking region of "NM_002035.4|KDSR|chr18|-".
    1.1.5 snpinfo.rds: Information of SNPs in "dosage.rds", including chromosome, position(hg19), counted allele and other allele.
 
#### 10.1.1.2 Output files
    1.2.1 extra.csv: Information of prediction R-squared, p value, number of SNPs in model.
    1.2.2 weights.csv: Weigts, effect allele, other allele of SNPs in model.
    1.2.3 cov.txt.gz: Covariance of SNPs in model estimated using breast tissue samples.
    1.2.4 breast.db: R db file requested for SPrediXcan analyses.
 ```
Rscript scripts/buildmodel.r
  ```  


### 10.2. SPrediXcan, GCTA-COJO, and SPrediXcan after GCTA-COJO

#### 10.2.1 Input files for first SPrediXcan analyses
    10.2.1.1 breast.db: R db file generated by "buildmodel.r"
    
     10.2.1.2 gwas.csv: Summary statistics of SNPs in model, extracted from data downloaded from Zhang et al. Nat Genet 2020.
     
     10.2.1.3 cov.txt.gz: Covariance generated by "buildmodel.r"

#### 10.2.2 Input files for GCTA-COJO analyses
    10.2.2.1 /scratch/sbcs/yhyang/Refs/Chr18/chr18: Plink bfiles of chr18 of 1000 Genomes European participants
    
    10.2.2.2 rsid: rs number of SNPs in model and rs17743054, the GWAS-identified index SNP which is closest to "NM_002035.4|KDSR|chr     18|-". 18:60991409:AAAAG:A could not be found in 1000 Genomes data.
    
    10.2.2.3 predictingSNP_indexSNP.ss: summary statistics of SNPs in model and rs17743054
    
    10.2.2.4 cond.snp: rs17743054, whichi needs to be adjusted in GCTA-COJO analyses.

10.2.3 Output files
    10.2.3.1 assoc.csv: Results from first step SPrediXcan, i.e. Association of "NM_002035.4|KDSR|chr18|-" with breast cancer.
    
    10.2.3.2 cojores.cma.cojo: summary statistics of model SNPs conditioning on rs17743054.
    
    10.2.3.3 cojores.gwas: input file for second step SPrediXcan, generated according to cojores.cma.cojo.
    
    10.2.3.4 cojores.assoc.csv: Results from second step SPrediXcan, i.e. Association of "NM_002035.4|KDSR|chr18|-" with breast cancer    , conditioning on rs17743054.
  ```
  python scripts/gcta_cojo.py
  ```  

R scripts for each step can be found in the folder, "scricts/"; Example data for each step can be found in the folder, "Example/".

