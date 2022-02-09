# module load PLINK/1.9b_5.2 GCC/8.2.0 OpenMPI/3.1.4 R/3.6.0

library(glmnet)
set.seed(1024)

### config
args <- commandArgs(trailingOnly = TRUE)

tissue <- args[1] ## Breast
targetGene <- args[2] ### test: NM_000699.4|AMY2A|chr1|+ / NR_146324.2|MIB2|chr1|+ / NM_002117.6|HLA-C|chr6|-
geneid <- paste0(strsplit(targetGene, "\\|")[[1]][1], "_", strsplit(targetGene, "\\|")[[1]][2])

plink_path <- 'plink'

SNPs_in_sumstat <- ""

workDir <- paste0("APA_TWAS/3aQTL/", tissue)
setwd(workDir)

tmp_folder <- paste0(workDir, "/tmp/")
exp_folder <- paste0(workDir, "/exp/")
res_folder <- paste0(workDir, "/res/")

covs <- read.table(paste0("APA_TWAS/data/covs/", tissue, "_covs.txt"), sep = "\t", header = 1, as.is = T)

racesex <- "EUR"

if (tissue %in% c("Breast", "Ovary")) {
  racesex <- "EUR_Female"
} else if (tissue == "Prostate") {
  racesex <- "EUR_Male"
}

genotype_path <- paste0("GTEx_v8_plink/", racesex, "/GTEx_v8_", racesex, "_Blood_WGS_final")

### Some functions ###
loadRData <- function(fileName){
  load(fileName)
  get(ls()[ls() != "fileName"])
}

# Load APA annotation file
cat(' INFO loading APA position annotation ...\n')
apa_info <- read.table("hg38_refseq_extracted_3UTR.bed", header = F, sep = "\t", as.is = T)

# mkdir tmp folder. will be cleaned
cat(' INFO mkdir tmp folders ...\n')
options(warn = -1)
dir.create(tmp_folder)
dir.create(res_folder)

dir.create(paste0(tmp_folder, '/', tissue, '_', geneid))
tmp_folder <- paste0(tmp_folder,'/', tissue, '_', geneid)
options(warn = 0)

# Get chr pos for the gene
chr <- as.numeric(sub('^...', '', apa_info$V1[which(apa_info$V4 == targetGene)]))
pos_from <- apa_info$V2[which(apa_info$V4 == targetGene)]
pos_to <- apa_info$V3[which(apa_info$V4 == targetGene)]
pos_from_500K <- max(1, pos_from - 500000)
pos_to_500K <- pos_to + 500000

# Load expression
cat(' INFO loading expression data ...\n')
exptab <- loadRData(paste0("APA_TWAS/data/PDUI/", tissue, "_PDUI_QN.Rdata"))
exp.all <- as.data.frame(t(t(exptab[, targetGene])))
rownames(exp.all) <- rownames(exptab)
colnames(exp.all) <- "PDUI"

exp <- t(t(exp.all[which(!is.na(exp.all$PDUI)), ]))
rownames(exp) <- rownames(exp.all)[which(!is.na(exp.all$PDUI))]
colnames(exp) <- "PDUI"

exp <- apply(exp, 2, function(x) qnorm(rank(x, ties.method = "r") / (length(x) + 1)) )

## Calculate residual for model

exp_cov <- merge(exp, covs, by.x = 'row.names', by.y = 'ID')
colnames(exp_cov)[1] <- "SUBJID"

# Extract genotypes from plink file to dosage file (500Kb)
cat(' INFO generating dosage genotype data ...\n')

dosagecmd <- paste0(plink_path,' --bfile ', genotype_path, ' --extract ', SNPs_in_sumstat, ' --chr ', chr, ' --from-bp ', pos_from_500K, ' --to-bp ', pos_to_500K, ' --recode A --out ', tmp_folder, '/', geneid)
system(dosagecmd, ignore.stdout = T, ignore.stderr = T)
bedcmd <- paste0(plink_path,' --bfile ', genotype_path, ' --extract ', SNPs_in_sumstat, ' --chr ', chr, ' --from-bp ', pos_from_500K, ' --to-bp ', pos_to_500K, ' --make-bed --out ', tmp_folder, '/', geneid)
system(bedcmd, ignore.stdout = T, ignore.stderr = T)

# Load dosage file (500Kb)
dosage_500K <- try(read.table(paste0(tmp_folder, '/', geneid, '.raw'), header = T, stringsAsFactors = F))

if ('try-error' %in% class(dosage_500K)) {
  stop(paste0('no SNP available for ', targetGene))
}

dosage_500K <- dosage_500K[, -c(1, 3:6)]
colnames(dosage_500K) <- c('SUBJID', sapply(colnames(dosage_500K)[-1], function(x) strsplit(x, "[_]")[[1]][1])) #rm the counted allele from rsnum
dosage_500K[, -1] <- round(apply(dosage_500K[, -1], 2, function(x) ifelse(is.na(x), mean(x, na.rm = T), x)), 3) #post imputation imputation. NA replaced by mean

# Load Allele Info (500Kb)
snp_info_500K <- read.table(paste0(tmp_folder, '/', geneid, '.bim'), stringsAsFactors = F)
snp_info_500K$counted_allele <- snp_info_500K$V5
snp_info_500K$ref_allele <- snp_info_500K$V6
snp_info_500K$chr_bp <- paste0(snp_info_500K$V1, '_', snp_info_500K$V4)
colnames(snp_info_500K)[c(2,4)] <- c('rsid', 'bp')
snp_info_500K <- snp_info_500K[, c('rsid', 'chr_bp', 'bp', 'ref_allele', 'counted_allele')]

# 3aQTL model

#fit single tissue model to get proper window size and a lambda range
cat(' INFO fitting signle tissue prediction model \n')

exp_cov_dosage <- merge(exp_cov, dosage_500K, by = 'SUBJID')

aqtl <- t(sapply(colnames(dosage_500K)[-1], function(x) {
  lmfit <- lm(as.formula(paste0(colnames(exp_cov)[2], " ~ ", paste(colnames(exp_cov)[-c(1:2)], collapse = " + "), " + ", x)), data = exp_cov_dosage)
  return(summary(lmfit)$coef[nrow(summary(lmfit)$coef), ])
}))

aqtl <- as.data.frame(cbind(targetGene, rownames(aqtl), aqtl))
colnames(aqtl)[2] <- "SNP"

write.table(aqtl, file = paste0(res_folder, tissue, "_3aQTL_", geneid, ".txt"), sep = "\t", row.name = F, col.names = T, quote = F)

#cleaning
cat(' INFO cleaning the tmp folder ... \n')
cmd = paste0('rm -r ', tmp_folder) #will only clean the subfolder under the tmp folder
system(cmd, wait = T)

cat(' INFO done \n')
