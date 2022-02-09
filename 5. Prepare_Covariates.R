library(data.table)

gtexDir <- ""
peerDir <- ""
covDir <- ""
pduiDir <- ""

gtex_pheno <- fread(paste0(gtexDir, "phs000424.v8.pht002742.v8.p2.c1.GTEx_Subject_Phenotypes.GRU.txt"))

gtex_age_sex <- as.data.frame(gtex_pheno[, c("SUBJID", "AGE", "SEX")])

gtex_sample_attr <- fread(paste0(gtexDir, "GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt"))
gtex_sample_attr$ID <- sapply(gtex_sample_attr$SAMPID, function(x) paste0(strsplit(x, "-")[[1]][1], "-", strsplit(x, "-")[[1]][2]) )
gtex_sample_attr <- subset(gtex_sample_attr, SMAFRZE == "RNASEQ")

gtex_covs <- t(read.table(paste0(gtexDir, "GTEx.v8.all.covariates.txt"), sep = '\t', head = T, row.name = 1, as.is = T))
rownames(gtex_covs) <- gsub("\\.", "-", rownames(gtex_covs))

gtex_age_sex_platform <- merge(gtex_age_sex, gtex_covs, by.x = "SUBJID", by.y = "row.names")
gtex_age_sex_platform <- gtex_age_sex_platform[, c("SUBJID", "AGE", "SEX", "platform")]

tissues <- c("Breast", "Colon_Transverse", "Ovary", "Lung", "Pancreas", "Prostate")

for (n in 1:length(tissues)) {
    tissue <- tissues[n]
	print(tissue)
    racesex <- "EUR"

    if (tissue %in% c("Breast", "Ovary")) {
        racesex <- "EUR_Female"
    } else if (tissue == "Prostate") {
        racesex <- "EUR_Male"
    }

    if (tissue == "Breast") {
        tissueType <- unique(gtex_sample_attr$SMTSD)[9]
    } else if (tissue == "Colon_Transverse") {
        tissueType <- unique(gtex_sample_attr$SMTSD)[24]
    } else {
        tissueType <- tissue
    }

    genome_pc <- fread(paste0(gtexDir, racesex, "/GTEx_v8_", racesex, "_Blood_WGS_final.posID_hg19.eigenvec"))
    genome_pc <- genome_pc[, c(2:5)]
    colnames(genome_pc) <- c("SUBJID", "PC1", "PC2", "PC3")

    gtex_cov <- merge(gtex_age_sex_platform, genome_pc, by = "SUBJID")

    tissue_peer_factor <- fread(paste0(peerDir, tissue, "_complete_data_PEER_factors.txt"), header = F)
    colnames(tissue_peer_factor)[1] <- "SUBJID"
    colnames(tissue_peer_factor)[2:ncol(tissue_peer_factor)] <- paste0("PEER", 1:(ncol(tissue_peer_factor) - 1))
    tissue_peer_factor$SUBJID <- sapply(tissue_peer_factor$SUBJID, function(x) paste0(strsplit(x, "-")[[1]][1], "-", strsplit(x, "-")[[1]][2]))

    tissue_cov <- merge(gtex_cov, tissue_peer_factor, by = "SUBJID")

    tissue_sample_attr <- subset(gtex_sample_attr, SMTSD == tissueType)

    tissue_cov <- merge(tissue_sample_attr[, c("ID", "SMRIN")], tissue_cov, by.x = "ID", by.y = "SUBJID")

    save(tissue_cov, file = paste0(covDir, tissue, "_covs.Rdata"))
    fwrite(tissue_cov, file = paste0(covDir, tissue, "_covs.txt"), sep = "\t")

}
