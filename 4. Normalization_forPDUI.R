### Quantile Normalization
library(preprocessCore)

workDir <- ""
setwd(workDir)

tissueTypes <- c("Breast", "Colon_Transverse", "Lung", "Ovary", "Pancreas", "Prostate")

for (n in 1:length(tissueTypes)) {
	print(tissueTypes[n])
    PDUI <- read.table(paste0(tissueTypes[n], "_PDUI.txt"), sep = "\t", header = T, row.names = 1, as.is = T)

    PDUI <- PDUI[, -c(1:3)]
    colnames(PDUI) <- sapply(colnames(PDUI), function(x) strsplit(x, "BedGraph.")[[1]][2])
    colnames(PDUI) <- gsub("\\.", "-", colnames(PDUI))

	na.no <- apply(PDUI, 1, function(x) length(which(is.na(x))))
	PDUI_new <- PDUI[which(na.no < ncol(PDUI) * 0.5), ]

    PDUI_QN <- normalize.quantiles(as.matrix(PDUI_new))
    colnames(PDUI_QN) <- colnames(PDUI_new)
    rownames(PDUI_QN) <- rownames(PDUI_new)

	PDUI_QN <- t(PDUI_QN)
    PDUI_QN <- as.data.frame(cbind(rownames(PDUI_QN), PDUI_QN))
	colnames(PDUI_QN)[1] <- "SUBJID"

    write.table(PDUI_QN, file = paste0(tissueTypes[n], "_PDUI_QN.txt"), sep = "\t", quote = F, row.names = F, col.names = T, eol = "\n")
}

### PEER Normalization

###### singularity exec /gpfs23/scratch/sbcs/pingj2/r-peer_1.3--r341h470a237_1.sif R

set.seed(1024)
library(peer)

args <- commandArgs(trailingOnly = T)

workDir <- ""
peerDir <- ""

tissueTypes <- c("Breast", "Colon_Transverse", "Lung", "Ovary", "Pancreas", "Prostate")

n <- as.numeric(args[1])

print(tissueTypes[n])

PDUI_QN <- read.table(paste0(workDir, tissueTypes[n], "_PDUI_QN.txt"), sep = "\t", header = T, row.names = 1, as.is = T)

PDUI_QN_complete <- PDUI_QN[, which(complete.cases(t(PDUI_QN)))]

model <- PEER()
PEER_setPhenoMean(model, as.matrix(PDUI_QN_complete))

peerN <- 15
# the num here need to be decided per no. of subjects in the tissue of interest, see https://www.gtexportal.org/home/documentationPage
# for eQTL: 15 factors for N<150, 30 factors for 150<= N <250, 45 factors for 250<= N <350, and 60 factors for N>=350
# as a result of optimizing for the number of eGenes discovered
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

save(factors, file = paste0(peerDir, tissueTypes[n], "_complete_data_PEER_factors.rda"))
write.table(factors, file = paste0(peerDir, tissueTypes[n], "_complete_data_PEER_factors.txt"), sep = "\t", row.names = T, quote = F)

pdf(paste0(peerDir, tissueTypes[n], "_complete_data_diagnostics_peer.pdf"))
PEER_plotModel(model)
dev.off()
