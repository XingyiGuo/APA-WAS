### Generate residual for modeling

rm(list = ls())
pacman::p_load(GenABEL, preprocessCore)
options("na.action")

args <- commandArgs(trailingOnly = TRUE)
pdui <- read.table("Prostate_PDUI.txt", header = T, row.names = 1, check.names = F, stringsAsFactors = F)
peer <- read.table("Prostate_complete_data_PEER_factors.txt", header = T, row.names = 1, check.names = F, stringsAsFactors = F) 
cov <- read.table("Prostate_covs.txt", header = T, row.names = 1, check.names = F, stringsAsFactors = F)
ind <- read.table("hg38.hg19.txt", header = F, stringsAsFactors = F)
rownames(ind) <- ind$V1

loci <- data.frame(cbind(rownames(pdui), pdui[,3]), stringsAsFactors = F)
colnames(loci) <- c("Gene", "Loci")
rownames(loci) <- loci$Loci
loci$hg19 <- ind[rownames(loci), ]$V2

write.table(file = 'gene.hg38.hg19', loci, sep = "\t", quote = F, col.names = T, row.names = F)

pdui <- pdui[, 4:ncol(pdui)]
cname <- c()
for (i in 1:ncol(pdui)) {
    cname <- c(cname, paste('GTEX', strsplit(colnames(pdui)[i], '-')[[1]][2], sep = '-'))
}
colnames(pdui) <- cname
rname <- c()
for (i in 1:nrow(peer)) {
    rname <- c(rname, paste('GTEX', strsplit(rownames(peer)[i], '-')[[1]][2], sep = '-'))
}
rownames(peer) <- rname
cov <- cov[,1:7]
cov <- cov[,-3]
tmp_cov <- apply(as.matrix(cov), 2, as.numeric)
rownames(tmp_cov) <- rownames(cov)
cov <- tmp_cov
tmp_peer <- apply(as.matrix(peer), 2, as.numeric)
rownames(tmp_peer) <- rownames(peer)
peer <- tmp_peer

qc <- rep(1,nrow(pdui))
for (i in 1:nrow(pdui)) { 
    if ((length(which(is.na(pdui[i,])))/ncol(pdui)) > 0.05) {
        qc[i] <- 0
    }
}
pdui <- pdui[which(qc == 1),]
qn_pdui <- normalize.quantiles(as.matrix(pdui))
colnames(qn_pdui) <- colnames(pdui)
rownames(qn_pdui) <- rownames(pdui)

table(colnames(qn_pdui) == rownames(cov))
table(colnames(qn_pdui) == rownames(peer))

output <- c("Gene", colnames(qn_pdui))
sp <- c("Gene", "Pval")
for (i in 1:nrow(qn_pdui)) {
  skip_to_next <- FALSE
  tryCatch(
    { rn_pdui <- as.numeric(rntransform(as.numeric(qn_pdui[i, ]))) },
    error = function(e) { skip_to_next <<- TRUE } )
  if (skip_to_next) next  
  tmp_d <- data.frame(cbind(rn_pdui, cov, peer[, 1:30]), stringsAsFactors = F)
  res <- residuals(lm(rn_pdui ~ ., tmp_d, na.action = "na.exclude"))
  rn_res <- as.numeric(rntransform(res))
  output <- rbind(output, c(rownames(qn_pdui)[i], rn_res))
}
colnames(output) <- output[1, ]
output <- output[-1, ]
write.table(file = 'pdui.res', output, sep = "\t", quote = F, col.names = T, row.names = F)

### Model building

rm(list = ls())
pacman::p_load(glmnet, foreach, doParallel)
args <- commandArgs(trailingOnly=TRUE)

### Genotype
gt_fd <- 'Geno/Prostate/'
bim_f <- paste(gt_fd, 'Chr', args[1], '/chr', args[1], '.bim', sep='')
fam_f <- paste(gt_fd, 'Chr', args[1], '/chr', args[1], '.fam', sep='')
gt_f <- paste(gt_fd, 'Chr', args[1], '/chr', args[1], '.traw', sep='')
bim <- read.table(bim_f, sep="\t", header=F, stringsAsFactors=F)
rownames(bim) <- bim$V2
fam <- read.table(fam_f, sep="\t", header=F, stringsAsFactors=F)
rownames(fam) <- fam$V1
gt <- read.table(gt_f, header=T, stringsAsFactors=F, row.names=2, check.names=F)
snp_info <- gt[, c("CHR","POS","COUNTED","ALT")]
gt <- data.matrix(t(gt[,6:ncol(gt)]))
rownames(gt) <- unlist(strsplit(rownames(gt), "_"))[2*(1:nrow(gt))-1]
gt <- apply(gt, 2, function(x) ifelse(is.na(x),mean(x,na.rm=T),x))

### PDUI
pdui <- read.table('pdui.res', header=T, row.names=1, stringsAsFactors=F, check.names=F)
gt <- gt[colnames(pdui),]
mycores <- detectCores() - 20
mycores <- 64

### Pos
pos <- read.table('gene.hg38.hg19', header=T, row.names=1, stringsAsFactors=F, check.names=F)
pos <- pos[rownames(pdui),]
chrpat <- paste("chr",args[1],":",sep='')
pos <- pos[(which(grepl(chrpat, pos$hg19))),]
pdui <- pdui[rownames(pos),]

ext_head <- c("Gene", "cvm", "lambda.iter", "lambda.min", "n.snps", "R2", "Pval")
wt_head <- c("Gene", "SNP", "EffectA", "RefA", "beta")
cov_head <- c('Gene', 'RSID2','RSID1','VALUE')

func <- function(x, m = pdui, g = gt, mc = mycores) {
  every <- ceiling(nrow(m)/mc)
  st <- (x-1) * every + 1
  en <- st + every - 1
  if (en > nrow(m)) { en <- nrow(m) }
  extra_file <- paste('chr', args[1], '_', x, '.extra.txt', sep=''); write(ext_head,file=extra_file,ncolumns=7,sep="\t")
  weights_file <- paste('chr', args[1], '_', x, '.weights.txt', sep=''); write(wt_head,file=weights_file,ncolumns=5,sep="\t")
  cov_file <- paste('chr', args[1], '_', x, '.cov.txt', sep=''); write(cov_head,file=cov_file,ncolumns=4,sep="\t")

  for (j in st:en) {
    pdui_val <- as.numeric(m[j,])
    pdui_val <- scale(pdui_val, center=T, scale=T)
    pdui_val[is.na(pdui_val)] <- 0
    rownames(pdui_val) <- colnames(m) 
    left <- as.numeric(strsplit(strsplit(pos[rownames(pdui)[j],]$hg19,':')[[1]][2],'-')[[1]][1])-500000
    right <- as.numeric(strsplit(strsplit(pos[rownames(pdui)[j],]$hg19,':')[[1]][2],'-')[[1]][2])+500000
    snps <- subset(bim, ((V4 >= left) & (V4 <= right)))$V2
    geno <- g[,snps]
    #minorsnps <- subset(colMeans(geno), colMeans(geno,na.rm=TRUE)>0)
    #minorsnps <- names(minorsnps); geno <- geno[,minorsnps]
    if (is.null(dim(geno)) | (!is.null(dim(geno)) && dim(geno)[2] == 0)) next
    set.seed(1024)
    fit <- cv.glmnet(geno, pdui_val, nfolds = 10, alpha = 0.5, keep = T, parallel = F)
    fit.df <- data.frame(fit$cvm, fit$lambda,1:length(fit$cvm))
    best.lam <- fit.df[which.min(fit.df[,1]),]
    cvm.best <- best.lam[,1]
    lambda.best <- best.lam[,2]
    nrow.best <- best.lam[,3]
    ret <- as.data.frame(fit$glmnet.fit$beta[,nrow.best]) # get betas from best lambda
    ret[ret == 0.0] <- NA
    bestbetas <- as.vector(ret[which(!is.na(ret)),]) # vector of non-zero betas
    names(bestbetas) <- rownames(ret)[which(!is.na(ret))]
    pred.mat <- fit$fit.preval[,nrow.best] # pull out predictions at best lambda
    if(length(bestbetas) > 0) {
      res <- cor.test(pdui_val, pred.mat)
      rval <- res$estimate[[1]]
      pval <- res$p.value
      if (rval > 0.1 & pval < 0.05) {
        sum_res <- c(rownames(m)[j], cvm.best, nrow.best, lambda.best, length(bestbetas), rval**2, pval)
        write(sum_res,file = extra_file,ncolumns=8,append=T,sep="\t")
      
        bestbetalist <- names(bestbetas)
        bestbetainfo <- snp_info[bestbetalist,]
        betatable <- as.matrix(cbind(bestbetainfo,bestbetas))
        betafile <- cbind(rownames(m)[j],rownames(betatable),betatable[,3],betatable[,4],betatable[,5]) ##output "gene","SNP","effectAllele","refAllele","beta"
        write(t(betafile), file = weights_file, ncolumns = 5, append = T, sep = "\t") # t() necessary for correct output from write() function

        if (length(bestbetalist)>1) {
          dsg <- g[,bestbetalist]
          cov <- cov(as.matrix(dsg))
          cov[upper.tri(cov)] <- NA
          cov <- cbind(expand.grid(dimnames(cov)), value = as.vector(cov))
          colnames(cov) <- c('RSID2','RSID1','VALUE')
          cov <- cov[!is.na(cov$VALUE),]
          cov$GENE <- rownames(m)[j]
          cov = cov[,c('GENE','RSID1','RSID2','VALUE')]
          write(t(cov),file=cov_file,ncolumns=4,append=T,sep="\t")
        } 
      }
    }
  }
}
cl <- makeCluster(mycores)
registerDoParallel(cl)
x <- foreach(x=1:mycores,.packages='glmnet') %dopar% func(x)
stopCluster(cl)

### Make DB file for association analysis ###

library(sqldf)

db <- dbConnect(SQLite(), dbname="prst.db")

dbWriteTable(conn = db, name = "construction", value = "construction.csv", row.names = FALSE, header = TRUE)
dbWriteTable(conn = db, name = "extra", value = "extra.csv", row.names = FALSE, header = TRUE)
dbWriteTable(conn = db, name = "sample_info", value = "sample_info.csv", row.names = FALSE, header = TRUE)
dbWriteTable(conn = db, name = "weights", value = "weights.csv", row.names = FALSE, header = TRUE)
