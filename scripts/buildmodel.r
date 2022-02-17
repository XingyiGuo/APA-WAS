# R version 3.6.0
library(GenABEL)
library(glmnet)
library(sqldf)

qn_pdui <- readRDS("qn_pdui.rds")
peer <- read.table("Breast_complete_data_PEER_factors.txt", header=T, row.names=1, check.names=F, stringsAsFactors=F) 
cov <- read.table("Breast_covs.txt", header=T, row.names=1, check.names=F, stringsAsFactors=F)

cov <- cov[,1:7]; cov <- cov[,-3]; tmp_cov <- apply(as.matrix(cov), 2, as.numeric); rownames(tmp_cov) <- rownames(cov); cov <- tmp_cov
rname <- c(); for (i in 1:nrow(peer)) {rname <- c(rname, paste('GTEX', strsplit(rownames(peer)[i],'-')[[1]][2],sep='-'))}; rownames(peer) <- rname
tmp_peer <- apply(as.matrix(peer), 2, as.numeric); rownames(tmp_peer) <- rownames(peer); peer <- tmp_peer

rn_pdui <- as.numeric(rntransform(as.numeric(qn_pdui[,1])))
tmp_d <- data.frame(cbind(rn_pdui, cov, peer[,1:15]), stringsAsFactors=F)
res <- residuals(lm(rn_pdui ~ ., tmp_d, na.action="na.exclude"))
rn_res <- as.numeric(rntransform(res))
rn_res <- data.frame(cbind(rownames(qn_pdui), rn_res), stringsAsFactors=F)

pdui_val <- scale(as.numeric(rn_res$rn_res), center=T, scale=T); pdui_val[is.na(pdui_val)] <- 0; rownames(pdui_val) <- as.character(rn_res$V1)
geno <- readRDS("dosage.rds"); geno <- apply(geno, 2, function(x) ifelse(is.na(x),mean(x,na.rm=T),x))
snp_info <- readRDS("snpinfo.rds"); rownames(snp_info) <- as.character(snp_info$SNP)

ext_head <- c("gene","genename","pred.perf.R2","n.snps.in.model","pred.perf.pval","pred.perf.qval")
wt_head <- c("gene", "rsid", "eff_allele", "ref_allele", "weight")
cov_head <- c('GENE', 'RSID2','RSID1','VALUE')

write(ext_head,file="extra.csv",ncolumns=6,sep=",")
write(wt_head,file="weights.csv",ncolumns=5,sep=",")
write(cov_head,file="cov.txt",ncolumns=4,sep="\t")

set.seed(05162020)
fit <- cv.glmnet(geno,pdui_val,nfolds=10,alpha=0.5,keep=T,parallel=F); fit.df <- data.frame(fit$cvm,fit$lambda,1:length(fit$cvm))
best.lam <- fit.df[which.min(fit.df[,1]),]; cvm.best = best.lam[,1]; lambda.best = best.lam[,2]; nrow.best = best.lam[,3]
ret <- as.data.frame(fit$glmnet.fit$beta[,nrow.best]) # get betas from best lambda
ret[ret == 0.0] <- NA
bestbetas = as.vector(ret[which(!is.na(ret)),]) # vector of non-zero betas
names(bestbetas) = rownames(ret)[which(!is.na(ret))]
pred.mat <- fit$fit.preval[,nrow.best] # pull out predictions at best lambda
res <- cor.test(pdui_val, pred.mat); rval <- res$estimate[[1]]; pval <- res$p.value
sum_res <- c("NM_002035.4|KDSR|chr18|-", "NM_002035.4|KDSR|chr18|", rval**2, length(bestbetas), pval, "0")
write(sum_res,file="extra.csv",ncolumns=6,append=T,sep=",")
bestbetalist <- names(bestbetas)
bestbetainfo <- snp_info[bestbetalist,]
betatable<-as.matrix(cbind(bestbetainfo,bestbetas))
betafile <- cbind("NM_002035.4|KDSR|chr18|-",paste("rs",rownames(betatable),sep=""),betatable[,3],betatable[,4],as.numeric(betatable[,6]))
write(t(betafile),file="weights.csv",ncolumns=5,append=T,sep=",")
dsg <- geno[,bestbetalist]; cov = cov(as.matrix(dsg)); cov[upper.tri(cov)] <- NA; cov = cbind(expand.grid(dimnames(cov)), value = as.vector(cov))
colnames(cov) <- c('RSID2','RSID1','VALUE'); cov = cov[!is.na(cov$VALUE),]; cov$GENE <- "NM_002035.4|KDSR|chr18|-"; cov = cov[,c('GENE','RSID1','RSID2','VALUE')]
cov$RSID1 <- paste("rs", cov$RSID1, sep=""); cov$RSID2 <- paste("rs", cov$RSID2, sep="")
write(t(cov),file="cov.txt",ncolumns=4,append=T,sep="\t")

if ("cov.txt.gz" %in% dir()) {system("rm cov.txt.gz")}
system("gzip cov.txt")
write(c("chr","cv.seed"),file="construction.csv",ncolumns=2,sep=",")
write(c("NM_002035.4|KDSR|chr18|-","05162020"),file="construction.csv",ncolumns=2,sep=",",append=T)
write("n.samples\n114\n",file="sample_info.csv",ncolumns=1)

if ("breast.db" %in% dir()) {system("rm breast.db")}
db <- dbConnect(SQLite(), dbname="breast.db")
dbWriteTable(conn = db, name = "construction", value = "construction.csv", row.names = FALSE, header = TRUE)
dbWriteTable(conn = db, name = "extra", value = "extra.csv", row.names = FALSE, header = TRUE)
dbWriteTable(conn = db, name = "sample_info", value = "sample_info.csv", row.names = FALSE, header = TRUE)
dbWriteTable(conn = db, name = "weights", value = "weights.csv", row.names = FALSE, header = TRUE)

