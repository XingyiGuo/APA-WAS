# Python version 3.6.3
import os, subprocess

if "assoc.csv\n" in os.popen("ls"): os.popen("rm assoc.csv")   
subprocess.call("SPrediXcan.py --model_db_path breast.db --covariance cov.txt.gz --gwas_file gwas.csv --snp_column snp --effect_allele_column effA --non_effect_allele_column otA --beta_column beta --se_column se --pvalue_column pval --separator , --output_file assoc.csv".split())

cojo = ['plink2 --bfile /scratch/sbcs/yhyang/Refs/Chr18/chr18 --extract rsid --make-bed --out ref1KG', 'echo "rs17743054" > cond.snp', 'gcta64 --bfile ref1KG --cojo-cond cond.snp --cojo-file predictingSNP_indexSNP.ss --cojo-collinear 0.99 --out cojores']
for c in cojo:
    subprocess.call(c.split())

d = {}
for line in open("snp.rsid"):
    snp, rsid = line.strip().split()
    d[rsid] = snp

output = open("cojores.gwas","w")
output.write("snp,effA,otA,beta,se,pval\n")
for line in open("cojores.cma.cojo"):
    if line.startswith("Chr"):continue
    l = line.strip().split()
    snp = d[l[1]]; effA = l[3]
    a1, a2 = snp.split(":")[2:]; otA = a2 if effA == a1 else a1
    res = [snp, effA, otA] + l[-3:]
    output.write(",".join(res)+"\n")
output.close()

if "cojores.assoc.csv\n" in os.popen("ls"): os.popen("rm cojores.assoc.csv")   

adjcojo = 'SPrediXcan.py --model_db_path breast.db --covariance cov.txt.gz --gwas_file cojores.gwas --snp_column snp --effect_allele_column effA --non_effect_allele_column otA --beta_column beta --se_column se --pvalue_column pval --separator , --output_file cojores.assoc.csv' 
subprocess.call(adjcojo.split())
