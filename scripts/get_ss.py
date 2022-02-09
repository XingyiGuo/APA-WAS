import sys, math

ss = {}
for line in open("Prostate/SumStat/meta_v3_onco_euro_overall_ChrAll_1_release.txt"):
    l = line.strip().split("\t")
    snpname, rsid, snp, c, p, eff_A, ot_A, frq1, frqse, minfrq, maxfrq, beta, se, pval, direct, r2 = l
    eff_A = eff_A.upper(); ot_A = ot_A.upper()
    tag = ":".join([c,p,eff_A,ot_A]); ss[tag] = [snpname, eff_A, ot_A, beta, se, pval]
    tag1 = ":".join([c,p,ot_A,eff_A]); ss[tag1] = [snpname, eff_A, ot_A, beta, se, pval]

print("snp,effA,otA,beta,se,pval")
for line in open("../weights.uniq.snp"):
    snp = line.strip().lstrip('rs')
    res = ss[snp]; res[0] = line.strip(); print(",".join(res))
