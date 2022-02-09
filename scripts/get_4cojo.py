import sys

def get_frq(chrom):
    frqd = {}
    for line in open("chr%s.afreq" % chrom):  # chr%s.afreq are allele frequency files from 1000 Genome European data.
        if line.startswith("#CHROM"): continue
        c, p, rsid, refA, altA, altAF, obs_ct = line.strip().split("\t")
        tag1 = ":".join([c,p,altA,refA]); frqd[tag1] = [rsid, altAF]
        tag2 = ":".join([c,p,refA,altA]); frqd[tag2] = [rsid, str(1-float(altAF))]
    return frqd

def main():    
    d = {}
    for line in open('weights.csv'): # weights file from model building
        if line.startswith("gene,rsid"): continue
        l = line.strip().split(",")
        if l[0] not in d:
            d[l[0]] = []
        d[l[0]].append(l[1])

    ss = {}
    for line in open("meta_v3_onco_euro_overall_ChrAll_1_release.txt"):
        if line.startswith("MarkerName"): continue
        l = [i.strip() for i in line.strip().split("\t")]
        chrom, pos, eff_A, ot_A = l[3:7]; eff_A = eff_A.upper(); ot_A = ot_A.upper()
        tag1=':'.join([chrom, pos, eff_A, ot_A]); tag2 =':'.join([chrom, pos, ot_A, eff_A])
        ss[tag1] = [eff_A, ot_A] + l[-5:-2]; ss[tag2] = [ot_A, eff_A, str(-float(l[-5]))] + l[-4:-2]

    for line in open("sig.index"): # significant associations with proximal index SNPs
        l = line.strip().split(",")
        outf = open("%s.4cojo" % l[0].split("|")[0],"w"); outf.write("SNP A1 A2 freq b se p N \n")
        snps = d[l[0]]; snps.append(l[-1].split("_")[-1])
        c = l[-2].lstrip("chr").split(":")[0]; frq = get_frq(c)
        for s in snps:
            sumstat = ss[s] if s in ss else ["NA"]*5
            rs, af = frq[s] if s in frq else [s,"NA"]
            res = [rs] + sumstat[:2] + [af] + sumstat[2:] + ['140306']
            outf.write(' '.join(res)+'\n')
        outf.close()

main()
