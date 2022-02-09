#/bin/sh

### Prostate as example
### sumstat was download from http://practical.icr.ac.uk/blog/wp-content/uploads/uploadedfiles/oncoarray/MetaSummaryData/meta_v3_onco_euro_overall_ChrAll_1_release.zip
wget meta_v3_onco_euro_overall_ChrAll_1_release.txt
unzip meta_v3_onco_euro_overall_ChrAll_1_release.zip
less meta_v3_onco_euro_overall_ChrAll_1_release.txt | perl -lane '$a1 = uc $F[5]; $a2 = uc $F[6]; print "$F[3]:$F[4]:$a1:$a2"; print "$F[3]:$F[4]:$a2:$a1"' > prst.snp
less prst.snp |grep -v $':A:T$'|grep -v $':T:A$'| grep -v $':C:G$'|grep -v $':G:C$' > prst.snp.noAmbiguous

less Prostate_PDUI.txt | head -n 1 | sed 's/\t/\n/g' | less | tail -n +5 | sed 's/.*GTEX/GTEX/' | cut -f 1-2 -d '-' > pdui.id

less gtex.eur.male.qc.fam |cut -f 1 > geno.id # QC-ed genotype data of European males in GTEx. "--geno 0.05 --maf 0.05 --hwe 1E-04" used on data of all samples. Then, extract males, apply "--maf 0.05". 

cat geno.id pdui.id |sort|uniq -d > overlap.id

plink2 --bfile gtex.eur.male.qc --extract  prst.snp.noAmbiguous --keep-fam overlap.id --rm-dup force-first --make-bed --out prst

## Split by chromosome, generated bfiles and .traw failes.
for i in {1..22}
do
  mkdir Chr${i}
  plink2 --bfile prst --chr ${i} --make-bed --out Chr${i}/chr${i}
  plink2 --bfile Chr${i}/chr${i} --export A-transpose --out Chr${i}/chr${i}
done
