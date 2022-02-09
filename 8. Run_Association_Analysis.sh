#!/bin/sh

less weights.rsid.csv | cut -f 2 -d ','| tail -n +2 | sort | uniq >  weights.uniq.snp

python scripts/get_ss.py > gwas.csv

MetaXcan.py \
    --model_db_path prst.db \
    --covariance cov.rsid.txt.gz \
    --gwas_folder ./ \
    --gwas_file_pattern "gwas.csv" \
    --separator , \
    --snp_column snp \
    --effect_allele_column effA \
    --non_effect_allele_column otA \
    --beta_column beta \
    --se_column se \
    --output_file res.csv
