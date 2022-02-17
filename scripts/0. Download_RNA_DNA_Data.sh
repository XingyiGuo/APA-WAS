#!/bin/sh

### Download GTEx RNA-Seq data from AnVIL

gen3clientDir=
downDir=
manifestFile=

# Breast (n=114) / Ovary (n=127) / Prostate (n=170) / Lung (n=170) / Colon Transverse (n=285) / Pancreas (n=231)

${gen3clientDir} download-multiple --profile=AnVIL --manifest=${manifestFile} --download-path=${downDir} --protocol=s3

### Download GTEx DNA-Seq data from dbGap # phs000424.GTEx.v8.p2.c1
### need permission to download from dbGap
