#!/bin/sh

### Generate DaPars2 configure file ###

Annotated_3UTR_File=
outDir=
Coverage_threshold=10
Num_Threads=16
depthDir=

for tissue in Breast Colon_Transverse Lung Ovary Pancreas Prostate
do
	mkdir ${outDir}/${tissue}
    sequencing_depth_file=${depthDir}/${tissue}.mapping_bam_location_with_depth.txt

    wigfilelist=$(awk '{print $1}' ${sequencing_depth_file} | tr "\n" "," | sed 's/.$//' )

    echo "Annotated_3UTR=${Annotated_3UTR_File}" > ${outDir}/${tissue}.configure

    echo "Aligned_Wig_files=${wigfilelist}" >> ${outDir}/${tissue}.configure

    echo "Output_directory=${tissue}/" >> ${outDir}/${tissue}.configure
    echo "Output_result_file=${tissue}" >> ${outDir}/${tissue}.configure
    echo "Coverage_threshold=${Coverage_threshold}" >> ${outDir}/${tissue}.configure
    echo "Num_Threads=${Num_Threads}" >> ${outDir}/${tissue}.configure
    echo "sequencing_depth_file=${sequencing_depth_file}" >> ${outDir}/${tissue}.configure
done

### Run DaPars2 by tissue type ###

pduiDir=

tissue=

ml load GCC/6.4.0-2.28 OpenMPI/2.1.1 Python/2.7.14 numpy/1.13.1-Python-2.7.14 scipy/0.19.1-Python-2.7.14 R/3.4.3

dapar2py=DaPars2_Multi_Sample_Multi_Chr.py

cd ${pduiDir}/${tissue}

python2 ${dapar2py} ${pduiDir}/${tissue}.configure

### Merge DaPars results ###

outDir=

for tissue in Breast Colon_Transverse Lung Ovary Pancreas Prostate
do
    echo $tissue
    head -n 1 ${outDir}/${tissue}/${tissue}_chr1/${tissue}_result_temp.chr1.txt > ${outDir}/${tissue}_PDUI.txt
	for n in {1..22}
	do
	    tail -n +2 ${outDir}/${tissue}/${tissue}_chr${n}/${tissue}_result_temp.chr${n}.txt >> ${outDir}/${tissue}_PDUI.txt
    done
done
