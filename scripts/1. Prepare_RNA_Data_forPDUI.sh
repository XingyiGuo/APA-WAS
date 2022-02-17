#!/bin/sh

### Generate BedGraph from BAM files ###

gtexDir=
chrsizeFile=

tissue=

bamDir=${gtexDir}/${tissue}

while read sample
do
    genomeCoverageBed -bg -ibam ${bamDir}/${sample}.Aligned.sortedByCoord.out.patched.md.bam -g ${chrsizeFile} -split > ${bamDir}/BedGraph/${sample}.bedgraph
done < ${bamDir}/${tissue}.samples.txt

### Calculate read depth ###

gtexDir=
depthDir=

tissue=

bamDir=${gtexDir}/${tissue}
[ -d "${depthDir}/${tissue}" ] && mkdir ${depthDir}/${tissue}

while read sample
do
    bamLoc=${bamDir}/${sample}.Aligned.sortedByCoord.out.patched.md.bam
    wigLoc=${bamDir}/BedGraph/${sample}.bedgraph
    echo "wigLoc" >> ${depthDir}/${tissue}/${sample}.depth
    samtools view -c -F 260 ${bamLoc} >> ${depthDir}/${tissue}/${sample}.depth
done < ${bamDir}/${tissue}.samples.txt

### Generate BAM location with depth ###

tissue=

for f in $(ls ${depthDir}/${tissue}/*.depth)
do
    cat ${f} | awk '{printf $0"\t"}' | awk '{print $1"\t"$2}' >> ${depthDir}/${tissue}.mapping_bam_location_with_depth.txt
done
