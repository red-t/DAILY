#!/bin/bash

samtools sort -n -@ 32 /data/tusers.ds/zhongrenhu/for_SMS/dna_pipline/TEMP2_result/HG002.hs37d5.300x_chr20/HG002.hs37d5.300x_chr20.bam -o HG002.hs37d5.300x_chr20_sorted.bam
samtools fastq -@ 32 HG002.hs37d5.300x_chr20_sorted.bam -1 HG002.hs37d5.300x_chr20_r1.fastq.gz -2 HG002.hs37d5.300x_chr20_r2.fastq.gz -0 /dev/null -s /dev/null -n

GENOME=hs37d5
PATH_PROG=`dirname ${0}` && PATH_ANNO="/data/tusers/yutianx/tongji2/GitHuB/piSet/annotation/${GENOME}"
[ -z ${BWA_INDEX} ] && BWA_INDEX=${PATH_ANNO}/BWAIndex/genome
[ -z ${GENOME_FA} ] && GENOME_FA=${PATH_ANNO}/${GENOME}.fa
[ -z ${TE_FA} ] && TE_FA=${PATH_ANNO}/ALUL1SVA.fa
[ -z ${TE_BED6} ] && TE_BED6=${PATH_ANNO}/ALUL1SVA.bed
[ -z ${OUT_PATH} ] && OUT_PATH=./TEMP2_result
[ -z ${CPU} ] && CPU=32

PREFIX='HG002.hs37d5.300x_chr20'
bwa mem -t ${CPU} ${BWA_INDEX} HG002.hs37d5.300x_chr20_r1.fastq.gz HG002.hs37d5.300x_chr20_r2.fastq.gz > ${PREFIX}.sam 2>${PREFIX}.bwamem.log || \
		{ echo0 0 "Error: bwa mem failed, please check ./${PREFIX}.bwamem.log. Exiting..." && exit 1; }
echo0 2 "transform sam to sorted bam and index it"
samtools view -bhS -@ ${CPU} ${PREFIX}.sam > ${PREFIX}.bam
samtools sort -@ ${CPU} -o ${PREFIX}.sorted.bam ${PREFIX}.bam
rm ${PREFIX}.sam ${PREFIX}.bam
samtools index -@ ${CPU} ${PREFIX}.sorted.bam
BAM=${PREFIX}.sorted.bam
samtools flagstat -@ ${CPU} ${PREFIX}.sorted.bam

samtools view -H ${BAM} > ${PREFIX}.tmp.header
samtools view -@ ${CPU} -h -f 0X2 -F 0X800 ${BAM} | awk -v uniq_ratio=0.8 'BEGIN{FS=OFS="\t"} {for(i=12;i<=100;i++){if($i~/AS:i:/){as=substr($i,6)};if($i~/XS:i:/){xs=substr($i,6);break}};if(xs/1<as*uniq_ratio){print $0}}' - | grep "SA:Z:" | awk '(and($2,16)==0 && $6~/^[0-9]+S/) || (and($2,16)==16 && $6~/S$/)' - | cat ${PREFIX}.tmp.header - | samtools view -@ ${CPU} -bhS - | samtools sort -@ ${CPU} -o ${PREFIX}.pair.uniq.split.bam -
#samtools view -@ ${CPU} -h -f 0X2 -F 0X800 ${BAM} > ${PREFIX}.tmp.1.sam
#awk -v uniq_ratio=0.8 'BEGIN{FS=OFS="\t"} {for(i=12;i<=100;i++){if($i~/AS:i:/){as=substr($i,6)};if($i~/XS:i:/){xs=substr($i,6);break}};if(xs/1<as*uniq_ratio){print $0}}' ${PREFIX}.tmp.1.sam > ${PREFIX}.tmp.2.sam