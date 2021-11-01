#!/bin/bash

# function
function checkTools(){
    if [ `which` $1 ];then
        echo0 3 "$1: `which $1`"
    else
        echo0 0 "$1 not found, please install it or add it into your PATH!"
        exit 1
    fi
}

function delet(){
    [ -f $1 ] && rm -rf $1
}

# help information
function help_info(){
    echo `basename $0`" -l '<left fastqs>' -r '<right fastqs>' -G <genome> -c <cpu>"
    echo `basename $0`" -i <bams> -G <genome> -c <cpu>"
    echo -e "\t-l <left fastq(s)>\tleft fastq file(s)."
    echo -e "\t-r <right fastq(s)>\tright fastq file(s)."
    echo -e "\t-i <bams>\tmapping result by BWA MEM in BAM format."
    echo -e "\t-I <BWA Index>\tBWA index for genome. Or using -G option."
    echo -e "\t-g <genome.fasta>\treference genome sequence in FASTA format."
    echo -e "\t-G <genome>\tname of reference genome. Or using -I option to specify BWA index."
    echo -e "\t-R <te.fasta>\transposon concensus sequence in FASTQ format."
    echo -e "\t-t <bed>\tendougenous transposon copies in BED format. Produced by repeatmasker."
    echo -e "\t-o <folder>\tnoutput directory."
    echo -e "\t-c <cpu>\tnummber of CPU to use."
    echo -e "\t-h \tShow this information"
}

if [ $# -lt 1 ];then
    help_info && exit 1
fi

# default parmeters
#CPU=8
#GENOME=dm3


while getopts ":l:r:i:I:g:G:R:t:o:c:h" OPTION; do
    case $OPTION in
        l)  LEFT=($OPTARG);;
        r)  RIGHT=($OPTARG);;
        i)  MAP_BAM=($OPTARG);;
        I)  BWA_INDEX=$OPTARG;;
        g)  GENOME_FA=$OPTARG;;
        G)  GENOME=$OPTARG;;
        R)  TE_FA=$OPTARG;;
        t)  TE_BED6=$OPTARG;;
        o)  OUT_PATH=$OPTARG;;
        c)  CPU=$OPTARG;;
        h)  help_info && exit 1;;
        *)  help_info && exit 1;;
    esac
done

PATH_PROG=`dirname ${0}` && PATH_ANNO="/data/tusers/yutianx/tongji2/GitHuB/piSet/annotation/${GENOME}"
[ -z ${BWA_INDEX} ] && BWA_INDEX=${PATH_ANNO}/BWAIndex/genome
[ -z ${GENOME_FA} ] && GENOME_FA=${PATH_ANNO}/${GENOME}.fa
[ -z ${TE_FA} ] && TE_FA=${PATH_ANNO}/ALUL1SVA.fa
[ -z ${TE_BED6} ] && TE_BED6=${PATH_ANNO}/ALUL1SVA.bed
[ -z ${OUT_PATH} ] && OUT_PATH=./TEMP2_result
[ -z ${CPU} ] && CPU=8


# check tools
echo0 1 "check tools......"
checkTools TEMP2
checkTools bwa
checkTools samtools
echo0 1 "OK"



# check parameters
echo0 1 "check parameters"
[ ! -d ${PATH_ANNO} ] && echo0 4 "genome ${GENOME} not installed" && exit 1

if [ -n ${LEFT} ];then
    for i in ${LEFT[*]}
    do
        [ ! -f "${i}" ] && echo0 0 "Error: No fastq file in ${i}" && exit 1
    done

    for i in ${RIGHT[*]}
    do
        [ ! -f "${i}" ] && echo0 0 "Error: No fastq file in ${i}" && exit 1
    done

    if [ -z ${PREFIX} ];then
        NUM=0
        for i in ${LEFT[*]}
        do
            PREFIX[$NUM]=`basename ${i%[._][r12]*f*q*}`
            NUM=$(($NUM + 1))
        done
    fi
fi

if [ -n ${MAP_BAM} ];then
    if [ -z ${PREFIX} ];then
        NUM=0
        for i in ${MAP_BAM[*]}
        do
            PREFIX[$NUM]=`basename ${i%.bam}`
            NUM=$(($NUM + 1))
        done
    fi
fi

echo0 1 "OK"


#############
## process ##

[ ! -d ${OUT_PATH} ] && mkdir ${OUT_PATH}

if [ -n ${MAP_BAM} ];then
    SAMPLE_INDEX=0
    for BAM in ${MAP_BAM[*]}
    do
        echo0 1 "processing for ${PREFIX[${SAMPLE_INDEX}]}" && date
        OUT_DIR=${OUT_PATH}/${PREFIX[${SAMPLE_INDEX}]} && mkdir ${OUT_DIR}
        cp `dirname ${BAM}`/${PREFIX[${SAMPLE_INDEX}]}*bam* ${OUT_DIR} && cp `dirname ${BAM}`/${PREFIX[${SAMPLE_INDEX}]}*f*q* ${OUT_DIR}
        BAM=${OUT_DIR}/${PREFIX[${SAMPLE_INDEX}]}.bam && TEMP_LEFT=${OUT_DIR}/${PREFIX[${SAMPLE_INDEX}]}[_.]1.f*q* && TEMP_RIGHT=${OUT_DIR}/${PREFIX[${SAMPLE_INDEX}]}[_.]2.f*q*
        TEMP2 insertion -l ${TEMP_LEFT} -r ${TEMP_RIGHT} -i ${BAM} -g ${GENOME_FA} -R ${TE_FA} -t ${TE_BED6} -o ${OUT_DIR} -p ${PREFIX[${SAMPLE_INDEX}]} -c ${CPU} -M 2 -m 10 -U 0.95 -S -A #-d
        SAMPLE_INDEX=$((${SAMPLE_INDEX} + 1))
        # mv ${OUT_DIR} /data/tusers/zhongrenhu/for_SMS/dna_pipline/TEMP2_result && date
    done
elif [ -z ${MAP_BAM} ];then
    SAMPLE_INDEX=0
    for TEMP_LEFT in ${LEFT[*]}
    do
        echo0 1 "processing for ${PREFIX[${SAMPLE_INDEX}]}" && date
        OUT_DIR=${OUT_PATH}/${PREFIX[${SAMPLE_INDEX}]}
        cp `dirname ${TEMP_LEFT}`/${PREFIX[${SAMPLE_INDEX}]}*f*q* ${OUT_DIR}
        TEMP_LEFT=${OUT_DIR}/${PREFIX[${SAMPLE_INDEX}]}[_.]1.f*q* && TEMP_RIGHT=${OUT_DIR}/${PREFIX[${SAMPLE_INDEX}]}[_.]2.f*q*
        TEMP2 insertion -l ${TEMP_LEFT} -r ${TEMP_RIGHT} -I ${BWA_INDEX} -g ${GENOME_FA} -R ${TE_FA} -t ${TE_BED6} -o ${OUT_DIR} -p ${PREFIX[${SAMPLE_INDEX}]} -c ${CPU} -M 2 -m 10 -U 0.95 -S -A #-d
        SAMPLE_INDEX=$((${SAMPLE_INDEX} + 1))
        # mv ${OUT_DIR} /data/tusers/zhongrenhu/for_SMS/dna_pipline/TEMP2_result && date
        date
    done
fi


###finished
echo0 4 "------finished🍺🍺🍺------"


#nohup /data/tusers/zhongrenhu/for_SMS/bin/TEMP2_local.sh -l /data/tusers.ds/zhongrenhu/for_SMS/dna/GIAB_HG002/HG002.hs37d5.300x_chr20_r1.fastq.gz -r /data/tusers.ds/zhongrenhu/for_SMS/dna/GIAB_HG002/HG002.hs37d5.300x_chr20_r2.fastq.gz -i /data/tusers.ds/zhongrenhu/for_SMS/dna/GIAB_HG002/HG002.hs37d5.300x_chr20.bam -G hs37d5 -c 32 -o ./ &

# nohup ../../bin/gp_rnaseq.sh -l "./SRR6189*" -g dm3 -o ./result -L miss --transposon-analysis &
# nohup ../../bin/ssnv_dnaseq.sh -l "*_1.fastq" -r "*_2.fastq" -G hs37d5 -c 20 &

# nohup tldr -b ../HG002.Sequel.10kb.pbmm2.hs37d5.whatshap.haplotag.RTG.10x.trio.bam -e /data/tusers/zhongrenhu/Software/tldr/ref/teref.ont.human.fa