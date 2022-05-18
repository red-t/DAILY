#!/bin/bash

######Help Information########
function help_info(){
    echo `basename $0`
    echo -e "\t-b <bam(s)>\tBAM files that used for downsampling."
    echo -e "\t-c <int>\tnumber of CPUs to use."
    echo -e "\t-d <depth(s)>\tsequencing depth of each BAM files."
    echo -e "\t-t <depth(s)>\ttarget depth to downsample."
    echo -e "\t-o <outdir>\tdirectory to output."
    echo -e "\t-h \tShow this information"
}

if [ $# -lt 1 ];then
    help_info && exit 1
fi

while getopts ":b:c:d:t:o:h" OPTION; do
    case $OPTION in
        b)  BAMs=($OPTARG);;
        c)  CPU=$OPTARG;;
        d)  R_DEPTH=($OPTARG);;
        t)  T_DEPTH=($OPTARG);;
        o)  OUT_PATH=$OPTARG;;
        h)  help_info && exit 1;;
        *)  help_info && exit 1;;
    esac
done


#############
## process ##
[ -z ${CPU} ] && CPU=20
[ -z ${OUT_PATH} ] && OUT_PATH="./"
[ ! -d ${OUT_PATH} ] && mkdir ${OUT_PATH}
SAMPLE_INDEX=0

for bam in ${BAMs[*]}
do 
    echo0 1 "downsampling for ${bam}" && date
    id=`basename ${bam}` && prefix=${id%.bam}

    ## downsampling
    DEPTH_INDEX=0
    for tdepth in ${T_DEPTH[*]}
    do
        fraction=$(awk -v rdepth=${R_DEPTH[${SAMPLE_INDEX}]} -v tdepth=${tdepth} 'BEGIN{printf "%u", (tdepth/rdepth)*100}')
        samtools view -b -1 -@ ${CPU} -s ${SAMPLE_INDEX}.${fraction} -o ${OUT_PATH}/${prefix}.${tdepth}X.bam ${bam}
        DEPTH_INDEX=$((${DEPTH_INDEX} + 1))
    done

    SAMPLE_INDEX=$((${SAMPLE_INDEX} + 1))
    date
done
