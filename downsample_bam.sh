#!/bin/bash

######Help Information########
function help_info(){
    echo `basename $0`
    echo -e "\t-b <bam(s)>\tBAM files that used for downsampling."
    echo -e "\t-c <int>\tnumber of CPUs to use."
    echo -e "\t-d <depth(s)>\tsequencing depth of each BAM files."
    echo -e "\t-f <format>\tbam/fastq/fasta."
    echo -e "\t-t <depth(s)>\ttarget depth to downsample."
    echo -e "\t-o <outdir>\tdirectory to output."
    echo -e "\t-h \tShow this information"
}

if [ $# -lt 1 ];then
    help_info && exit 1
fi

while getopts ":b:c:d:f:t:o:h" OPTION; do
    case $OPTION in
        b)  BAMs=($OPTARG);;
        c)  CPU=$OPTARG;;
        d)  R_DEPTH=($OPTARG);;
        f)  FORMAT=$OPTARG;;
        t)  T_DEPTH=($OPTARG);;
        o)  OUT_PATH=$OPTARG;;
        h)  help_info && exit 1;;
        *)  help_info && exit 1;;
    esac
done


#############
## process ##
[ -z ${CPU} ] && CPU=20
[ -z ${FORMAT} ] && FORMAT="bam"
[ -z ${OUT_PATH} ] && OUT_PATH="./"
[ ! -d ${OUT_PATH} ] && mkdir ${OUT_PATH}

SAMPLE_INDEX=0
for bam in ${BAMs[*]}
do 
    id=`basename ${bam}` && prefix=${id%.bam}

    # downsampling
    DEPTH_INDEX=0
    for tdepth in ${T_DEPTH[*]}
    do
        echo "downsampling ${bam} to ${tdepth}X start" && date
        fraction=`awk -v rdepth=${R_DEPTH[${SAMPLE_INDEX}]} -v tdepth=${tdepth} 'BEGIN{printf "%u", (tdepth/rdepth)*100}'`
        samtools view -h -b -@ ${CPU} -s ${SAMPLE_INDEX}.${fraction} -o ${OUT_PATH}/${prefix}.${tdepth}X.bam ${bam}
        samtools index ${OUT_PATH}/${prefix}.${tdepth}X.bam
        
        # transfer BAM format into FASTQ/FASTA
        if [ $FORMAT == 'fastq' ];then
            samtools sort -n -@ ${CPU} -o tmp.bam ${OUT_PATH}/${prefix}.${tdepth}X.bam
            samtools fastq -@ ${CPU} -1 ${OUT_PATH}/${prefix}.${tdepth}X_1.fq -2 ${OUT_PATH}/${prefix}.${tdepth}X_2.fq -0 /dev/null -s /dev/null -n tmp.bam
            rm tmp.bam
        elif [ $FORMAT == 'fasta' ];then
            samtools fasta -@ ${CPU} ${OUT_PATH}/${prefix}.${tdepth}X.bam > ${OUT_PATH}/${prefix}.${tdepth}X.fa
            samtools faidx ${OUT_PATH}/${prefix}.${tdepth}X.fa
        fi

        DEPTH_INDEX=$((${DEPTH_INDEX} + 1))
        echo "downsampling ${bam} to ${tdepth}X end" && date
    done

    SAMPLE_INDEX=$((${SAMPLE_INDEX} + 1))
done


# downsample_bam.sh -b line_28_pacbio.bam -d 50 -f fasta -t "10 20 30 40" -o ./