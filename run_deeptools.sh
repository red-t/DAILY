#!/bin/bash

# help information
function help_info(){
    echo -e `basename $0` " [options] -b 'bam(s)'"
    echo -e "\t-b <BAM(s)>\tBAM file(s) to process"
    echo -e "\t-T <toolname>\twhich deepTools tool to use. possible choices: bigwig, bedgraph"
    echo -e "\t-t <threads>\tnumber of threads, default:25"
}

if [ $# -lt 1 ];then
    help_info && exit 1
fi

while getopts ":b:h" OPTION; do
    case $OPTION in
        b)  BAMs=($OPTARG);;
        h)  help_info && exit 1;;
        *)  help_info && exit 1;;
    esac
done

for BAM in ${BAMs[*]}
do
    PREFIX=`basename ${BAM%.bam}`
    bamCoverage -b $BAM -o ./${PREFIX}.bw -of bigwig --ignoreDuplicates --normalizeUsing BPM --numberOfProcessors 25
done