#!/bin/bash

# help information
function help_info(){
    echo `basename $0`" [options] -R reference.fa -Q query.fastq"
    echo -e "\t-R <ref.fa>\treference FASTA"
    echo -e "\t-Q <query.fastq/fasta>\tquery FASTQ/FASTQ"
    echo -e "\t-a output in the SAM format (PAF by default)"
    echo -e "\t-t <threads>\tnumber of threads"
    echo -e "Preset:
    -x STR       preset (always applied before other options; see minimap2.1 for details) []
                 - map-pb/map-ont - PacBio CLR/Nanopore vs reference mapping
                 - map-hifi - PacBio HiFi reads vs reference mapping
                 - ava-pb/ava-ont - PacBio/Nanopore read overlap
                 - asm5/asm10/asm20 - asm-to-ref mapping, for ~0.1/1/5% sequence divergence
                 - splice/splice:hq - long-read/Pacbio-CCS spliced alignment
                 - sr - genomic short-read mapping"
    echo -e "\t-h \tShow this information"
}

if [ $# -lt 1 ];then
    help_info && exit 1
fi

while getopts ":R:Q:at:x:h" OPTION; do
    case $OPTION in
        R)  REF_FASTA=($OPTARG);;
        Q)  QUERY=($OPTARG);;
        a)  SAM_FORMAT=1;;
        t)  THREADS=$OPTARG;;
        x)  PRESET=$OPTARG;;
        h)  help_info && exit 1;;
        *)  help_info && exit 1;;
    esac
done

PREFIX=`basename ${QUERY%.f*}`
if [ -n $SAM_FORMAT ];then
    minimap2 -a -x $PRESET --MD -t $THREADS $REF_FASTA $QUERY > $PREFIX.sam
    samtools view -bhS -@ $THREADS $PREFIX.sam | samtools sort -@ $THREADS -o $PREFIX.bam -
    else
    minimap2 -x $PRESET --MD -t $THREADS $REF_FASTA $QUERY > $PREFIX.paf
fi

# minimap2_local.sh -R /data/tusers.ds/zhongrenhu/for_SMS/reference/phaCin/phaCin.fa -Q /data/tusers.ds/zhongrenhu/for_SMS/dna/koala/Koala2_20191231_Adelaide.fastq -t 32 -x map-ont
# nohup minimap2_local.sh -R /data/tusers/zhongrenhu/for_SMS/reference/dm3/dm3.fa -Q line_28_pacbio.fasta -t 20 -x map-pb -a &