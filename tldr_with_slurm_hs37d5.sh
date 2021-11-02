#!/bin/bash
#SBATCH --nodes=1
#SBATCH --time=12:00:00
#SBATCH --mem=220G
#SBATCH -c 64
#SBATCH --array=1-2
#SBATCH --partition=12hours
#SBATCH --output=/data/tusers/zhongrenhu/for_SMS/dna_pipline/TLDR_result/logs/tldr-log-%A-%a.out

# Print a little info for debugging
echo "HOSTNAME: " $(hostname)
echo "SLURM_JOB_NODELIST: " $SLURM_JOB_NODELIST
echo "SLURM_JOBID: " $SLURM_JOBID
echo "SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID
echo "SLURM_ARRAY_JOB_ID: " $SLURM_ARRAY_JOB_ID
date
echo ""


######Help Information########
function help_info(){
    echo `basename $0`
    echo -e "\t-e <elts> \tReference elements in .fasta format. The header for each should be formatted as'>Superfamily:Subfamily' e.g. >ALU:AluYb9."
    echo -e "\t-b <bams>\tMultiple .bam files can provided in a comma-delimited list."
    echo -e "\t-r <ref>\tReference genome .fasta, expects a samtools index i.e. samtools faidx."
    echo -e "\t-o <outdir>\toutput path."
    echo -e "\t-h \tShow this information"
}


while getopts ":f:o:h" OPTION; do
    case $OPTION in
        e)  ELTS=$OPTARG;;
        r)  REF=$OPTARG;;
        f)  FILELIST=$OPTARG;;
        o)  OUTPUTDIR=$OPTARG;;
        h)  help_info && exit 1;;
        *)  help_info && exit 1;;
    esac
done


# Parameter initialization
# List of input files, I make this with ls /path/to/folder/*.suffix | sort | pr -1aT -s' ' > filelist.txt
[ -z $FILELIST ] && FILELIST=/data/tusers/zhongrenhu/for_SMS/dna_pipline/TLDR_result/tldr_filelist
[ -z $ELTS ] && ELTS=/data/tusers/zhongrenhu/Software/tldr/ref/teref.ont.human.fa
[ -z $REF ] && REF=/data/tusers/zhongrenhu/for_SMS/reference/hs37d5/hs37d5.fa
[ -z $NON_REF ] && NON_REF=/data/tusers/zhongrenhu/Software/tldr/ref/nonref.collection.hg19.bed.gz


# Get files from FILELIST
echo "Reading filelist..."
str=`sed "${SLURM_ARRAY_TASK_ID}q;d" $FILELIST`
# Read the SLURM_ARRAY_TASK_ID'th line from filelist
# Split line into array
toks=($str)
INPUT=${toks[0]}

# Create temporary directory for work
echo "Creating temporary working dir..."
TMPDIR=$(mktemp -d -t zhongrenhu-tmp-XXXXXXXX)
cd $TMPDIR
echo "Changing wd: " $(pwd)
echo ""

# Copy inputs to temp dir
echo "Copying input files to temp directory..."
cp $INPUT $INPUT.bai .
echo "Done."
echo ""

# Get filenames in temp by cutting of the path
BAM=$(basename $INPUT)
PREFIX=${BAM%.bam}
[ -z $OUTPUTDIR ] && OUTPUTDIR=/data/tusers/zhongrenhu/for_SMS/dna_pipline/TLDR_result/$PREFIX && mkdir $OUTPUTDIR

# exporting environment variable and default parameter
export PATH=/data/tusers/zhongrenhu/Software/anaconda3/bin:$PATH
export PATH=$PATH:/data/tusers/zhongrenhu/Software/minimap2/
export PATH=$PATH:/data/tusers/zhongrenhu/Software/exonerate/bin/
export PATH=$PATH:/data/tusers/zhongrenhu/Software/mafft-7.487-without-extensions/bin
export PATH=$PATH:/data/tusers/zhongrenhu/Software/tldr/tldr/

# Process the data
echo "Running TLDR insertion2 for ${BAM}"
tldr -p 24 -b ${BAM} -e ${ELTS} -r ${REF} --color_consensus -n ${NON_REF}
echo ""

# Copy files to output dir
echo "Copying results to destination..."
ls -lh
rm ./*bam*
cp -r ${PREFIX}* $OUTPUTDIR/
echo "Done."
echo ""


# Clean up
echo "Cleaning up..."
cd $HOME
echo "Deleting temp dir: " $TMPDIR
rm -rd $TMPDIR
#echo "/tmp contents:"
#ls -lh /tmp
echo ""
echo "Script complete."
date

# nohup tldr -b ../HG002.Sequel.10kb.pbmm2.hs37d5.whatshap.haplotag.RTG.10x.trio.bam -e /data/tusers/zhongrenhu/Software/tldr/ref/teref.ont.human.fa
