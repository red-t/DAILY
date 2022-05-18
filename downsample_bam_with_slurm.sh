#!/bin/bash
#SBATCH --nodes=1
#SBATCH --time=4:00:00
#SBATCH --mem=90G
#SBATCH -c 64
#SBATCH --array=1-3
#SBATCH --partition=4hours
#SBATCH --output=/data/tusers/zhongrenhu/for_TE_and_sSNV/rawdata/dna/slurm-logs/downsample-log-%A-%a.out

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
    echo -e "\t-f <filelist>\tfilelist of the BAM files to downsample, with sequencing depth."
    echo -e "\t-o <outdir>\toutput path."
    echo -e "\t-s <seed(s)>\tseed used for downsampling."
    echo -e "\t-h \tShow this information"
}

if [ $# -lt 1 ];then
    help_info && exit 1
fi

while getopts ":f:o:s:h" OPTION; do
    case $OPTION in
        f)  FILELIST=$OPTARG;;
        o)  OUTPUTDIR=$OPTARG;;
        s)  SEEDS=($OPTARG);;
        h)  help_info && exit 1;;
        *)  help_info && exit 1;;
    esac
done


# List of input files
[ -z $FILELIST ] && FILELIST=/data/tusers/zhongrenhu/for_TE_and_sSNV/rawdata/dna/BAM/downsample/downsample_filelist.txt
[ -z $OUTPUTDIR ] && OUTPUTDIR=/data/tusers/zhongrenhu/for_TE_and_sSNV/rawdata/dna/BAM/downsample

# Get files from FILELIST
echo "Reading filelist..."
str=`sed "${SLURM_ARRAY_TASK_ID}q;d" $FILELIST`
# Read the SLURM_ARRAY_TASK_ID'th line from filelist
# Split line into array
toks=($str)
INPUT=${toks[0]}
DEPTH=${toks[1]}
echo ""

# Create temporary directory for work
echo "Creating temporary working dir..."
TMPDIR=$(mktemp -d -t zhongrenhu-tmp-XXXXXXXX)
cd $TMPDIR
echo "Changing wd: " $(pwd)
echo ""

# Copy inputs to temp dir
echo "Copying input files to temp directory..."
cp $INPUT .
echo "Done."
echo ""

# Get filenames in temp by cutting of the path
BAM=$(basename $INPUT)
PREFIX=${BAM%.bam}

# exporting environment variable
export PATH=/data/tusers/zhongrenhu/Software/TEMP2/bin:$PATH

# Process the data
echo "downsample for ${BAM}..."
for seed in ${SEEDS[*]}
do
    echo "downsampling 10X..." && date
    fraction=$(awk -v depth=${DEPTH} 'BEGIN{printf "%u", (10/depth)*100}')
    samtools view -b -1 -@ 56 -s ${seed}.${fraction} -o ./${PREFIX}.10X.seed${seed}.bam ${BAM}
    echo ""

    echo "downsampling 20X..." && date
    fraction=$(awk -v depth=${DEPTH} 'BEGIN{printf "%u", (20/depth)*100}')
    samtools view -b -1 -@ 56 -s ${seed}.${fraction} -o ./${PREFIX}.20X.seed${seed}.bam ${BAM}
    echo ""

    echo "downsampling 30X..." && date
    fraction=$(awk -v depth=${DEPTH} 'BEGIN{printf "%u", (30/depth)*100}')
    samtools view -b -1 -@ 56 -s ${seed}.${fraction} -o ./${PREFIX}.30X.seed${seed}.bam ${BAM}
    echo ""
done
echo "Done." && date
echo ""

# Copy files to output dir
echo "Copying results to destination..."
ls -lh
cp *seed* $OUTPUTDIR/
echo "Done."
echo ""


# Clean up
echo "Cleaning up..."
cd $HOME
echo "Deleting temp dir: " $TMPDIR
rm -rd $TMPDIR
echo "/tmp contents:"
ls -lh /tmp
echo ""
echo "Script complete."
date