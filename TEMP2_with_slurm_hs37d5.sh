#!/bin/bash
#SBATCH --nodes=1
#SBATCH --time=120:00:00
#SBATCH --mem=200G
#SBATCH -c 64
#SBATCH --array=1
#SBATCH --partition=5days
#SBATCH --output=/data/tusers/zhongrenhu/for_SMS/dna_pipline/TEMP2_result/logs/TEMP2-log-%A-%a.out

# Print a little info for debugging
echo "HOSTNAME: " $(hostname)
echo "SLURM_JOB_NODELIST: " $SLURM_JOB_NODELIST
echo "SLURM_JOBID: " $SLURM_JOBID
echo "SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID
echo "SLURM_ARRAY_JOB_ID: " $SLURM_ARRAY_JOB_ID
date
echo ""


# List of input files
[ -z $FILELIST ] && FILELIST=/data/tusers/zhongrenhu/for_SMS/dna_pipline/TEMP2_result/TEMP2_filelist1
[ -z $OUTPUTDIR ] && OUTPUTDIR=/data/tusers/zhongrenhu/for_SMS/dna_pipline/TEMP2_result

# Get files from FILELIST
echo "Reading filelist..."
str=`sed "${SLURM_ARRAY_TASK_ID}q;d" $FILELIST`
# Read the SLURM_ARRAY_TASK_ID'th line from filelist
# Split line into array
toks=($str)
LEFT=${toks[0]}
RIGHT=${toks[1]}

# Create temporary directory for work
echo "Creating temporary working dir..."
TMPDIR=$(mktemp -d -t zhongrenhu-tmp-XXXXXXXX)
cd $TMPDIR
echo "Changing wd: " $(pwd)
echo ""

# Copy inputs to temp dir
echo "Copying input files to temp directory..."
PREFIX=`basename ${LEFT%_1.fastq.gz}` && INPUT_DIR=`dirname ${INPUT}`
cp $LEFT $RIGHT .
echo "Done."
echo ""


# exporting environment variable and default parameter
export PATH=/data/tusers/zhongrenhu/Software/gawk-5.1.0/opt/bin:$PATH
export PATH=/data/tusers/zhongrenhu/Software/TEMP2:$PATH
export PATH=/data/tusers/zhongrenhu/Software/TEMP2/bin:$PATH
export PATH=/data/tusers/zhongrenhu/Software/anaconda3/envs/python27/bin/:$PATH
PATH_ANNO="/data/tusers/yutianx/tongji2/GitHuB/piSet/annotation/hs37d5"
BWA_INDEX=${PATH_ANNO}/BWAIndex/genome
GENOME_FA=${PATH_ANNO}/hs37d5.fa
TE_FA=${PATH_ANNO}/ALUL1SVA.fa
TE_BED6=${PATH_ANNO}/ALUL1SVA.bed

# Process the data
echo "Running TEMP2 insertion for ${BAM}"
mkdir ${PREFIX}
TEMP2 insertion -l ${LEFT} -r ${RIGHT} -I ${BWA_INDEX} -g ${GENOME_FA} -R ${TE_FA} -t ${TE_BED6} -o ${PREFIX} -p ${PREFIX} -c 32 -M 2 -m 10 -U 0.95 -S -A -d
echo ""

# Copy files to output dir
echo "Copying results to destination..."
ls -lh
rm ${LEFT} ${RIGHT}
cp -r ${PREFIX} $OUTPUTDIR/
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
