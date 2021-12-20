#!/bin/bash
#SBATCH --nodes=1
#SBATCH --time=12:00:00
#SBATCH --mem=220G
#SBATCH -c 64
#SBATCH --array=1-4
#SBATCH --partition=12hours
#SBATCH --output=/data/tusers/zhongrenhu/for_SMS/dna_pipline/NANOVAR_result/logs/nanovar-log-%A-%a.out

# Print a little info for debugging
echo "HOSTNAME: " $(hostname)
echo "SLURM_JOB_NODELIST: " $SLURM_JOB_NODELIST
echo "SLURM_JOBID: " $SLURM_JOBID
echo "SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID
echo "SLURM_ARRAY_JOB_ID: " $SLURM_ARRAY_JOB_ID
date
echo ""


# Parameter initialization
# List of input files, I make this with ls /path/to/folder/*.suffix | sort | pr -1aT -s' ' > filelist.txt
[ -z $FILELIST ] && FILELIST=/data/tusers/zhongrenhu/for_SMS/dna_pipline/NANOVAR_result/nanovar_filelist
[ -z $THREADS ] && THREADS=24


# Get files from FILELIST
echo "Reading filelist..."
str=`sed "${SLURM_ARRAY_TASK_ID}q;d" $FILELIST`
# Read the SLURM_ARRAY_TASK_ID'th line from filelist
# Split line into array
toks=($str)
INPUT=${toks[0]}
GAP_FILE=${toks[1]}
REF_FASTA=${toks[2]}


# Create temporary directory for work
echo "Creating temporary working dir..."
TMPDIR=$(mktemp -d -t zhongrenhu-tmp-XXXXXXXX)
cd $TMPDIR
echo "Changing wd: " $(pwd)
echo ""


# Copy inputs to temp dir
echo "Copying input file ${INPUT} to temp directory..."
cp $INPUT .
echo "Done."
echo ""


# Get filenames in temp by cutting of the path
PREFIX=`basename ${INPUT%.f*q*}` && LOCAL_INPUT=`basename $INPUT`
[ -z $OUTPUTDIR ] && OUTPUTDIR=/data/tusers/zhongrenhu/for_SMS/dna_pipline/NANOVAR_result/$PREFIX && mkdir $OUTPUTDIR


# exporting environment variable and default parameter
export PATH=/data/tusers/zhongrenhu/Software/anaconda3/bin:$PATH
export PATH=/data/tusers/zhongrenhu/Software/bin:$PATH


# Process the data
echo "Running nanovar for ${PREFIX}"
nanovar -f ${GAP_FILE} -t ${THREADS} ${LOCAL_INPUT} ${REF_FASTA}
echo ""

# Copy files to output dir
echo "Copying results to destination..."
ls -lh
cp -r ./* $OUTPUTDIR/
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