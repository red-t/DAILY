#!/bin/bash
#SBATCH --nodes=1
#SBATCH --time=120:00:00
#SBATCH --mem=220G
#SBATCH -c 64
#SBATCH --array=1-2
#SBATCH --partition=5days
#SBATCH --output=/data/tusers/zhongrenhu/for_SMS/dna_pipline/TLDR_result/logs/minimap2-log-%A-%a.out

# Print a little info for debugging
echo "HOSTNAME: " $(hostname)
echo "SLURM_JOB_NODELIST: " $SLURM_JOB_NODELIST
echo "SLURM_JOBID: " $SLURM_JOBID
echo "SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID
echo "SLURM_ARRAY_JOB_ID: " $SLURM_ARRAY_JOB_ID
date
echo ""


######Help Information########
# help information
function help_info(){
    echo `basename $0`" [options] -R reference.fa -Q query.fastq"
    echo -e "\t-R <ref.fa>\treference FASTA"
    echo -e "\t-Q <query.fastq>\tquery FASTQ"
    echo -e "\t-t <threads>\tnumber of threads [3]"
    echo -e "\t-a output in the SAM format (PAF by default)"
    echo -e "\t-t <threads>\tnumber of threads [3]"
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

while getopts ":R:Q:at:x:h" OPTION; do
    case $OPTION in
        R)  REF_FASTA=($OPTARG);;
        Q)  QUERY_FASTQ=($OPTARG);;
        a)  SAM_FORMAT=1;;
        t)  THREADS=$OPTARG;;
        x)  PRESET=$OPTARG;;
        h)  help_info && exit 1;;
        *)  help_info && exit 1;;
    esac
done

# Parameter initialization
# List of input files, I make this with ls /path/to/folder/*.suffix | sort | pr -1aT -s' ' > filelist.txt
[ -z $FILELIST ] && FILELIST=/data/tusers.ds/zhongrenhu/for_SMS/dna_pipline/TLDR_result/minimap2_filelist
[ -z $SAM_FORMAT ] && SAM_FORMAT=1
[ -z $THREADS ] && THREADS=32

# Get files from FILELIST
echo "Reading filelist..."
str=`sed "${SLURM_ARRAY_TASK_ID}q;d" $FILELIST`
# Read the SLURM_ARRAY_TASK_ID'th line from filelist
# Split line into array
toks=($str)
REF_FASTA=${toks[0]}
FASTQ=${toks[1]}
PRESET=${toks[2]}

# Create temporary directory for work
echo "Creating temporary working dir..."
TMPDIR=$(mktemp -d -t zhongrenhu-tmp-XXXXXXXX)
cd $TMPDIR
echo "Changing wd: " $(pwd)
echo ""

# Copy inputs to temp dir
echo "Copying input file ${FASTQ} to temp directory..."
cp $FASTQ .
echo "Done."
echo ""

# Get filenames in temp by cutting of the path
PREFIX=`basename ${FASTQ%.f*q*}` && QUERY_FASTQ=`basename $FASTQ`
[ -z $OUTPUTDIR ] && OUTPUTDIR=`dirname $FASTQ`

# exporting environment variable and default parameter
export PATH=/data/tusers/zhongrenhu/Software/anaconda3/bin:$PATH
export PATH=$PATH:/data/tusers/zhongrenhu/Software/minimap2/
export PATH=$PATH:/data/tusers/zhongrenhu/Software/exonerate/bin/
export PATH=$PATH:/data/tusers/zhongrenhu/Software/mafft-7.487-without-extensions/bin
export PATH=$PATH:/data/tusers/zhongrenhu/Software/tldr/tldr/

# Process the data
echo "Running minimap2 for ${PREFIX}"
if [ -n $SAM_FORMAT ];then
    minimap2 -ax $PRESET --MD -t $THREADS $REF_FASTA $QUERY_FASTQ > $PREFIX.sam
    samtools view -bhS -@ $THREADS $PREFIX.sam | samtools sort -@ $THREADS -o $PREFIX.bam -
    else
    minimap2 -x $PRESET --MD -t $THREADS $REF_FASTA $QUERY_FASTQ > $PREFIX.paf
fi
echo ""

# Copy files to output dir
echo "Copying results to destination..."
ls -lh
rm ./*fastq*
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