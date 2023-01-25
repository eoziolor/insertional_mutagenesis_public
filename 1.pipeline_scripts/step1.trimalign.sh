#!/bin/bash

module load eb/2017  GCC/5.4.0-2.27  OpenMPI/2.0.0 cutadapt/1.9.1-Python-2.7.12

set -e
set -o pipefail

#############################################################################
# Creating a usage sentence
usage="To use this stringent alignment pipeline please specify the following:
$(basename $0) [-g cyno|mouse] [-i path-to-data-dir] [-e fastq|fq|fastq.gz|fq.gz] [-r single|paired] [-o path-do-output-dir] [-w this is genewerk data] [-f kick off full wg or tse pipeline] [-h] usage info"

## Printing it when there are no arguments
if [ "$#" == "0" ] || [ "$1" == "--help" ]; then
printf -- "\033[33m$usage. \033[0m\n" >&2
exit
fi


while getopts 'g:i:e:r:o:wf' OPTION; do
case "$OPTION" in
g) species="$OPTARG"
printf -- "\033[32m The genome you've chosen is $OPTARG. \033[0m\n"
;;
i) indir="$OPTARG"
printf -- "\033[32m Your input directory is: $OPTARG. \033[0m\n"
;;
e) suffix="$OPTARG"
printf -- "\033[32m The suffix of your files is $OPTARG. \033[0m\n"
;;
r) reads="$OPTARG"
printf -- "\033[32m Your readtype is $OPTARG end. \033[0m\n"
;;
o) outdir="$OPTARG"
printf -- "\033[32m Your output directory is $OPTARG. \033[0m\n"
;;
w) gw=1
printf -- "\033[32m You've indicated that this is genewerk data. \033[0m\n"
;;
f) full=1
printf -- "\033[32m You've selected to kick off the full gw/tse pipeline. \033[0m\n"
esac
done

# Checking for missings

if [ "x" == "x$species" ]; then
printf -- "\033[31m [-g cyno|mouse] option is required. Please select the species from which this data was sequenced.\033[0m\n"
exit
fi

if [ "x" == "x$indir" ]; then
printf -- "\033[31m [-i path-to-data-dir] option is required. Please select the directory that contains your fastq files.\033[0m\n"
exit
fi

if [ "x" == "x$suffix" ]; then
printf -- "\033[31m [-e fastq|fq|fastq.gz|fq.gz] option is required. Please select the suffix of your fastq files.\033[0m\n"
exit
fi

if [ "x" == "x$reads" ]; then
printf -- "\033[31m [-r single|paired] option is required. Please select the read type of your sequencng.\033[0m\n"
exit
fi


if [ "x" == "x$outdir" ]; then
printf -- "\033[31m [-o path-do-output-dir] option is required. Please select the output directory.\033[0m\n"
exit
fi

# Starting the process

mkdir -p ${outdir}

## Defining genome
if [ "$species" == "cyno" ]; then
genome=/path-to-genome/Macaca_fascicularis.Macaca_fascicularis_5.0.dna.toplevel.fa
elif [ "$species" == "human" ]; then
genome=/path-to-genome/Homo_sapiens.GRCh38.dna.toplevel.fa
else
  genome=/path-to-genome/Mus_musculus.GRCm38.dna.toplevel.fa
fi

# Defining tools
my_bwa=/path-to-program/bwa-0.7.17/bwa
my_stools=/path-to-program/samtools-1.9/samtools
my_sblaster=/path-to-program/samblaster/samblaster
my_bedtools=/path-to-program/bedtools2/bin/bedtools

## Creating a sample list

cd $indir

if [ "$reads" == "paired" ]; then
ls -1 *${suffix} | sed -e "s/\_1.*//g" | sed -e "s/\_2.*//g" |  uniq > $indir/sample_list.txt
else
  ls -1 *${suffix} | sed "s/\.$suffix//g"| uniq > $indir/sample_list.txt
fi

list=$indir/sample_list.txt

# Checking the total number of samples

tot=$(cat $indir/sample_list.txt | wc -l)

if [ "$tot" -eq "0" ]; then
printf -- "\033[31m There were no samples of your specified suffix detected in your data directory. Please check if you have your fastq.gz or fq.gz files in the sample directory. If so, please make sure that you have specified fastq or fq in the -e option above.\033]0m\n"
fi

# Trimming and aligning samples to correct host

for f in $(cat $list); do
# Defining RG code
rg=$(echo \"\@RG\\tID:${species}_${f}\\tPL:Illumina\\tPU:x\\tLB:$reads\\tSM:${species}_${f}\")
    
  if [ "$reads" == "single" ]; then
    
    # Defining fastq to use
    fq1=${f}.${suffix}
    
    if [ "$gw" == 1 ]; then
    
    # Running pipeline
    bsub -q medium -app large -J "v.int.ta-${f}" -n 16,16 -M 64GB -R "span[hosts=1] rusage[mem=64GB]" \
    -o "${outdir}/${f}.v.int.1.log" -e "${outdir}/${f}.v.int.1.err" \
    "cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC  -m 40 -q 30,30 ${indir}/$fq1 | \
     ${my_bwa} mem ${genome} -t 16 -R ${rg} - | \
     ${my_stools} view -h -u -q 30 -F 4 - | \
     ${my_stools} sort -nT ${outdir}/${f}.prefix.bam - | \
     ${my_stools} fixmate -m - - | \
     ${my_stools} sort -T ${outdir}/${f}.postfix.bam - | \
     ${my_stools} markdup -r - -O BAM ${outdir}/${f}.bam"
    
      else
    
      # Running pipeline
       bsub -q long -app large -J "v.int.ta-${f}" -n 16,16 -M 64GB -R "span[hosts=1] rusage[mem=64GB]" \
       -o "${outdir}/${f}.v.int.1.log" -e "${outdir}/${f}.v.int.1.err" \
       "cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC  -m 40 -q 30,30 ${indir}/$fq1 | \
     ${my_bwa} mem ${genome} -t 16 -R ${rg} - | \
     ${my_stools} view -h -u -q 30 -F 4 - | \
     ${my_stools} sort -nT ${outdir}/${f}.prefix.bam - | \
     ${my_stools} fixmate -m - - | \
     ${my_stools} sort -T ${outdir}/${f}.postfix.bam - | \
     ${my_stools} markdup -r - -O BAM ${outdir}/${f}.bam"
    fi
    
    else
    
    # Submitting paired end script
    
    bsub -q long -app large -J "v.int.ta-${f}" -n 16,16 -M 64GB -R "span[hosts=1] rusage[mem=64GB]" \
    -o "${outdir}/${f}.v.int.1.log" -e "${outdir}/${f}.v.int.1.err" \
    "/hpc/grid/dsrd-invtox/workspace/scripts/v.int/step1.2.paired.sh \
     -g ${genome} \
     -i ${indir} \
     -e ${suffix} \
     -o ${outdir} \
     -f ${f} >& $outdir/err1.2.txt"
  fi
done

if [ "$full" == 1 ]; then
  bsub -q express -app large -w 'ended("v.int.ta-*")' -J "v.int.step2" -n 4,4 -M 32GB -R "span[hosts=1] rusage[mem=32GB]" \
  -o "${outdir}/st1.to.st2.log" -e "${outdir}/st1.to.st2.err" \
  "/hpc/grid/dsrd-invtox/workspace/scripts/v.int/step2.ms.extract.sh \
     -i $outdir -g $species >& $outdir/err2.txt"
fi

# Exit
printf -- '\n';
exit 0;
