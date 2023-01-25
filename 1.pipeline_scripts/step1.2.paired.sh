#!/bin/bash

module load eb/2017  GCC/5.4.0-2.27  OpenMPI/2.0.0 cutadapt/1.9.1-Python-2.7.12

set -e
set -o pipefail

#############################################################################
# Creating a usage sentence
usage="To use this stringent alignment pipeline please specify the following:
$(basename $0) [-g cyno|mouse] [-i path-to-data-dir] [-e fastq|fq|fastq.gz|fq.gz] [-o path-do-output-dir] [-f sample name] [-r rg tag] [-h] usage info"

## Printing it when there are no arguments
if [ "$#" == "0" ] || [ "$1" == "--help" ]; then
printf -- "\033[33m$usage. \033[0m\n" >&2
exit
fi


while getopts 'g:i:e:o:f:r:' OPTION; do
case "$OPTION" in
g) genome="$OPTARG"
printf -- "\033[32m The genome you've chosen is $OPTARG. \033[0m\n"
;;
i) indir="$OPTARG"
printf -- "\033[32m Your input directory is: $OPTARG. \033[0m\n"
;;
e) suffix="$OPTARG"
printf -- "\033[32m The suffix of your files is $OPTARG. \033[0m\n"
;;
o) outdir="$OPTARG"
printf -- "\033[32m Your output directory is $OPTARG. \033[0m\n"
;;
f) f="$OPTARG"
printf -- "\033[32m The sample name is $OPTARG. \033[0m\n"
;;
esac
done

# Checking for missings

if [ "x" == "x$genome" ]; then
printf -- "\033[31m [-g cyno|mouse] option is required. Please select the genome to use.\033[0m\n"
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

if [ "x" == "x$outdir" ]; then
printf -- "\033[31m [-o path-do-output-dir] option is required. Please select the output directory.\033[0m\n"
exit
fi


if [ "x" == "x$f" ]; then
printf -- "\033[31m [-f sample name] is required. Please put in a sample name. \033[0m\n"
exit
fi


# Defining tools
my_bwa=/path-to-program/bwa-0.7.17/bwa
my_stools=/path-to-program/samtools-1.9/samtools

# Defining RG code
rg=$(echo \@RG\\tID:${f}\\tPL:Illumina\\tPU:x\\tLB:$reads\\tSM:${f})

# Defining fqs

fq1=${f}_1.${suffix}
fq2=${f}_2.${suffix}

# Echo them all

# Code to run

paste <(zcat ${indir}/$fq1 | paste - - - -) \
<(zcat ${indir}/$fq2 | paste - - - -) | \
tr '\t' '\n' | \
cutadapt --interleaved -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTA -m 40 -q 30,30 - | \
${my_bwa} mem ${genome} -t 16 -R ${rg} -p - | \
${my_stools} view -h -u - | \
${my_stools} sort -nT ${outdir}/${f}.prefix.bam - | \
${my_stools} fixmate -m - - | \
${my_stools} sort -T ${outdir}/${f}.postfix.bam - | \
${my_stools} markdup -r - -O BAM ${outdir}/${f}.bam

# Exit
printf -- '\n';
exit 0;