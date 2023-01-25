#!/bin/bash

module load eb/2017  GCC/5.4.0-2.27  OpenMPI/2.0.0 cutadapt/1.9.1-Python-2.7.12

set -e
set -o pipefail

#############################################################################
# Creating a usage sentence
usage="To use this viral alignment for initial sample, mate and softclip reads:
$(basename $0) [-i path-to-orig-fastq] [-e fastq|fq|fastq.gz|fq.gz] [-r single|paired] [-m path-to-folder-containing-mate-and-softclip-folders] 
[-v path-to-viral-genome] [-w this is genewerk data] [-h] usage info"

## Printing it when there are no arguments
if [ "$#" == "0" ] || [ "$1" == "--help" ]; then
        printf -- "\033[33m$usage. \033[0m\n" >&2
        exit
fi

while getopts 'i:e:r:m:v:w' OPTION; do
	case "$OPTION" in
		i) indir="$OPTARG"
		printf -- "\033[32m The original fastq files are at $OPTARG. \033[0m\n"
		;;
		e) suffix="$OPTARG"
		printf -- "\033[32m The suffix of your files is $OPTARG. \033[0m\n"
		;;
		r) reads="$OPTARG"
		printf -- "\033[32m Your readtype is $OPTARG end. \033[0m\n"
		;;
		m) readdir="$OPTARG"
		printf -- "\033[32m The samples aligned to host genome are at $OPTARG. \033[0m\n"
		;;
		v) vir_genome="$OPTARG"
		printf -- "\033[32m The viral genome is $OPTARG. \033[0m\n"
		;;
		w) gw=1
		printf -- "\033[32m You've indicated that this is genewerk data. \033[0m\n"
		;;
	esac
done

# Checking for missings

if [ "x" == "x$indir" ]; then
  printf -- "\033[31m [-i path-to-data-dir] option is required. Please select the directory that contains your fastq files.\033[0m\n"
  exit
fi

if [ "x" == "x$suffix" ]; then
  printf -- "\033[31m [-e fastq|fq|fastq.gz|fq.gz] option is required. Please indicate the suffix of your original fastq files.\033[0m\n"
  exit
fi

if [ "x" == "x$reads" ]; then
  printf -- "\033[31m [-r single|paired] option is required. Please select the read type of your sequencng.\033[0m\n"
  exit
fi

if [ "x" == "x$readdir" ]; then
  printf -- "\033[31m [-m path-to-folder-containing-mate-and-softclip-folders] option is required. Please select the parent folder of the mate and softclip folders.\033[0m\n"
  exit
fi

if [ "x" == "x$vir_genome" ]; then
  printf -- "\033[31m [-v path-to-viral-genome] option is required. Please direct this option to a bwa indexed viral vector genome.\033[0m\n"
  exit
fi


# Defining tools
my_bwa=/path-to-program/bwa-0.7.17/bwa
my_stools=/path-to-program/samtools-1.9/samtools
my_sblaster=/path-to-program/samblaster/samblaster
my_bedtools=/path-to-program/bedtools2/bin/bedtools
## Here we're using my_clipm to match at least 30 bp of viral sequence
my_clipm=/path-to-scripts/v.int/samclip2_Pfmod

# Making the new full viralign diractory
mkdir -p ${readdir}/full_viral/

if [ "$gw" == 1 ]; then

  cd ${indir}
  
else
  # Mapping mate
  cd ${readdir}/mate/
  
  ls -1 *.fq | sed 's/\_unaligned\.fq//g' | uniq > $readdir/mate/sample_list.txt
  for i in $(cat $readdir/mate/sample_list.txt); do
    # sample and read group
    fq1=${readdir}/mate/${i}_unaligned.fq
    rg=$(echo \"\@RG\\tID:${i}\\tPL:Illumina\\tPU:x\\tLB:single\\tSM:${i}\")
  
    # Submitting alignment
    bsub -q medium -app large -J "v.int.mate-${i}" -n 4,4 -M 16GB -R "span[hosts=1] rusage[mem=16GB]" \
    -o "${readdir}/mate/${i}.v.int.2.log" -e "${readdir}/mate/${i}.v.int.2.err" \
      "${my_bwa} mem ${vir_genome} -t 4 -R ${rg} ${fq1} | \
      ${my_stools} view -h -F 4 -q 30 - | \
      ${my_clipm} \
      --max 120 \
      --match 30 \
      --ref ${vir_genome} | \
      ${my_stools} sort -nT ${readdir}/mate/${i}_unaligned.prefix.bam - | \
      ${my_stools} fixmate -m - - | \
      ${my_stools} sort -T ${readdir}/mate/${i}_unaligned.postfix.bam - | \
      ${my_stools} markdup -r - -O BAM ${readdir}/mate/${i}_unaligned_mapped_filtered.bam"
      done
      
    # Mapping softclip reads
    cd ${readdir}/softclip/
    
    ls -1 *.fq | sed 's/\_softclipped\.fq//g' | uniq > $readdir/softclip/sample_list.txt

    for i in $(cat $readdir/softclip/sample_list.txt); do
      # sample and read group
      fq1=${readdir}/softclip/${i}_softclipped.fq
      rg=$(echo \"\@RG\\tID:${i}\\tPL:Illumina\\tPU:x\\tLB:single\\tSM:${i}\")
  
      # Submitting alignment
      bsub -q medium -app large -J "v.int.softclip-${i}" -n 4,4 -M 16GB -R "span[hosts=1] rusage[mem=16GB]" \
      -o "${readdir}/softclip/${i}.v.int.2.log" -e "${readdir}/softclip/${i}.v.int.2.err" \
      "${my_bwa} mem ${vir_genome} -t 4 -R ${rg} ${fq1} | \
      ${my_stools} view -h -F 4 -q 30 - | \
      ${my_clipm} \
      --max 120 \
      --match 30 \
      --ref ${vir_genome} |\
      ${my_stools} sort -nT ${readdir}/softclip/${i}_softclipped.prefix.bam - | \
      ${my_stools} fixmate -m - - | \
      ${my_stools} sort -T ${readdir}/softclip/${i}_softclipped.postfix.bam - | \
      ${my_stools} markdup -r - -O BAM ${readdir}/softclip/${i}_softclipped_mapped_filtered.bam"
  done
fi

# Doing alignments of raw sequence
## Creating a sample list

cd $indir

if [ "$reads" == "paired" ]; then
  ls -1 *${suffix} | sed -e "s/\_1.*//g" | sed -e "s/\_2.*//g" |  uniq > $indir/sample_list.txt
  else
  ls -1 *${suffix} | sed "s/\.$suffix//g" | uniq > $indir/sample_list.txt
fi

list=$indir/sample_list.txt

## Checking the total number of samples

tot=$(cat $indir/sample_list.txt | wc -l)

if [ "$tot" -eq "0" ]; then
	printf -- "\033[31m There were no samples of your specified suffix detected in your data directory. Please check if you have your fastq.gz or fq.gz files in the sample directory. If so, please make sure that you have specified fastq or fq in the -e option above.\033]0m\n"
fi


for f in $(cat $list); do
    
  if [ "$reads" == "single" ]; then
    
    # Defining fastq to use
    fq1=${f}.${suffix}
    
    if [ "$gw" == 1 ]; then
    # Defining RG code
    rg=$(echo \"\@RG\\tID:${f}\\tPL:Illumina\\tPU:x\\tLB:single\\tSM:${f}\")
    
    # Running pipeline
    bsub -q medium -app large -J "v.int.viralign-${f}" -n 4,4 -M 16GB -R "span[hosts=1] rusage[mem=16GB]" \
    -o "${readdir}/full_viral/${f}.v.int.viralign.log" -e "${readdir}/full_viral/${f}.v.int.viralign.err" \
    "cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC  -m 40 -q 30,30 ${indir}/$fq1 | \
    ${my_bwa} mem ${vir_genome} -t 4 -R ${rg} - | \
    ${my_stools} view -h -q 30 -F 4 - | \
    ${my_clipm} --max 350 --match 30 --ref ${vir_genome} | \
    ${my_stools} sort -nT ${readdir}/full_viral/${f}.prefix.bam - | \
    ${my_stools} fixmate -m - - | \
    ${my_stools} sort -T ${readdir}/full_viral/${f}.postfix.bam - | \
    ${my_stools} markdup -r - -O BAM ${readdir}/full_viral/${f}.bam"
    
      else
      
      # Defining RG code
      rg=$(echo \"\@RG\\tID:${f}\\tPL:Illumina\\tPU:x\\tLB:single\\tSM:${f}\")
    
      # Running pipeline
       bsub -q medium -app large -J "v.int.viralign-${f}" -n 4,4 -M 16GB -R "span[hosts=1] rusage[mem=16GB]" \
       -o "${readdir}/full_viral/${f}.v.int.viralign.log" -e "${readdir}/full_viral/${f}.v.int.viralign.err" \
       "cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC  -m 40 -q 30,30 ${indir}/$fq1 | \
       ${my_bwa} mem ${vir_genome} -t 4 -R ${rg} - | \
       ${my_stools} view -h -q 30 -F 4 - | \
       ${my_clipm} --max 120 --match 30 --ref ${vir_genome} | \
       ${my_stools} sort -nT ${readdir}/full_viral/${f}.prefix.bam - | \
       ${my_stools} fixmate -m - - | \
       ${my_stools} sort -T ${readdir}/full_viral/${f}.postfix.bam - | \
       ${my_stools} markdup -r - -O BAM ${readdir}/full_viral/${f}.bam"
    fi
    
    else
    
    # Submitting paired end script
    
    bsub -q medium -app large -J "v.int.viralign-${f}" -n 4,4 -M 16GB -R "span[hosts=1] rusage[mem=16GB]" \
    -o "${readdir}/full_viral/${f}.v.int.viralign.log" -e "${readdir}/full_viral/${f}.v.int.viralign.err" \
    "/hpc/grid/dsrd-invtox/workspace/scripts/v.int/step3.2.paired.sh \
    -g ${vir_genome} \
    -i ${indir} \
    -e ${suffix} \
    -o ${readdir}/full_viral \
    -f ${f} >& ${readdir}/full_viral/${f}.3.2.txt"
  fi
done

# Exit
printf -- '\n';
exit 0;