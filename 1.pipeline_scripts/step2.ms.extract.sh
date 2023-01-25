#!/bin/bash

module load eb/2017  GCC/5.4.0-2.27  OpenMPI/2.0.0 cutadapt/1.9.1-Python-2.7.12

set -e
set -o pipefail

#############################################################################
# Creating a usage sentence
usage="To use this post-alignment mate and softclip extraction pipeline please specify the following:
$(basename $0) [-i path-to-data-dir] [-g cyno|mouse] [-h] usage info"

## Printing it when there are no arguments
if [ "$#" == "0" ] || [ "$1" == "--help" ]; then
        printf -- "\033[33m$usage. \033[0m\n" >&2
        exit
fi


while getopts 'i:g:w' OPTION; do
	case "$OPTION" in
		i) indir="$OPTARG"
		printf -- "\033[32m Your input directory is: $OPTARG. \033[0m\n"
		;;
		g) species="$OPTARG"
		printf -- "\033[32m The genome you've chosen is $OPTARG. \033[0m\n"
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

if [ "x" == "x$species" ]; then
  printf -- "\033[31m [-g cyno|mouse] option is required. Please select the species from which this data was sequenced.\033[0m\n"
  exit
fi

# Starting the process
cd $indir

## Defining genome
if [ "$species" == "cyno" ]; then
  genome=/path-to-genome/Macaca_fascicularis.Macaca_fascicularis_5.0.dna.toplevel.fa
  elif [ "$species" == "human" ]; then
  genome=/path-to-genome/Homo_sapiens.GRCh38.dna.toplevel.fa
  else
  genome=/path-to-genome/Mus_musculus.GRCm38.dna.toplevel.fa
fi

if [ "$gw" == 1 ]; then

  cd ${indir}
  
else

  ## Defining tools
  my_stools=/path-to-program/samtools-1.9/samtools
  my_bedtools=/path-to-program/bedtools2/bin/bedtools
  # Here my-clip is used to extract reads reads that have clips, so appropriate to use the softclip filter, rather than the match filter
  my_clipm=/path-to-scripts/v.int/samclip2_Pfmod

  # Making directories
  mkdir -p ${indir}/mate
  mkdir -p ${indir}/softclip

  # Defining directories
  mate_out=${indir}/mate
  soft_out=${indir}/softclip

  # Extracting mated reads

  for sample in $(ls -1 *bam | grep -v "sorted" | grep -v "merged" | sed 's/\.bam//g'); do
    # Singly unmapped reads (remember to never include q30 flag here)
    bsub -q medium -app large -J "v.int.mates${sample}" -n 4,4 -M 16GB -R "span[hosts=1] rusage[mem=16GB]" \
    -o "${mate_out}/${sample}.v.int.2.log" -e "${mate_out}/${sample}.v.int.2.err" \
    "${my_stools} view -h -f 4 -F 8 ${indir}/${sample}.bam |\
    ${my_bedtools} bamtofastq -i - -fq ${mate_out}/${sample}_unaligned.fq"
    
    # Mates of unaligned reads
    bsub -q medium -app large -J "v.int.mates${sample}" -n 4,4 -M 16GB -R "span[hosts=1] rusage[mem=16GB]" \
    -o "${mate_out}/${sample}.v.int.2.log" -e "${mate_out}/${sample}.v.int.2.err" \
    "${my_stools} view -h -f 8 -F 4 -q 30 ${indir}/${sample}.bam > ${mate_out}/mapped_mates_${sample}.bam"
    
    # Softclipped reads
    bsub -q medium -app large -J "v.int.soft${sample}" -n 4,4 -M 16GB -R "span[hosts=1] rusage[mem=16GB]" \
    -o "${soft_out}/${sample}.v.int.2.log" -e "${soft_out}/${sample}.v.int.2.err" \
    "$my_stools view -h ${indir}/${sample}.bam |\
    $my_clipm \
    --match 30 \
    --invert \
    --max 30 \
    --ref $genome |\
    $my_stools sort > ${soft_out}/${sample}_softclipped.bam
    
    # Converting bam to fq
    ${my_bedtools} bamtofastq -i ${soft_out}/${sample}_softclipped.bam -fq ${soft_out}/${sample}_softclipped.fq
    "
    
  done
fi

# Exit
printf -- '\n';
exit 0;
