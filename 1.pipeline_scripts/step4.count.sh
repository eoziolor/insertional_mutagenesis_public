#!/bin/bash

module load eb/2017  GCC/5.4.0-2.27  OpenMPI/2.0.0 cutadapt/1.9.1-Python-2.7.12

set -e
set -o pipefail

#############################################################################
# Creating a usage sentence
usage="To use this viral alignment for initial sample, mate and softclip reads:
$(basename $0) [-g cyno|mouse|human] [-i path-to-folder-containing-mate-and-softclip-folders] [-w this is genewerk data] [-h] usage info"

## Printing it when there are no arguments
if [ "$#" == "0" ] || [ "$1" == "--help" ]; then
        printf -- "\033[33m$usage. \033[0m\n" >&2
        exit
fi

while getopts 'g:i:w' OPTION; do
	case "$OPTION" in
		g) species="$OPTARG"
		printf -- "\033[32m The genome you've chosen is $OPTARG. \033[0m\n"
		;;
		i) indir="$OPTARG"
		printf -- "\033[32m The sample directory is $OPTARG. \033[0m\n"
		;;
		w) gw=1
		printf -- "\033[32m You've indicated that this is genewerk data. \033[0m\n"
		;;
	esac
done

# Checking for missings
if [ "x" == "x$species" ]; then
  printf -- "\033[31m [-g cyno|mouse|human] option is required. Please select the species from which this data was sequenced.\033[0m\n"
  exit
fi

if [ "x" == "x$indir" ]; then
  printf -- "\033[31m [-i path-to-data-dir] option is required. Please select the directory that contains your fastq files.\033[0m\n"
  exit
fi

## Defining genome
if [ "$species" == "cyno" ]; then
  genome=/path-to-genome/Macaca_fascicularis.Macaca_fascicularis_5.0.dna.toplevel.fa
  elif [ "$species" == "human" ]; then
  genome=/path-to-genome/Homo_sapiens.GRCh38.dna.toplevel.fa
  else
  genome=/path-to-genome/Mus_musculus.GRCm38.dna.toplevel.fa
fi

# Defining tools
my_stools=/path-to-program/samtools-1.9/samtools
my_sblaster=/path-to-program/samblaster/samblaster
my_bedtools=/path-to-program/bedtools2/bin/bedtools
my_clipm=/path-to-scripts/v.int/samclip2_Pfmod
my_extract=/path-to-program/extract_reads.py

# first part is extracting the reads

# Starting with mates
cd ${indir}/mate

## Resorting and filtering mate files

for f in $(ls -1 mapped_mates*); do
nam=$(echo $f | sed 's/mapped\_mates\_//g' | sed 's/\.bam//g')

bsub -q medium -app large -J "${nam}-mate-resort" -n 4,4 -M 16GB -R "span[hosts=1] rusage[mem=16GB]" -o "${indir}/mate/${nam}_mate_resort.log" -e "${indir}/mate/${nam}_mate_resort.err" \
"$my_stools view -h ${indir}/mate/${f} | \
${my_clipm} --max 120 --match 30 --ref ${genome} | \
$my_stools sort > ${indir}/mate/${nam}_sorted_mates.bam"

done

## extracting virally maped reads and filtering

for f in $(ls -1 *unaligned_mapped_filtered.bam); do
nam=$(echo $f | sed 's/\.bam//g')

bsub -q express -app large -w 'ended("*mate-resort")' -J "mate_read_extract_${f}" -n 4,4 -M 16GB -R "span[hosts=1] rusage[mem=16GB]" \
    -o "${indir}/mate/${nam}.mate_read_extract.log" -e "${indir}/mate/${nam}.mate_read_extract.err" \
    "${my_stools} view ${indir}/mate/${f} | cut -f 1 | sort | uniq > ${indir}/mate/${nam}_reads.txt"

bsub -q express -app large -w 'ended("*mate-resort")' -J "mate_read_info_extract_${f}" -n 4,4 -M 16GB -R "span[hosts=1] rusage[mem=16GB]" \
    -o "${indir}/mate/${nam}.mate_read_info_extract.log" -e "${indir}/mate/${nam}.mate_read_info_extract.err" \
    "${my_stools} view ${indir}/mate/${f} | cut -f 1,3,4,6 | sort | uniq > ${indir}/mate/${nam}_reads_info.txt"
done

## extracitng viral mates in a bam

for f in $(ls -1 mapped_mates*); do
nam=$(echo $f | sed 's/mapped\_mates\_//g' | sed 's/\.bam//g')

bsub -q medium -app large -w 'ended("mate_read_extract*")' -J "mate_bam_extract_${nam}" -n 4,4 -M 16GB -R "span[hosts=1] rusage[mem=16GB]" \
    -o "${indir}/mate/${nam}.mate_bam_extract.log" -e "${indir}/mate/${nam}.mate_bam_extract.err" \
    "python ${my_extract} -b ${indir}/mate/${nam}_sorted_mates.bam -n ${indir}/mate/${nam}_unaligned_mapped_filtered_reads.txt -o ${indir}/mate/${nam}_viral_mate.bam"

done

# Extracting softclips
cd ${indir}/softclip/

for f in $(ls -1 *softclipped_mapped_filtered.bam); do
nam=$(echo $f | sed 's/\.bam//g')

bsub -q express -app large -J "softclip_read_extract_${f}" -n 4,4 -M 16GB -R "span[hosts=1] rusage[mem=16GB]" \
    -o "${indir}/softclip/${nam}.softclip_read_extract.log" -e "${indir}/softclip/${nam}.softclip_read_extract.err" \
    "$my_stools view ${indir}/softclip/${f} | cut -f 1 > ${indir}/softclip/${nam}_reads.txt"


bsub -q express -app large -J "softclip_read_info_extract_${f}" -n 4,4 -M 16GB -R "span[hosts=1] rusage[mem=16GB]" \
    -o "${indir}/softclip/${nam}.softclip_read_info_extract.log" -e "${indir}/softclip/${nam}.softclip_read_info_extract.err" \
    "${my_stools} view ${indir}/softclip/${f} | cut -f 1,3,4,6 | sort | uniq > ${indir}/softclip/${nam}_reads_info.txt"
done

for f in $(ls -1 *softclipped.bam); do
nam=$(echo $f | sed 's/\_softclipped\.bam//g')

bsub -q medium -app large -w 'ended("softclip_read_extract*")' -J "softclip_bam_extract_${f}" -n 4,4 -M 16GB -R "span[hosts=1] rusage[mem=16GB]" \
    -o "${indir}/softclip/${nam}.softclip_bam_extract.log" -e "${indir}/softclip/${nam}.softclip_bam_extract.err" \
    "python $my_extract -b ${indir}/softclip/${f} -n ${indir}/softclip/${nam}_softclipped_mapped_filtered_reads.txt -o ${indir}/softclip/${nam}_viral_softclip.bam"
done

# Merging the two while still in the softclip directory
cd ${indir}/softclip/

for f in $(ls -1 *softclipped.fq); do

base=$(echo $f | sed 's/_softclipped.fq//g')

bsub -q medium -app large -w 'ended("*_extract_*")' -J "merge_ms_${f}" -n 4,4 -M 16GB -R "span[hosts=1] rusage[mem=16GB]" \
    -o "${indir}/${base}.merge.ms.log" -e "${indir}/${base}.merge.ms.err" \
    "$my_stools merge -f ${indir}/${base}_merged.bam ${indir}/mate/${base}_viral_mate.bam ${indir}/softclip/${base}_viral_softclip.bam
    $my_stools sort ${indir}/${base}_merged.bam -o ${indir}/${base}_sorted.bam"

bsub -q express -app large -w 'ended("merge_*")' -J "host_read_info_extract_${f}" -n 4,4 -M 16GB -R "span[hosts=1] rusage[mem=16GB]" \
    -o "${indir}/${nam}.host_read_info_extract.log" -e "${indir}/${nam}.host_read_info_extract.err" \
    "${my_stools} view ${indir}/${base}_sorted.bam | cut -f 1,3,4,6 | sort | uniq > ${indir}/${base}_reads_info.txt"
done



# Exit
printf -- '\n';
exit 0;
