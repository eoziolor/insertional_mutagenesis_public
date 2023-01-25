# Viral integration processing pipeline

This pipeline is intended to take fastq files from either S-EPTS, whole genome (WG), or target-enrichment sequencing (TES) and test for integration between a host and viral vector. Most of the workflow has been validated with TES as it is the preferred method with highest sensitivity and lowest wasted sequence.

### Input needs

The pipeline is broken down into a few steps. They are all intended to work within one distinctly named directory, which will contain all files. The general structure of the directory is as follows.

* viral_genome
+ fastas - directory that holds individual fastas for any elements of the viral genome you'd like to include; can include contaminating regions (repcap/helper/backbone etc.)
    + *Note: fasta files must end in the suffix '.fa'*
    + *All the fasta files will be contanenated into a single file named \<project\>_viral_genome.fa in the project directory. The resulting genome is indexed by _bwa index_ and _samtools faidx_*
* meta
  + contains metadata csv file with specified headers (list)
* raw_fastq
  + contains raw fastq files that need to be named by __sample__ name in metadata and finish in _1 or _1/_2 depending on wether they are single or paired end
  
### Created directories

* trimalign - directory created by running the pipeline
  + *.bam - alignment to host genome
  + *_merged.bam - merge of evidence of IS from softclip and mate directories
  + *_sorted.bam - re-sort based on name for parsing
  + *_reads_info.txt - file containing read name, host chromosome, chromosome location, cigar string
  + mate:
    + *unaligned.fq - reads that did not align properly to host
    + mapped_mates*.bam - the properly aligned mates to unaligned reads
    + *sorted_mates.bam - re-sorted mapped_mates by name (for search purposes)
    + *unaligned_mapped_filtered_reads.txt - reads that didn't align to host, but aligned to viral genome
+ *viral_mate.bam - the mates that aligned to host of reads that didn't align to host, but did align to virus
  + softclip:
    + *softclipped.bam/fq - reads that aligned properly to host (min 30M), but also had a large softclipped portion (min 30M)
    + *softclipped_mapped_filtered.bam - reads from above files that aligned properly (min 30M) to viral genome
    + *softclipped_mapped_filtered_reads.txt - read names for those reads
    + *viral_softclip.bam - host locations of the reads that were softclipped in host and aligned to viral genome
  + full_viral:
    + *.bam - alignment of original raw fastq files to viral genome (to understand viral load overall)
+ secondary - directory for read analysis post extraction
  + reads_info
    + contains reads info files for host/softclip/mate read origin locations - used to merge and get a file that contains host-virus relationships for each insertional site
  + vir_orig
    + contains bam files for viral origins from softclip and mate reads directories
  + is
    + contains *_sorted.bam files from trimalign dir

### Workflow

__Pre-processing: Compile & Index Viral Genome Sequences__

Extract the fasta files for each viral element against which TES was created and merge them in a multi-fasta for viral genome. Then index with bwa and samtools.

_Concatenating viral elements_ - make sure your elements are in individual fastas in the _fastas_ directory below

```{bash}
# Directories
parent=<path-to-base-dir>
project=<name-of-specific-project>

# Software
my_bwa=/hpc/grid/dsrd-invtox/workspace/program/bwa-0.7.17/bwa
my_stools=/hpc/grid/dsrd-invtox/workspace/program/samtools-1.9/bin/samtools

# Concatenating fastas
cat ${parent}/${project}/viral_genome/fastas/*.fa >> ${parent}/${project}/${project}_viral_genome.fa

# Indexing that genome
$my_bwa index ${parent}/${project}/${project}_viral_genome.fa
$my_stools faidx ${parent}/${project}/${project}_viral_genome.fa

```

__Step1: Identification of insertional sites (IS)__

This step uses local aligner (bwa mem) as well as duplication removal (samtools markdup) to align reads to host genome. After alignment it uses samtools to extract two sets of reads 1) reads that _did not_ align to the host genome, but whose mates _did_ properly align. Both the forementioned reads (unaligned and aligned mates) go in a sub-directory of mates. 2) reads that only partially aligned to the host and had at minimum 30 bp that are _softclipped_ for lack of alignment. These go in a _softclip_ subdirectory.

_Alignment to host_ - make sure to place your fastq files following directions above

```{bash}
# Creating fastq directory
mkdir -p ${parent}/${project}/raw_fastq

# Running alignment
my_step1=/hpc/grid/dsrd-invtox/workspace/scripts/v.int/step1.trimalign.sh

${my_step1} \
[-g cyno|mouse|human] enter as a string the option of which host genome you want to use \
[-i path-to-data-dir] path to fastq files (suggested: ${parent}/${project}/raw_fastq) \
[-e fastq|fq|fastq.gz|fq.gz] enter as string the suffix of your fastq files \
[-r single|paired] as string the type of reads that you have \
[-o path-to-output-dir] output directory (suggested: ${parent}/${project}/trimalign) \
[-w] genewerk data flag; parameters slightly tweaked to handle genewerk/s-epts data \
[-f] full wg or tse pipeline flag; this will enable automatically pushing step2 after step1 is finished \
[-h] help menu

```

__*Step2: Extraction of putative IS (reads that do not align properly to host)*__  

*Note: This next part is usually automated after the above if __-f__ is selected and does not need to be run*

```{bash}
my_step2=/hpc/grid/dsrd-invtox/workspace/scripts/v.int/step2.ms.extract.sh

${my_step2} \
[-i path-to-data-dir] the directory where alignments are (suggested: ${parent}/${project}/trimalign) \
[-g cyno|mouse|human] enter as a string the option of which host genome you want to use \
[-h] help menu

```

__Step3: Alignment to viral genome - alignment of putative IS to viral genome__  

```{bash}
my_step3=/hpc/grid/dsrd-invtox/workspace/scripts/v.int/step3.viralign.sh 

${my_step3} \
[-i path-to-orig-fastq] path to fastq files (suggested: ${parent}/${project}/raw_fastq) \
[-e fastq|fq|fastq.gz|fq.gz] enter as string the suffix of your fastq files \
[-r single|paired] as string the type of reads that you have \
[-m path-to-folder-with-mate-and-softclip] step1 output directory (suggested: ${parent}/${project}/trimalign) \
[-v path-to-viral-genome] viral genome fasta (suggested: ${parent}/${project}/${project}_viral_genome.fa) \
[-w] genewerk data flag; parameters slightly tweaked to handle genewerk/s-epts data \
[-h] help menu

```

__Step4: Creating final IS output - takes the final outputs and extracts filtered IS reads and information__  

```{bash}
my_step4=/hpc/grid/dsrd-invtox/workspace/scripts/v.int/step4.count.sh 

${my_step4} \
[-g cyno|mouse|human] enter as a string the option of which host genome you want to use \
[-i path-to-folder-with-mate-and-softclip] step1 output directory (suggested: ${parent}/${project}/trimalign) \
[-w] genewerk data flag; parameters slightly tweaked to handle genewerk/s-epts data \
[-h] help menu

```

## Secondary analysis

_Moving IS file info_ - putting IS info into sub-folders for secondary analysis

```{bash}
# Moving to sub-directories
indir=${parent}/${project}/trimalign
outdir=${parent}/${project}/secondary

## Sorted IS files to secondary analysis
mkdir -p ${outdir}/is
cp ${indir}/*sorted.bam ${outdir}/is/

## Moving viral info reads
mkdir -p ${outdir}/reads_info
cp ${indir}/*reads_info.txt ${outdir}/reads_info
cp ${indir}/mate/*reads_info.txt ${outdir}/reads_info
cp ${indir}/softclip/*reads_info.txt ${outdir}/reads_info

## Vir origins
mkdir -p ${outdir}/vir_orig/
cp ${indir}/mate/*unaligned_mapped_filtered.bam ${outdir}/vir_orig/
cp ${indir}/softclip/*softclipped_mapped_filtered.bam ${outdir}/vir_orig/

```

_Generating IS info data frames_ - using metadata csv to load and reformat data into secondary analysis format

```{r}
# Functions
source("/hpc/grid/dsrd-invtox/workspace/oziole/CompTox_insertmuta_DNA/functions/cigarParse.R")
source("/hpc/grid/dsrd-invtox/workspace/oziole/CompTox_insertmuta_DNA/functions/loadIS.R")

# Folders and files
parent <- <path-to-base-dir>
project <- <name-of-specific-project>
meta.data <-read_csv(paste0(parent,"/",project,"/meta/",project,"_meta.csv"))

# Generating secondary files
loadIS(meta = meta.data,
       base.dir = paste0(parent,"/",project,"/secondary/"), # directory of secondary analysis
       pattern.is = "sorted.bam$", # suffix pattern of is bam files
       set = project, # any string name you want to give to this project
       do.is = TRUE,
       do.vir.origin = TRUE, 
       do.read.info = TRUE)
```
