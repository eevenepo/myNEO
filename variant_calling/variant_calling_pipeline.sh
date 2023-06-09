#!/bin/bash

####################################################
# Cell line variant calling (illustrated on CCK81) # 
####################################################

# NOTE: Significant intermediary data is generated. To avoid cluttering backup space, we recommend running this script in non-backed up directories

### Example data
cell_line="CCK81"
n_threads=16

read1_url="ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR339/ERR3396992/Sample_CCK81_T3_R1_fastq.gz"
read2_url="ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR339/ERR3396992/Sample_CCK81_T3_R2_fastq.gz"

ref_genome_url="ftp://ftp.ensembl.org/pub/release-94/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"
gtf_file_url="https://ftp.ensembl.org/pub/release-94/gtf/homo_sapiens/Homo_sapiens.GRCh38.94.gtf.gz"


### Directory structure

mkdir -p $cell_line \
         $cell_line/raw_data/ \
         $cell_line/raw_data/$cell_line/ \
         $cell_line/raw_data/hg38/ \
         $cell_line/quality_control/ \
         $cell_line/alignment/ \
         $cell_line/variant_calls/

cd $cell_line


### Data

# Download compressed cell line sequencing data  
read1="raw_data/"$cell_line"/read_1.fq"
read2="raw_data/"$cell_line"/read_2.fq"

wget -O $read1".gz" $read1_url
wget -O $read2".gz" $read2_url

gunzip $read1".gz"
gunzip $read2".gz"

# Download hg38 reference genome from Ensembl
ref_genome=raw_data/hg38/"primary_assembly.fa"
gtf_file=raw_data/hg38/"annotation_info.gtf"

wget -O $ref_genome".gz" $ref_genome_url
wget -O $gtf_file".gz" $gtf_file_url

gunzip $ref_genome".gz"
gunzip $gtf_file".gz"


### Quality check
fastqc $read1 $read2 \
    -o alignment \
    -t $n_threads
  
  
### Alignment
cd alignment

# Index reference genome
STAR --runMode genomeGenerate \
     --genomeFastaFiles $ref_genome \
     --sjdbGTFfile $gtf_file
     --genomeDir alignment \
     --runThreadN $n_threads
     
# Align reads to reference

STAR --runMode alignReads \
     --genomeDir alignment \
     --readFilesIn $read1 $read2 \
     --outSAMtype BAM SortedByCoordinate \
     --outFileNamePrefix $cell_line  \
     --runThreadN $n_threads

output_bam=$cell_line".sortedByCoord.out.bam"

### Variant calling
cd ../variant_calls

gatk Mutect2 \
    -R $ref_genome \
    -I $output_bam \
    -O $cell_line_"variants.vcf"

### Free space 
rm -rf $cell_line/raw_data

