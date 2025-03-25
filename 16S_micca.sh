#!/bin/bash

#PBS -N micca_16S
#PBS -l select=1:ncpus=20:mem=20G
#PBS -q workq
#PBS -l walltime=120:00:00
#PBS -M paola.maragno98@gmail.com
#PBS -m e

source /data/apps/anaconda3/bin/activate /data/bioinfo/FC/16S_analysis

/data/apps/singularity/3.11.3/bin/singularity exec --bind /data/bioinfo/FC/:/data/bioinfo/FC/ /data/bioinfo/FC/singularity_cache/micca_1.7.2.sif /bin/bash /data/bioinfo/FC/pmaragno/16S/micca_pipeline.sh

cd /data/bioinfo/FC/pmaragno/16S/fastq

gunzip *_R*.gz

/usr/local/bin/micca mergepairs -i *_R1*.fastq -o merged.fastq -l 32 -d 8

/usr/local/bin/micca stats -i merged.fastq -o ./stats

/usr/local/bin/micca trim -i merged.fastq -o trimmed.fastq -w CCTACGGGNGGCWGCAG -r GACTACHVGGGTATCTAATCC -W -R -c

/usr/local/bin/micca filterstats -i trimmed.fastq -o ./filterstats

/usr/local/bin/micca filter -i trimmed.fastq -o filtered.fasta -e 1.0 -m 400 -t

/usr/local/bin/micca otu -i filtered.fasta -o ./denovo_unoise_sv -m denovo_unoise -c -t 18

/usr/local/bin/micca classify -m rdp -i ./denovo_unoise_sv/otus.fasta -o ./denovo_unoise_sv/taxa.txt

wget ftp://ftp.fmach.it/metagenomics/micca/dbs/core_set.tar.gz

tar -zxvf core_set.tar.gz

/usr/local/bin/micca msa -m nast -i ./denovo_unoise_sv/otus.fasta -o ./denovo_unoise_sv/msa.fasta --nast-template core_set_aligned.fasta.imputed  --nast-threads 18

/usr/local/bin/micca tree -i ./denovo_unoise_sv/msa.fasta -o ./denovo_unoise_sv/tree.tree

/usr/local/bin/micca root -i ./denovo_unoise_sv/tree.tree -o ./denovo_unoise_sv/tree_rooted.tree

/usr/local/bin/micca tablestats -i ./denovo_unoise_sv/otutable.txt -o ./denovo_unoise_sv/tablestats

/usr/local/bin/micca tablerare -i ./denovo_unoise_sv/otutable.txt -o ./denovo_unoise_sv/otutable_rare.txt -d 10455

/usr/local/bin/micca tabletotax -i ./denovo_unoise_sv/otutable_rare.txt -t ./denovo_unoise_sv/taxa.txt -o ./denovo_unoise_sv/taxtables

