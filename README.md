# Poplar CRISPR whole genome sequencing

## Setup

```bash
# data
mkdir data
ln -s /mnt/ada/projects/aspseq/mschmid/poplar-CRISPR-WGS/RNA-seq/raw data/RNASeq
ln -s /mnt/ada/projects/aspseq/mschmid/poplar-CRISPR-WGS/raw data/WGS

# tools
ln -s /mnt/picea/storage/singularity .

# reference
mkdir reference
ln -s /mnt/picea/storage/reference/rRNA/sortmerna/v4.3.4 reference/sortmerna
ln -s /mnt/picea/storage/reference/Populus-tremula/v2.2/indices/salmon1.6.0/ reference/salmon
ln -s /mnt/picea/storage/reference/Populus-tremula/v2.2/annotation reference/annotation
ln -s /mnt/picea/storage/reference/Illumina/adapters/TruSeq3-PE-2.fa reference/trimmomatic

# results
ln -s /mnt/ada/projects/aspseq/mschmid/poplar-CRISPR-WGS/RNA-seq/preprocessed analysis
```

## Description

There are two CRISPR lines with different targets:

* For the Line 26, the target locus is `Potri.006G174000/Potra001877g14982`
* For the Line 03 it is `Potri.018G096200/Potra001273g10998`

## Tasks

### RNA-Seq

Run the BioQA and the DE to compare line 3 and 26 to T89 and versus each-other.

### WGS

The aim is to look for off-target effect of CRISPR in other genes of the same family. Primarily from the genes in the SM and LSM families.

