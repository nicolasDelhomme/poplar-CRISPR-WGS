# The splicing genes SmEa and SmEb regulate plant development during vegetative growth in poplar

__Goretti D.1,2&, Collani S.1,2, Marcon A.3, Nilsson O.3, Schmid M.1,4&__

1. Umeå Plant Science Centre, Department of Plant Physiology, Umeå University, SE-90187 Umeå, Sweden
2. DISSTE, University of Eastern Piedmont, Vercelli, 13100, Italy
3. Umeå Plant Science Centre, Department of Forest Genetics and Plant Physiology, Swedish University of Agricultural Sciences, SE-90183 Umeå, Sweden
4. Department of Plant Biology, Linnean Center for Plant Biology, Swedish University of Agricultural Sciences, S-75007 Uppsala, Sweden

& corresponding author

## Abstract

Spliceosomes, the molecular machinery that catalyzes RNA splicing, contain at their core heptameric rings of Sm (or LSm) proteins. Loss of Sm gene function can have detrimental consequences on plant growth and development, as exemplified by the strong temperature-sensitive phenotype of the sme1/pcp mutant in the model plant species Arabidopsis thaliana. In this work, we characterized the SmE genes from poplar. We found that poplar contains two paralogs SmE genes, PtSmEa and PtSmEb, that encode for identical proteins. Mutagenesis using CRISPR/Cas9 revealed both overlapping and specific functions of the two poplar SmE genes, which was further substantiated by transcriptome and alternative splicing analyses using RNA-seq. Our study provides a first glimpse at the function of Sm genes in a woody species and suggests a key role of PtSmEs in regulating vegetative growth in trees.

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

* For the Line 26, the target locus is `Potri.006G174000/Potra001877g14982` (Potra2n6c13821)
* For the Line 03 it is `Potri.018G096200/Potra001273g10998` ( 	Potra2n6c13821)

## Tasks

### RNA-Seq

Run the BioQA and the DE to compare line 3 and 26 to T89 and versus each-other.

### WGS

The aim is to look for off-target effect of CRISPR in other genes of the same family. Primarily from the genes in the SM and LSM families.

The guide RNAs are:

Line3_CRISPR_4108_gRNA1: GACTAGACGTACAATGGGTT
Line3_CRISPR_4108_gRNA2: GCACATAAGCAGATACGCTC

Line26_CRISPR_4110_gRNA1: GGTCATAATACGCTGGACTT
Line26_CRISPR_4110_gRNA2: TATAAAGAGCAAGAATTGAC

1. Ask Daniella on which genome she designed the guides. (Nico)
2. Align Potra2n6c13821 to the new assembly, figure out how many loci are present. Blastn. (Kristina)
3. WGS, look for indels / SNPs, in particular for the loci identified in the DE report (1.2.2 - other genes of interest). BWA -> GATK. (Kristina)
4. Align the gRNA to the assembly (blastn) and intersect with 3. To identify and ignore effects due to the transformation or other artefacts. If there are some relevant cases, we could check if there are within a gene. (Kristina)
5. align RNA-Seq with STAR and look at the loci targeted by the gRNA. (Fai)

SNPs / Indels (Yes) -> Overlap gRNA (Yes) -> Is it within / next to a gene (Yes) => putative off-target.

For SNPs, extend to a 100bp total and check for overlap with putative gRNA target sites.

