#at the project root directory
mkdir data
mkdir reference
mkdir reference/trimmomatic
mkdir analysis

#soft links
ln -s /mnt/ada/projects/aspseq/mschmid/poplar-CRISPR-WGS/RNA-seq/raw data/RNASeq
ln -s /mnt/ada/projects/aspseq/mschmid/poplar-CRISPR-WGS/raw data/WGS
ln -s /mnt/picea/storage/singularity .
ln -s /mnt/picea/storage/reference/rRNA/sortmerna/v4.3.4 reference/sortmerna
ln -s /mnt/picea/storage/reference/Populus-tremula/v2.2/indices/salmon1.6.0/ reference/salmon
ln -s /mnt/picea/storage/reference/Illumina/adapters/TruSeq3-PE-2.fa reference/trimmomatic/
ln -s /mnt/ada/projects/aspseq/mschmid/poplar-CRISPR-WGS/RNA-seq/preprocessed analysis/