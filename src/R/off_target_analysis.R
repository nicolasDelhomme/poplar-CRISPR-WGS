
suppressPackageStartupMessages({
  library(data.table)
  library(splitstackshape)
  library(tidyverse)
  library(VariantAnnotation)
  library(RColorBrewer)
  library(dplyr)
  library(ggplot2)
  library(R.utils)
})


extract_DP <- function(column) {
  dp_entries <- str_extract_all(column, "DP=\\d+")
  return(ifelse(dp_entries != "character(0)", dp_entries, ""))
}

#' Function to filter a VCF file from GATK Joint Genotyping workflow. 
#' The input is required to have the following structure: 
#' #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT sample1 sample2 sample3 ...

filterVCF <- function(vcf,file_prefix,DP_lower,DP_upper,sample_names) {
  
  #'  # Load VCF file into a data table
  vcf_file <- fread(vcf, sep="\t", skip = "#CHROM")
  
  #'  # Split the columns with multiple information fields into sepearate columns
  snps_split <- cSplit(vcf_file, 8, sep = ";", "wide")
  
  #'  # Apply the function to INFO_4 and INFO_5 columns (either of which contain the DP information) and create a new column DP
  
  snps_split <- snps_split %>%
    mutate(DP_4 = extract_DP(INFO_04),
           DP_5 = extract_DP(INFO_05)) %>% 
    unite(DP,DP_4,DP_5) %>%
    mutate(DP=str_replace_all(DP,"_","")) 
  
  snps <- as_tibble(snps_split) %>% dplyr::select(., `#CHROM`, POS, REF, ALT, QUAL, INFO_01, INFO_02, INFO_03, DP, any_of(sample_names))
  
  #'  # Rename the new columns
  snps <- snps %>% rename_with(~gsub("INFO_01", "AC", .), ends_with("_01")) %>%
    rename_with(~gsub("INFO_02", "AF", .), ends_with("_02")) %>%
    rename_with(~gsub("INFO_03", "AN", .), ends_with("_03")) %>%
    mutate(AC = str_replace(AC, "AC=", "")) %>%
    mutate(AF = str_replace(AF, "AF=", "")) %>%
    mutate(AN = str_replace(AN, "AN=", "")) %>%
    mutate(DP = str_replace(DP, "DP=", ""))
  
  snps$DP <- as.numeric(snps$DP)
  snps$AF <- as.numeric(snps$AF)
  
  snps_filtered <- snps %>% filter(DP >= DP_lower) %>% filter(DP <= DP_upper)
  snps_filtered_homo <- snps %>% filter(DP >= DP_lower) %>% filter(DP <= DP_upper) %>% filter(AF == 1.0)
  
  #'  # Write output to tsv files
  write_tsv(snps,paste0(file_prefix, "_all_SNPs.tsv"))
  write_tsv(snps_filtered,paste0(file_prefix, "_filtered_SNPs.tsv"))
  write_tsv(snps_filtered_homo,paste0(file_prefix, "_filtered_homozygous_SNPs.tsv"))
  
  return(snps) 
}

#' # Load VCF 
vcf <- "data/WGS/gatk/gDNA_line26_EKDN230011267-1A_HVCGCDSX5_HYCW3DSX5_L3_L4_merged.sorted_mkdup.snps.indels.vcf.gz"
system(paste("srun -A u2023008 apptainer exec -B /mnt:/mnt /mnt/picea/storage/singularity/vcftools_v0.1.16.sif vcftools --gzvcf", vcf,"--site-depth 2>/dev/null"))
data <- read.table("out.ldepth", header = TRUE)

#' # Basic variant filtering 
#' 
#' ## Mean read depth - decide on DP threshold
x <- as.data.frame(table(data$SUM_DEPTH))
lower <- 0.75 * median(data$SUM_DEPTH)
upper <- median(data$SUM_DEPTH) + 1 * sd(data$SUM_DEPTH)
xupper <- ceiling(upper/100) * 100
ggplot(x, aes(x = as.numeric(as.character(Var1)), y = Freq)) + geom_line() + xlab("Depth") +
  ylab("bp") + xlim(0, xupper) + geom_vline(xintercept = lower, color = "red",
                                            linewidth = 1.3) + geom_vline(xintercept = upper, color = "red", linewidth = 1.3) + ggtitle("Example threshold: 0.8X median depth, median depth + 2sd")

#' ## Variant quality distribution
system(paste("srun -A u2023008 apptainer exec -B /mnt:/mnt /mnt/picea/storage/singularity/vcftools_v0.1.16.sif vcftools --gzvcf", vcf,"--site-quality 2>/dev/null"))
data <- read.table("out.lqual", header = TRUE)
ggplot(subset(data, QUAL < 1000), aes(x = QUAL)) + geom_histogram(fill = "white",
                                                                  color = "black", bins = 50) + xlab("Quality value") + ylab("Count") + geom_vline(xintercept = 30,
                                                                                                                                                   color = "red", size = 1.3) + ggtitle("Example threshold: Q30")
filterVCF(vcf,"results/gDNA_line26",round(lower,0),round(upper,0),c("gDNA_line3", "gDNA_line26"))


#' # Load VCF 
vcf <- "data/WGS/gatk/gDNA_line3_EKDN230011266-1A_HVGLKDSX5_L1.sorted_mkdup.filtered.snps.indels.vcf.gz"
system(paste("srun -A u2023008 apptainer exec -B /mnt:/mnt /mnt/picea/storage/singularity/vcftools_v0.1.16.sif vcftools --gzvcf", vcf,"--site-depth 2>/dev/null"))
data <- read.table("out.ldepth", header = TRUE)

#' # Basic variant filtering 
#' 
#' ## Mean read depth - decide on DP threshold
x <- as.data.frame(table(data$SUM_DEPTH))
lower <- 0.75 * median(data$SUM_DEPTH)
upper <- median(data$SUM_DEPTH) + 1 * sd(data$SUM_DEPTH)
xupper <- ceiling(upper/100) * 100
ggplot(x, aes(x = as.numeric(as.character(Var1)), y = Freq)) + geom_line() + xlab("Depth") +
  ylab("bp") + xlim(0, xupper) + geom_vline(xintercept = lower, color = "red",
                                            linewidth = 1.3) + geom_vline(xintercept = upper, color = "red", linewidth = 1.3) + ggtitle("Example threshold: 0.8X median depth, median depth + 2sd")

#' ## Variant quality distribution
system(paste("srun -A u2023008 apptainer exec -B /mnt:/mnt /mnt/picea/storage/singularity/vcftools_v0.1.16.sif vcftools --gzvcf", vcf,"--site-quality 2>/dev/null"))
data <- read.table("out.lqual", header = TRUE)
ggplot(subset(data, QUAL < 1000), aes(x = QUAL)) + geom_histogram(fill = "white",
                                                                  color = "black", bins = 50) + xlab("Quality value") + ylab("Count") + geom_vline(xintercept = 30,
                                                                                                                                                   color = "red", size = 1.3) + ggtitle("Example threshold: Q30")
filterVCF(vcf,"results/gDNA_line3",round(lower,0),round(upper,0),c("gDNA_line3", "gDNA_line26"))
