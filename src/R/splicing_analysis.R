library(ASpli)
library(GenomicFeatures)
library(readxl)
library(tidyverse)


#' GFF preprocessing
gtfFileName <- "reference/gff/Potra02_genes.gff"
genomeTxDb <- makeTxDbFromGFF( gtfFileName )

saveDb(genomeTxDb,file="Potra_genes.sqlite")

#'  Feature extraction
features <- binGenome( genomeTxDb )

#' BAMS and target file
bams <- list.files("/mnt/ada/projects/aspseq/mschmid/poplar-CRISPR-WGS/RNA-seq/preprocessed/STAR", pattern="*.bam$", full.names = TRUE)

targets <- data.frame(row.names = paste0(rep(c("L03", "L26", "T89"), each = 4), rep("_", times = 12), rep(1:4, times = 3)),
                      bam = bams,
                      genotype = c( 'line3','line3','line3','line3','line26','line26','line26','line26','control','control','control','control'),
                      stringsAsFactors = FALSE)

mBAMs <- data.frame( bam = list.files("/mnt/ada/projects/aspseq/mschmid/poplar-CRISPR-WGS/RNA-seq/preprocessed/STAR/merged", pattern="*.bam$", full.names = TRUE),
                     condition = c("line3","line26", "control"))

gbcounts <- gbCounts(features=features, targets=targets,
                     minReadLength = 100, maxISize = 50000, libType = "PE", strandMode = 2)
gbcounts

asd <- jCounts(counts=gbcounts, features=features, minReadLength=100)

asd

#' Differential splicing expression
gb_line3 <- gbDUreport(gbcounts, contrast = c(1,0,-1))
gb_line26 <- gbDUreport(gbcounts, contrast = c(0,1,-1))

#' Differential junction usage analysis
jdur_line3 <- jDUreport(asd, contrast=c(1,0,-1))
jdur_line26 <- jDUreport(asd, contrast=c(0,1,-1))

#' Bin and junction signal integration
sr_line3 <- splicingReport(gb_line3, jdur_line3, counts=gbcounts)
sr_line26 <- splicingReport(gb_line26, jdur_line26, counts=gbcounts)
sr <- splicingReport(gb_t89, jdur_t89, counts=gbcounts)

#' Summary of integration of splicing signals along genomic-regions
is_line3 <- integrateSignals(sr_line3,asd)
is_line26 <- integrateSignals(sr_line26,asd)

#' Export results
exportIntegratedSignals(is_line3,sr=sr_line3,
                        output.dir = "Line3_AS_stranded",
                        counts=gbcounts,features=features,asd=asd,
                        mergedBams = mBAMs)

exportIntegratedSignals(is_line26,sr=sr_line26,
                        output.dir = "Line26_AS_stranded",
                        counts=gbcounts,features=features,asd=asd,
                        mergedBams = mBAMs)


#' T89

#' # Line 3
#' 
#' #' Read in Excel
as_line3 <- read_excel("doc/line3_control_stranded.xlsx")
#' Get the frequency per splicing event
events_line3 <- as.data.frame(table(as_line3$Event))
#' Remove unassigned
events_line3 <- events_line3[-1,]
#' Get total number of splicing events
total_events_line3 <- sum(events_line3$Freq)
#' Calculate percentages
events_line3$Percentage <- (events_line3$Freq / total_events_line3) * 100
#' Combine annotations
events_line3$Annotation <- c("Alt 3'", "Alt 3'", "Alt 5'", "Alt 5'", "ASCE", "ES", "ES", "IoR", "IR", "IR", "Novel", "Novel", "Novel", "Novel", "Undefined")
#' Add annotation for line for plotting
events_line3$Line <- rep("Line3",15)
#' Filter some types of events
events_line3 <- events_line3[which(!events_line3$Annotation %in% c("ASCE", "Undefined", "IoR", "Novel")),]

#' # Line 26
#' 
#' #' Read in Excel
as_line26 <- read_excel("doc/line26_control_stranded.xlsx")
#' Get the frequency per splicing event
events_line26 <- as.data.frame(table(as_line26$Event))
#' Remove unassigned
events_line26 <- events_line26[-1,]
#' Get total number of splicing events
total_events_line26 <- sum(events_line26$Freq)
#' Calculate percentages
events_line26$Percentage <- (events_line26$Freq / total_events_line26) * 100
#' Combine annotations
events_line26$Annotation <- c("Alt 5'/3'","Alt 3'", "Alt 3'", "Alt 5'", "Alt 5'", "ASCE", "ES", "ES", "IoR", "IR", "IR", "Novel", "Novel", "Novel", "Novel", "Undefined")
#' Add annotation for line for plotting
events_line26$Line <- rep("Line26",16)
#' Filter some types of events
events_line26 <- events_line26[which(!events_line26$Annotation %in% c("ASCE", "Undefined", "IoR", "Novel")),]

#' Local AS events annotated in the Potra genome (from SUPPA2)
#' 

es_anno <- 2874
ir_anno <- 7215
three_anno <- 4128
five_anno <- 5080

total_as_annotated <- es_anno + ir_anno + three_anno + five_anno

events_as_anno <- data.frame(
  Var1 = c("ES","IR", "Alt 3'", "Alt 5'"),
  Freq = c(es_anno, ir_anno, three_anno, five_anno),
  Annotation = c("ES", "IR", "Alt 3'", "Alt 5'"),
  Percentage = c(es_anno, ir_anno, three_anno, five_anno) / total_as_annotated, 
  Line = rep("Total Annotated", 4)
)


#' Combine events for Line 3 and Total Annotated
events_line3_anno <- rbind(events_line3,events_as_anno)

#' Plot 
ggplot(data = events_line3_anno, aes(x = Line, y = Percentage, fill = Annotation,
                              label = sprintf("%.02f", Percentage))) + geom_bar(position = "fill", stat = "identity") +
  scale_fill_viridis_d() +  # Colorblind-friendly and print-friendly palette
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 14)) +
  labs(x = "", y = "% AS events")

#' Combine events for Line 3 and Total Annotated
events_line26_anno <- rbind(events_line26,events_as_anno)

#' Plot 
ggplot(data = events_line26_anno, aes(x = Line, y = Percentage, fill = Annotation,
                                     label = sprintf("%.02f", Percentage))) + geom_bar(position = "fill", stat = "identity") +
  scale_fill_viridis_d() +  # Colorblind-friendly and print-friendly palette
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 14)) +
  labs(x = "", y = "% AS events")


#' Compare with DE genes
#' 
#' Line 3

DEG_line3 <- read_csv("analysis/DE/Line03_genes.csv") 

colnames(DEG_line3)[1] <- "Locus"

common_genes_line3 <- intersect(as_line3$Locus, DEG_line3$Locus)

print(paste0("Number of Transcription and AS regulated genes: ", length(common_genes_line3)))

#' 
#' Line 26

DEG_line26 <- read_csv("analysis/DE/Line26_genes.csv") 

colnames(DEG_line26)[1] <- "Locus"

common_genes_line26 <- intersect(as_line26$Locus, DEG_line26$Locus)

print(paste0("Number of Transcription and AS regulated genes: ", length(common_genes_line26)))
