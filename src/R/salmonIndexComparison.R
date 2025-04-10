#' ---
#' title: "Salmon index comparison"
#' author: "Nicolas Delhomme"
#' date: "`r Sys.Date()`"
#' output:
#'  html_document:
#'    toc: true
#'    number_sections: true
#'    code_folding: hide
#' ---
#' # Setup
#' * Libraries
suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(here)
  library(tidyjson)
  library(tidyr)
  library(tximport)
})

#' * Functions
"compare_sets" <- function(set1,set2,
                           nam1="set1",
                           nam2="set2",
                           export=FALSE,
                           export.dir=here("data/validation"),
                           export.prefix="comparison_"){
  
  # sanity
  stopifnot(is.matrix(set1))
  stopifnot(is.matrix(set2))
  stopifnot(colnames(set1)==colnames(set2))
  
  # overlap
  cg <- intersect(rownames(set1),rownames(set2))

  # create the object
  tb <- bind_cols(set1[cg,] %>% as.data.frame %>% 
                    pivot_longer(ends_with("gz"),
                                 names_to="SampleID",
                                 values_to="Set1"),
                  set2[cg,] %>% as.data.frame %>% 
                    pivot_longer(ends_with("gz"),
                                 values_to="Set2") %>% 
                    select(Set2))
  
  # plot
  p <- ggplot(tb,aes(x=log1p(Set1),y=log1p(Set2))) + 
    # geom_density_2d_filled() + 
    geom_smooth() +
    scale_x_continuous(paste("log TPM",nam1)) +
    scale_y_continuous(paste("log TPM",nam2)) +
    geom_abline(slope=1,intercept=0,colour="red",linetype=2) +
    ggtitle(paste(nam1, "vs.", nam2)) +
    facet_wrap(~SampleID)
  
  print(p)
  
  # export
  if ( isTRUE(export) ) {
    ggsave(file.path(export.dir,paste0(export.prefix,nam1,"_",nam2,".pdf")),p)  
  }
  
  return(NULL)
  
}

#' # Data
reports <- list.files(here("data/RNASeq/validation"),pattern="meta_info.json",recursive=TRUE,full.names=TRUE)
quants <- list.files(here("data/RNASeq/validation"),pattern="quant.sf",recursive=TRUE,full.names=TRUE)

#' ## Percent mapping
pct <- lapply(reports,function(f){
  read_json(f) %>% 
    spread_all() %>% 
    select(percent_mapped) %>% 
    as_tibble()
  }) %>% purrr::list_rbind() %>% 
  cbind(tibble(Set=basename(dirname(dirname(dirname(reports)))),
                 Sample=basename(dirname(dirname(reports)))))

#' ## Quantification
tpms <- lapply(split(quants,basename(dirname(dirname(quants)))),
               function(l){
                 tx <- tximport(l,type="salmon",txOut=TRUE,
                          countsFromAbundance="scaledTPM")
                 cts <- tx$counts
                 colnames(cts) <- basename(dirname(l))
                 return(cts)
               })
               
#' # Visualization
#' 
#' ## Percent mapping
ggplot(pct,aes(x=Set,y=percent_mapped,fill=Set)) + 
  geom_violin() +
  ggtitle("Percentage of mapped reads by index") +
  scale_x_discrete(sprintf("Set (n=%d)",nrow(pct)/length(unique(pct$Set))))

pct <- pct %>% filter(Set!="Potrx2")

ggplot(pct,aes(x=Set,y=percent_mapped,fill=Set)) + 
  geom_violin() +
  ggtitle("Percentage of mapped reads by index") +
  scale_x_discrete(sprintf("Set (n=%d)",nrow(pct)/length(unique(pct$Set))))

#' ## Scatterplot

# This would need ortholog information

#' ### Potra2 vs. Potrx1
#dev.null <- compare_sets(tpms$Potra2,tpms$Potrx1,"Potra2","Potrx1")

#' ### Potra2 vs. Potrs1
#dev.null <- compare_sets(tpms$amRNA,tpms$damRNA,"without_decoy","with_decoy")

#' ### Potrx1 vs. Potrs1
#dev.null <- compare_sets(tpms$amRNA,tpms$lmRNA,"all_mRNA_plus_decoy","longest_mRNA_plus_decoy")

#' # Session Info
#' ```{r session info, echo=FALSE}
#' sessionInfo()
#' ```
