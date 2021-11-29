args <- commandArgs(TRUE)

defaultW <- getOption("warn")
options(warn = -1)

suppressPackageStartupMessages({
  library(tidyverse)
  library(rtracklayer)
})

options(warn = defaultW)

gtf_file <- args[1]
input_genes <- args[2]
prefix <- args[3]
outpath <- args[4]

ids <- str_split(input_genes, ";") %>% unlist()

gtf <- rtracklayer::import(gtf_file)
gtf = keepStandardChromosomes(gtf, species = "Homo_sapiens", pruning.mode="coarse")
gtf.df=as.data.frame(gtf,stringsAsFactor=F)


x <- gtf.df %>%
  dplyr::filter(type %in% "gene", gene_id %in% ids) %>%
  dplyr::select(seqnames, start, end, gene_id)

write.table(x, str_c(outpath, "/", prefix, "_geneLocus.bed"), row.names=FALSE, col.names=FALSE, sep="\t", quote=FALSE)
