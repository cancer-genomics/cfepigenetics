library(stringr)
library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg19)
library(tidyverse)
library(reshape2)

## Get SGE_ID for parallelizing this script over the 9 files

args <- commandArgs(trailingOnly = TRUE)
case <- as.numeric(args[1])

## Data locations

path <- "/dcl01/scharpf1/data/cristiano/projects/bedhighcov/granges"
outpath <- "../../../cfepigenetics_data/end_motifs_sheared"
files <- list.files(path)

if (!dir.exists(outpath)) {
  dir.create(outpath)
}

## Process files

file <- files[case]
name <- unlist(str_split(file, "[.]"))[1]
gr <- readRDS(file.path(path, file))
gr <- gr[,c()]
seqinfo(gr) <- seqinfo(Hsapiens)

# Extended fragments that now fall out of the ranges of the chromosomes, are being eliminated

gr_start <- gr
start(gr_start) <- start(gr_start) - 1
end(gr_start) <- start(gr_start) + 2
idx <- GenomicRanges:::get_out_of_bound_index(gr_start)
if (length(idx) != 0L) {
  gr_start <- gr_start[-idx]
}
strand(gr_start) <- "+"
tb_start <- tibble()
for (subdiv in c(1:100)) {
  print(subdiv)
  gr_start_temp <- gr_start[c((round((length(gr_start) / 100) * (subdiv-1))+1):round((length(gr_start) / 100) * subdiv))]
  tb_start_temp <- as_tibble(as.data.frame(table(as.character(getSeq(Hsapiens, gr_start_temp)))))
  tb_start <- rbind(tb_start, tb_start_temp)
  tb_start <- tb_start %>%
    group_by(Var1) %>%
    summarise(Freq = sum(Freq))
}
rm(gr_start)

gr_end <- gr
rm(gr)
start(gr_end) <- end(gr_end) - 1
end(gr_end) <- start(gr_end) + 2
idx <- GenomicRanges:::get_out_of_bound_index(gr_end)
if (length(idx) != 0L) {
  gr_end <- gr_end[-idx]
}
strand(gr_end) <- "-"
tb_end <- tibble()
for (subdiv in c(1:100)) {
  print(subdiv)
  gr_end_temp <- gr_end[c((round((length(gr_end) / 100) * (subdiv-1))+1):round((length(gr_end) / 100) * subdiv))]
  tb_end_temp <- as_tibble(as.data.frame(table(as.character(getSeq(Hsapiens, gr_end_temp)))))
  tb_end <- rbind(tb_end, tb_end_temp)
  tb_end <- tb_end %>%
    group_by(Var1) %>%
    summarise(Freq = sum(Freq))
}
rm(gr_end)

tb <- rbind(tb_start, tb_end)
tb <- tb %>%
  group_by(Var1) %>%
  summarise(Freq = sum(Freq))

saveRDS(tb, file.path(outpath, paste0(name,".rds")))


