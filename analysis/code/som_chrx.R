library("GenomicRanges")
library("tidyverse")

args = commandArgs(trailingOnly = TRUE)
case = as.integer(args[1])

val_files <- readRDS("/dcl02/leased/cglab/rscharpf/mnoe/PROJECTS/cfdna_methylation/data/selected_validation_files.rds")
luc_files <- readRDS("/dcl02/leased/cglab/rscharpf/mnoe/PROJECTS/cfdna_methylation/data/selected_lucas_files.rds")

val <- val_files[which(val_files$patient_type == "healthy"),]$pgdx_id
luc_h <- luc_files[which(luc_files$patient_type == 2),]$pgdx_id
h <- c(val, luc_h)
c <- luc_files[which(luc_files$patient_type == 3),]$pgdx_id

temp <- readRDS(h[case])

ls <- list.files("/dcs04/scharpf/data/mnoe/motif", full.names = T)
f <- tibble()
for (file in ls) {
  print(file)
  gr <- readRDS(file)
  motif <- sub(".rds", "", basename(file))
  start(gr[which(strand(gr) == "+")]) <- start(gr[which(strand(gr) == "+")]) + 1
  end(gr[which(strand(gr) == "-")]) <- end(gr[which(strand(gr) == "-")]) - 1
  
  st_som_is <- length(unique(queryHits(findOverlaps(temp, gr[which((strand(gr) == "+") & (seqnames(gr) %in% paste0("chr", c(1:22))) & (gr$group == "hg19_cpg_islands"))], type = "start"))))
  en_som_is <- length(unique(queryHits(findOverlaps(temp, gr[which((strand(gr) == "-") & (seqnames(gr) %in% paste0("chr", c(1:22))) & (gr$group == "hg19_cpg_islands"))], type = "end"))))
  st_x_is <- length(unique(queryHits(findOverlaps(temp, gr[which((strand(gr) == "+") & (seqnames(gr) == "chrX") & (gr$group == "hg19_cpg_islands"))], type = "start"))))
  en_x_is <- length(unique(queryHits(findOverlaps(temp, gr[which((strand(gr) == "-") & (seqnames(gr) == "chrX") & (gr$group == "hg19_cpg_islands"))], type = "end"))))
  
  st_som_sho <- length(unique(queryHits(findOverlaps(temp, gr[which((strand(gr) == "+") & (seqnames(gr) %in% paste0("chr", c(1:22))) & (gr$group == "hg19_cpg_shores"))], type = "start"))))
  en_som_sho <- length(unique(queryHits(findOverlaps(temp, gr[which((strand(gr) == "-") & (seqnames(gr) %in% paste0("chr", c(1:22))) & (gr$group == "hg19_cpg_shores"))], type = "end"))))
  st_x_sho <- length(unique(queryHits(findOverlaps(temp, gr[which((strand(gr) == "+") & (seqnames(gr) == "chrX") & (gr$group == "hg19_cpg_shores"))], type = "start"))))
  en_x_sho <- length(unique(queryHits(findOverlaps(temp, gr[which((strand(gr) == "-") & (seqnames(gr) == "chrX") & (gr$group == "hg19_cpg_shores"))], type = "end"))))
  
  st_som_she <- length(unique(queryHits(findOverlaps(temp, gr[which((strand(gr) == "+") & (seqnames(gr) %in% paste0("chr", c(1:22))) & (gr$group == "hg19_cpg_shelves"))], type = "start"))))
  en_som_she <- length(unique(queryHits(findOverlaps(temp, gr[which((strand(gr) == "-") & (seqnames(gr) %in% paste0("chr", c(1:22))) & (gr$group == "hg19_cpg_shelves"))], type = "end"))))
  st_x_she <- length(unique(queryHits(findOverlaps(temp, gr[which((strand(gr) == "+") & (seqnames(gr) == "chrX") & (gr$group == "hg19_cpg_shelves"))], type = "start"))))
  en_x_she <- length(unique(queryHits(findOverlaps(temp, gr[which((strand(gr) == "-") & (seqnames(gr) == "chrX") & (gr$group == "hg19_cpg_shelves"))], type = "end"))))
  
  st_som_in <- length(unique(queryHits(findOverlaps(temp, gr[which((strand(gr) == "+") & (seqnames(gr) %in% paste0("chr", c(1:22))) & (gr$group == "hg19_cpg_inter"))], type = "start"))))
  en_som_in <- length(unique(queryHits(findOverlaps(temp, gr[which((strand(gr) == "-") & (seqnames(gr) %in% paste0("chr", c(1:22))) & (gr$group == "hg19_cpg_inter"))], type = "end"))))
  st_x_in <- length(unique(queryHits(findOverlaps(temp, gr[which((strand(gr) == "+") & (seqnames(gr) == "chrX") & (gr$group == "hg19_cpg_inter"))], type = "start"))))
  en_x_in <- length(unique(queryHits(findOverlaps(temp, gr[which((strand(gr) == "-") & (seqnames(gr) == "chrX") & (gr$group == "hg19_cpg_inter"))], type = "end"))))
 
  end(gr[which(strand(gr) == "+")]) <- start(gr[which(strand(gr) == "+")]) + 50
  start(gr[which(strand(gr) == "+")]) <- start(gr[which(strand(gr) == "+")]) - 50
  start(gr[which(strand(gr) == "-")]) <- end(gr[which(strand(gr) == "-")]) - 50
  end(gr[which(strand(gr) == "-")]) <- end(gr[which(strand(gr) == "-")]) + 50
  
  over_som_is <- length(queryHits(findOverlaps(temp, gr[which((seqnames(gr) %in% paste0("chr", c(1:22))) & (gr$group == "hg19_cpg_islands"))])))
  over_x_is <- length(queryHits(findOverlaps(temp, gr[which((seqnames(gr) == "chrX") & (gr$group == "hg19_cpg_islands"))])))
  
  over_som_sho <- length(queryHits(findOverlaps(temp, gr[which((seqnames(gr) %in% paste0("chr", c(1:22))) & (gr$group == "hg19_cpg_shores"))])))
  over_x_sho <- length(queryHits(findOverlaps(temp, gr[which((seqnames(gr) == "chrX") & (gr$group == "hg19_cpg_shores"))])))

  over_som_she <- length(queryHits(findOverlaps(temp, gr[which((seqnames(gr) %in% paste0("chr", c(1:22))) & (gr$group == "hg19_cpg_shelves"))])))
  over_x_she <- length(queryHits(findOverlaps(temp, gr[which((seqnames(gr) == "chrX") & (gr$group == "hg19_cpg_shelves"))])))
  
  over_som_in <- length(queryHits(findOverlaps(temp, gr[which((seqnames(gr) %in% paste0("chr", c(1:22))) & (gr$group == "hg19_cpg_inter"))])))
  over_x_in <- length(queryHits(findOverlaps(temp, gr[which((seqnames(gr) == "chrX") & (gr$group == "hg19_cpg_inter"))])))
  
  rat <-  length(temp[which(seqnames(temp) == "chrX")]) / length(temp)
  
  f_temp <- tibble("case" = basename(h[case]),
                   "motif" = motif,
                   "chromosome" = c(rep("somatic", 4), rep("chrX", 4)),
                   "group" = rep(c("islands", "shores", "shelves", "open sea"), 2),
                   "ratio" = c(((st_som_is + en_som_is) / over_som_is),
                               ((st_som_sho + en_som_sho) / over_som_sho),
                               ((st_som_she + en_som_she) / over_som_she),
                               ((st_som_in + en_som_in) / over_som_in),
                               ((st_x_is + en_x_is) / over_x_is),
                               ((st_x_sho + en_x_sho) / over_x_sho),
                               ((st_x_she + en_x_she) / over_x_she),
                               ((st_x_in + en_x_in) / over_x_in)),
                   "rel_chrx" = rat)
  f <- rbind(f, f_temp)
}
saveRDS(f, paste0("/dcs04/scharpf/data/mnoe/som_chrx/", basename(h[case])))
