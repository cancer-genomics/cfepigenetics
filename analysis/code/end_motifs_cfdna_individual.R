library("GenomicRanges")
library("BSgenome.Hsapiens.UCSC.hg19")
library("BSgenome")
library("tidyverse")
library("here")

args = commandArgs(trailingOnly=TRUE)
case <- as.integer(args[1])

val_files <- readRDS("../../data/selected_validation_files.rds")
luc_files <- readRDS("../../data/selected_lucas_files.rds")
val <- val_files[which(val_files$patient_type == "healthy"),]$pgdx_id
luc_h <- luc_files[which(luc_files$patient_type == 2),]$pgdx_id
h <- c(val, luc_h)
c <- luc_files[which(luc_files$patient_type == 3),]$pgdx_id
all <- c(h,c)

file <- h[case]
temp <- readRDS(file)
start_gr <- temp
start(start_gr) <- start(start_gr) -1
end(start_gr) <- start(start_gr) + 2 
strand(start_gr) <- "+"
start_tb <- as_tibble(as.data.frame(table(as.character(getSeq(Hsapiens, start_gr)))))
end_gr <- temp
start(end_gr) <- end(end_gr) -1
end(end_gr) <- start(end_gr) + 2
strand(end_gr) <- "-"
end_tb <- as_tibble(as.data.frame(table(as.character(getSeq(Hsapiens, end_gr)))))

tb <- rbind(start_tb, end_tb)
tb <- tb %>%
  group_by(Var1) %>%
  summarize(Freq = sum(Freq))

if (!dir.exists("../../../cfepigenetics_data/end_motifs_cfdna_individual")) {
  dir.create("../../../cfepigenetics_data/end_motifs_cfdna_individual")
}
saveRDS(tb, paste0("../../../cfepigenetics_data/end_motifs_cfdna_individual/", basename(file)))
