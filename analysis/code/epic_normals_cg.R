library("minfi")
library("IlluminaHumanMethylationEPICmanifest")
library("BSgenome.Hsapiens.UCSC.hg19")
library("GenomicRanges")

path <- "../../cfdna_methylation_data/epic/Moss_et_al"
total <- readRDS(file.path(path, "results", "healthy.rds"))

val_files <- readRDS("../data/selected_validation_files.rds")
luc_files <- readRDS("../data/selected_lucas_files.rds")
val <- val_files[which(val_files$patient_type == "healthy"),]$pgdx_id
luc_h <- luc_files[which(luc_files$patient_type == 2),]$pgdx_id
h <- c(val, luc_h)
c <- luc_files[which(luc_files$patient_type == 3),]$pgdx_id
all <- c(h,c)

total_gr <- GRanges(seqnames = total$chr, IRanges(start = total$pos, end = total$pos + 1))
seq <- as.character(getSeq(Hsapiens, total_gr))
total <- total[which(seq == "CG"),]
total_gr <- total_gr[which(seq == "CG"),]
total_gr_large <- total_gr
start(total_gr_large) <- start(total_gr_large) - 50
end(total_gr_large) <- end(total_gr_large) + 50
total$start <- 0
total$end <- 0
total$over <- 0
for (file in h) {
  print(paste0(match(file,h), "/", length(h)))
  temp <- readRDS(file)
  over <- findOverlaps(temp, total_gr, type = "start")
  tover <- table(subjectHits(over))
  total[as.numeric(names(tover)),]$start <- total[as.numeric(names(tover)),]$start + as.numeric(tover)
  over <- findOverlaps(temp, total_gr, type = "end")
  tover <- table(subjectHits(over))
  total[as.numeric(names(tover)),]$end <- total[as.numeric(names(tover)),]$end + as.numeric(tover)
  over <- findOverlaps(temp, total_gr_large)
  tover <- table(subjectHits(over))
  total[as.numeric(names(tover)),]$over <- total[as.numeric(names(tover)),]$over + as.numeric(tover)
}
saveRDS(total, "../../cfdna_methylation_data/cg/start_over_epic.rds")

