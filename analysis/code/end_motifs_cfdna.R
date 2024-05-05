library("GenomicRanges")
library("BSgenome.Hsapiens.UCSC.hg19")
library("BSgenome")
library("here")

args = commandArgs(trailingOnly=TRUE)
case <- as.integer(args[1])

chromosomes <- paste0("chr", c(1:22))
tilegr <- unlist(tileGenome(seqlengths(Hsapiens), tilewidth = 20000000))
tilegr <- tilegr[which(seqnames(tilegr) %in% chromosomes)]
tilegr <- tilegr[case]
chr <- as.character(unique(seqnames(tilegr)))

val_files <- readRDS("../../data/selected_validation_files.rds")
luc_files <- readRDS("../../data/selected_lucas_files.rds")
val <- val_files[which(val_files$patient_type == "healthy"),]$pgdx_id
luc_h <- luc_files[which(luc_files$patient_type == 2),]$pgdx_id
h <- c(val, luc_h)
c <- luc_files[which(luc_files$patient_type == 3),]$pgdx_id
all <- c(h,c)

collect <- data.frame("seqnames" = chr,
                      "location" = c(start(tilegr) : end(tilegr)),
                      "start" = 0,
                      "end" = 0,
                      "over" = 0)
collect_gr <- GRanges(seqnames = chr, IRanges(start = collect$location - 50, end = collect$location + 50))
collect_gr <- collect_gr[which((start(collect_gr) >= 1) & (end(collect_gr) <= as.numeric(seqlengths(Hsapiens)[chr])))]

if (!dir.exists("../../../cfepigenetics_data/end_motifs_cfdna")) {
  dir.create("../../../cfepigenetics_data/end_motifs_cfdna")
}

for (file in h) {
  print(paste0(match(file, h), "/", length(h)))
  temp <- readRDS(file)
  temp <- temp[which(seqnames(temp) == chr)]
  temp <- temp[which((end(temp) >= (start(tilegr) - 51)) & (start(temp) <= (end(tilegr) + 51)))]
  
  start <- table(start(temp))
  start <- start[which(names(start) %in% collect$location)]
  collect[match(as.numeric(names(start)), collect$location),]$start <- collect[match(as.numeric(names(start)), collect$location),]$start + as.numeric(start)
  end <- table(end(temp))
  end <- end[which(names(end) %in% collect$location)]
  collect[match(as.numeric(names(end)), collect$location),]$end <- collect[match(as.numeric(names(end)), collect$location),]$end + as.numeric(end)
  table_over <- table(queryHits(findOverlaps(collect_gr, temp)))
  collect[match(start(collect_gr[as.numeric(names(table_over))]) + 50, collect$location),]$over <- collect[match(start(collect_gr[as.numeric(names(table_over))]) + 50, collect$location),]$over + as.numeric(table_over)
  if (match(file, h) %in% seq(10, 600, 10)) {
    prev <- list.files(path="../../../cfepigenetics_data/end_motifs_cfdna", pattern=paste0(chr, "_", case, "_"), full.names = T)
    saveRDS(collect, file=paste0("../../../cfepigenetics_data/end_motifs_cfdna/", chr, "_", case, "_", match(file, h), ".rds"))
    if (length(prev) != 0) {
      file.remove(prev)
    }
  }
}
prev <- list.files(path="../../../cfepigenetics_data/end_motifs_cfdna", pattern=paste0(chr, "_", case, "_"), full.names = T)
saveRDS(collect, file=paste0("../../../cfepigenetics_data/end_motifs_cfdna/", chr, "_", case, "_full.rds"))
if (length(prev) != 0) {
  file.remove(prev)
}
print("success")


