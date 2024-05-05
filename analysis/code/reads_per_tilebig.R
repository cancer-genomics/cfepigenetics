library("BSgenome.Hsapiens.UCSC.hg19")
library("GenomicRanges")

args = commandArgs(trailingOnly=TRUE)
case <- as.integer(args[1])

val_files <- readRDS("../../data/selected_validation_files.rds")
luc_files <- readRDS("../../data/selected_lucas_files.rds")
val <- val_files[which(val_files$patient_type == "healthy"),]$pgdx_id
luc_h <- luc_files[which(luc_files$patient_type == 2),]$pgdx_id
h <- c(val, luc_h)
c <- luc_files[which(luc_files$patient_type == 3),]$pgdx_id
all <- c(h,c)

tiles <- unlist(tileGenome(seqinfo(Hsapiens), tilewidth = 100000, cut.last.tile.in.chrom=FALSE))
tiles <- tiles[which(seqnames(tiles) %in% paste0("chr", c(1:22)))]
tiles_big <- unlist(tileGenome(seqinfo(Hsapiens), tilewidth = 5000000, cut.last.tile.in.chrom=FALSE))
tiles_big <- tiles_big[which(seqnames(tiles_big) %in% paste0("chr", c(1:22)))]
o <- as.data.frame(findOverlaps(tiles, tiles_big))
o <- o[!duplicated(o[,1]),]
tiles$group <- 0
tiles[o$queryHits,]$group <- o$subjectHits
tile <- tiles[which(tiles$group == case)]
tile_big <- GRanges(paste0(unique(seqnames(tile)), ":", min(start(tile)), "-", max(end(tile))))
start(tile_big) <- start(tile_big) - 100
end(tile_big) <- end(tile_big) + 100
reads <- GRanges()
for (file in h) {
  print(paste0(match(file, h)," / ", length(h)))
  temp <- readRDS(file)
  temp <- temp[subjectHits(findOverlaps(tile_big, temp))]
  temp <- temp[which((width(temp) >= 100) & (width(temp) <=220))]
  reads <- c(reads, temp[,c()])
}
saveRDS(reads, paste0("../../../cfepigenetics_data/reads_per_tilebig/", seqnames(tile_big), "_", min(start(tile)), "_", max(end(tile)), "_", case, ".rds"))

