library("BSgenome.Hsapiens.UCSC.hg19")
library("GenomicRanges")
library("tidyverse")
library("stringr")
library("pracma")

args = commandArgs(trailingOnly=TRUE)
case <- as.integer(args[1])

tiles <- unlist(tileGenome(seqinfo(Hsapiens), tilewidth = 100000, cut.last.tile.in.chrom=FALSE))
tiles <- tiles[which(seqnames(tiles) %in% paste0("chr", c(1:22)))]

ls <- list.files("../../../cfepigenetics_data/reads_per_tilebig", full.names=T)
df <- data.frame("loc" = ls)
df$nam <- sub(".rds", "", basename(df$loc))
df$seqnames <- sapply(str_split(df$nam, "_"), "[[", 1)
df$start <- as.numeric(sapply(str_split(df$nam, "_"), "[[", 2))
df$end <- as.numeric(sapply(str_split(df$nam, "_"), "[[", 3))
df$order <- as.numeric(sapply(str_split(df$nam, "_"), "[[", 4))
gr <- GRanges(paste0(df$seqnames, ":", df$start, "-", df$end))
gr$or <- sub("chr", "", seqnames(gr))
gr$or <- as.numeric(as.character(gr$or))
df <- df[order(gr$or, start(gr)),]
gr <- gr[order(gr$or, start(gr))]
tiles_big <- gr
o <- as.data.frame(findOverlaps(tiles, tiles_big, type="within"))
#o <- o[!duplicated(o[,1]),]
tiles$group <- 0
tiles[o$queryHits,]$group <- o$subjectHits
tile <- tiles[which(tiles$group == case)]
tile_big <- GRanges(paste0(unique(seqnames(tiles_big[case])), ":", start(tiles_big[case]), "-", end(tiles_big[case])))

temp <- readRDS(paste0("../../../cfepigenetics_data/reads_per_tilebig/", seqnames(tiles_big[case]), "_", start(tiles_big[case]), "_", end(tiles_big[case]), "_", df[case,]$order, ".rds"))
for (i in c(1:length(tile))) {
  print(i)
  tiletemp <- tile[i]
  tiletemp_start <- start(tiletemp)
  tiletemp_end <- end(tiletemp)
  tiletemp_wide <- GRanges(paste0(seqnames(tiletemp), ":", (tiletemp_start - 2000), "-", (tiletemp_end + 2000)))
  if (length(findOverlaps(tiletemp_wide, tiles_big)) > 1) {
    ot <- findOverlaps(tiletemp_wide, tiles_big)
    other_tiles_big <- tiles_big[subjectHits(ot)][-match(tiles_big[case], tiles_big[subjectHits(ot)])]
    for (v in 1:length(other_tiles_big)) {
      temp2 <- paste0("../../../cfepigenetics_data/reads_per_tilebig/", seqnames(other_tiles_big[v]), "_", start(other_tiles_big[v]), "_", end(other_tiles_big[v]))
      temp2_file <- readRDS(ls[grep(temp2, ls)])
      ot_temp2 <- findOverlaps(temp2_file, tiletemp_wide)
      temp <- c(temp,temp2_file[queryHits(ot_temp2)])
      rm(temp2_file)
    }
  }
  tiletemp <- GRanges(paste0(seqnames(tiletemp), ":", c((start(tiletemp)-2000):(end(tiletemp)+2000))))
  tiletemp_window <- tiletemp
  start(tiletemp_window) <- start(tiletemp_window) - 60
  end(tiletemp_window) <- end(tiletemp_window) + 60
  o <- as_tibble(findOverlaps(tiletemp, temp))
  o$width <- width(temp)[o$subjectHits]
  t <- table(o$queryHits)
  tiletemp$cov <- 0
  tiletemp$size <- 0
  tiletemp$wps <- 0
  tiletemp$wps_ratio <- 0
  tiletemp$start_reads <- 0
  tiletemp$end_reads <- 0
  tiletemp$overlap <- 0 
  tiletemp$seq <- "NA"
  tiletemp[as.numeric(as.character(names(t)))]$cov <- as.numeric(t)
  ow <- o %>%
    group_by(queryHits) %>%
    summarize(width_sum = sum(width))
  tiletemp[ow$queryHits]$size <- ow$width_sum
  tiletemp[which(tiletemp$cov !=0)]$size <- tiletemp[which(tiletemp$cov !=0)]$size / tiletemp[which(tiletemp$cov !=0)]$cov
  ow <- table(queryHits(findOverlaps(tiletemp_window, temp, type = "within")))
  oa <- table(queryHits(findOverlaps(tiletemp_window, temp)))
  tiletemp_window$within  <- 0
  tiletemp_window[as.numeric(as.character(names(ow)))]$within  <- as.numeric(ow)
  tiletemp_window$all  <- 0
  tiletemp_window[as.numeric(as.character(names(oa)))]$all  <- as.numeric(oa)
  tiletemp$wps <- tiletemp_window$within - tiletemp_window$all
  tiletemp$wps_ratio <- tiletemp_window$within / tiletemp_window$all
  if (length(tiletemp[which(is.na(tiletemp$wps_ratio))]$wps_ratio) > 0) {
    tiletemp[which(is.na(tiletemp$wps_ratio))]$wps_ratio <- 0
  }
  o <- as_tibble(findOverlaps(tiletemp, temp, type = "start"))
  t <- table(o$queryHits)
  tiletemp[as.numeric(as.character(names(t)))]$start_reads <- as.numeric(t)
  o <- as_tibble(findOverlaps(tiletemp, temp, type = "end"))
  t <- table(o$queryHits)
  tiletemp[as.numeric(as.character(names(t)))]$end_reads <- as.numeric(t)
  tiletemp_window <- tiletemp
  start(tiletemp_window) <- start(tiletemp_window) - 50
  end(tiletemp_window) <- end(tiletemp_window) + 50
  oa <- as_tibble(findOverlaps(tiletemp_window, temp))
  t <- table(oa$queryHits)
  tiletemp[as.numeric(as.character(names(t)))]$overlap <- as.numeric(t)
  tiletemp_window <- tiletemp
  start(tiletemp_window) <- start(tiletemp_window) - 2
  end(tiletemp_window) <- end(tiletemp_window) + 3
  seqlevels(tiletemp_window) <- seqlevels(Hsapiens)
  seqinfo(tiletemp_window) <- seqinfo(Hsapiens)
  tiletemp_window <- trim(tiletemp_window)
  tiletemp$seq <- getSeq(Hsapiens, tiletemp_window)
  twps<- savgol(tiletemp$wps, fl = 21, forder=2, dorder=0)
  tiletemp$wps_adj <- twps
  twide <- tiletemp
  start(twide) <- start(twide) - 500
  end(twide) <- end(twide) + 500
  o <- as_tibble(data.frame(findOverlaps(twide, tiletemp)))
  o$wps_adj <- tiletemp[o$subjectHits]$wps_adj
  o <- o %>%
    group_by(queryHits) %>%
    summarize(wps_adj_mean = mean(wps_adj))
  tiletemp$wps_adj_mean <- 0
  tiletemp[o$queryHits]$wps_adj_mean <- o$wps_adj_mean
  tiletemp$wps_adj_norm <- 0
  tiletemp$wps_adj_norm <- tiletemp$wps_adj - tiletemp$wps_adj_mean
  tiletemp <- tiletemp[which((start(tiletemp) >= tiletemp_start) & (start(tiletemp) <= tiletemp_end))]
  saveRDS(tiletemp, paste0("../../../cfepigenetics_data/reads_per_tilesmall/", seqnames(tile[i]), "_", start(tile[i]), "_", end(tile[i]), ".rds"))
}

