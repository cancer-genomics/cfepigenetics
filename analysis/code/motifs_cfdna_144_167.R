library("BSgenome.Hsapiens.UCSC.hg19")
library("GenomicRanges")
library("tidyverse")

args = commandArgs(trailingOnly=TRUE)
case <- as.integer(args[1])

loc <- "../../../cfepigenetics_data/reads_per_tilebig"
loc_out <- "../../../cfepigenetics_data/motifs_cfdna_144_167"
ls <- list.files(loc)

temp <- readRDS(file.path(loc,ls[case]))
temp144 <- temp[which(width(temp) == 144)]
temp167 <- temp[which(width(temp) == 167)]
rm(temp)
start(temp144) <- start(temp144) - 9
end(temp144) <- end(temp144) + 9
start(temp167) <- start(temp167) - 9
end(temp167) <- end(temp167) + 9

seq144 <- as.character(getSeq(Hsapiens, temp144))
seq167 <- as.character(getSeq(Hsapiens, temp167))

final <- tibble()
for (i in c(1:162)) {
  print(i)
  temp <- as_tibble(as.data.frame(table(substr(seq144, i, i))))
  temp$pos <- i
  temp$motif_length <- 1
  temp$group <- "seq144"
  final <- rbind(final, temp)
}
for (i in c(1:185)) {
  print(i)
  temp <- as_tibble(as.data.frame(table(substr(seq167, i, i))))
  temp$pos <- i
  temp$motif_length <- 1
  temp$group <- "seq167"
  final <- rbind(final, temp)
}
for (i in c(1:161)) {
  print(i)
  temp <- as_tibble(as.data.frame(table(substr(seq144, i, i+1))))
  temp$pos <- i
  temp$motif_length <- 2
  temp$group <- "seq144"
  final <- rbind(final, temp)
}
for (i in c(1:184)) {
  print(i)
  temp <- as_tibble(as.data.frame(table(substr(seq167, i, i+1))))
  temp$pos <- i
  temp$motif_length <- 2
  temp$group <- "seq167"
  final <- rbind(final, temp)
}
for (i in c(1:160)) {
  print(i)
  temp <- as_tibble(as.data.frame(table(substr(seq144, i, i+2))))
  temp$pos <- i
  temp$motif_length <- 3
  temp$group <- "seq144"
  final <- rbind(final, temp)
}
for (i in c(1:183)) {
  print(i)
  temp <- as_tibble(as.data.frame(table(substr(seq167, i, i+2))))
  temp$pos <- i
  temp$motif_length <- 3
  temp$group <- "seq167"
  final <- rbind(final, temp)
}
for (i in c(1:159)) {
  print(i)
  temp <- as_tibble(as.data.frame(table(substr(seq144, i, i+3))))
  temp$pos <- i
  temp$motif_length <- 4
  temp$group <- "seq144"
  final <- rbind(final, temp)
}
for (i in c(1:182)) {
  print(i)
  temp <- as_tibble(as.data.frame(table(substr(seq167, i, i+3))))
  temp$pos <- i
  temp$motif_length <- 4
  temp$group <- "seq167"
  final <- rbind(final, temp)
}
saveRDS(final, file.path(loc_out, ls[case]))
