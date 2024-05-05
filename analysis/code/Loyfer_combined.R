library("GenomicRanges")
library("BSgenome.Hsapiens.UCSC.hg19")
library("tidyverse")
library("Biostrings")

args = commandArgs(trailingOnly=TRUE)
case <- as.integer(args[1])

cg <- DNAString("CG")
cg.loc <- vmatchPattern(cg, Hsapiens)
cg.loc <- cg.loc[which((seqnames(cg.loc) %in% c(paste0("chr", c(1:22, "X", "Y", "M")))) & (strand(cg.loc) == "+"))]

ls <- list.files("../../data", pattern = "CpG_Loyfer")
l <- tibble()
for (file in ls) {
  temp <- readRDS(file.path("../../data", file))
  l <- rbind(l, temp)
}
cg.loc <- cg.loc[l$index]
cg.loc$beta <- l$beta

loc <- "../../../cfepigenetics_data/reads_per_tilesmall"
ls <- list.files(loc)
g <- sub(".rds", "", ls)
g <- sub("_", ":", g)
g <- sub("_", "-", g)
gr <- GRanges(g)
ls <- ls[order(gr)]
gr <- gr[order(gr)]

seqlevels(cg.loc) <- seqlevels(gr)
seqinfo(cg.loc) <- seqinfo(gr)
o <- data.frame(findOverlaps(cg.loc, gr))
o <- o[order(o$subjectHits),]
ou <- unique(o$subjectHits)
oui <- round((1:length(ou))/(length(ou)/150)+0.5)

end(cg.loc) <- start(cg.loc)

final <- tibble()
for (i in ou[which(oui == case)]) {
  print(i)
  temp <- readRDS(file.path(loc, ls[i]))
  cg.loc_temp <- cg.loc[queryHits(findOverlaps(cg.loc, temp))]
  tempall <- temp
  if (i > 1) {
    if (end(gr[i-1]) == (start(gr[i]) - 1)) {
      tempmin <- readRDS(file.path(loc, ls[i-1]))
      tempall <- c(tempmin, tempall)
    }
  }
  if (i < length(ls)) {
    if (start(gr[i+1]) == (end(gr[i]) + 1)) {
      tempplus <- readRDS(file.path(loc,ls[i+1]))
      tempall <- c(tempall,tempplus)
    }
  }
  temp <- tempall
  temp_seq <- temp
  start(temp_seq) <- start(temp_seq) - 2
  end(temp_seq) <- end(temp_seq) + 3
  seqlevels(temp_seq) <- seqlevels(Hsapiens)
  seqinfo(temp_seq) <- seqinfo(Hsapiens)
  temp_seq <- trim(temp_seq)
  temp$seq <- as.character(getSeq(Hsapiens, temp_seq))
  temp <- temp[which(nchar(temp$seq) == 6)]
  otemp <- as_tibble(as.data.frame(findOverlaps(cg.loc_temp, temp)))
  otemp$pos <- paste0(seqnames(cg.loc_temp[otemp$queryHits]), ":", start(cg.loc_temp[otemp$queryHits]))
  otemp$beta <- cg.loc_temp[otemp$queryHits]$beta
  otemp$seq <- substr(temp[otemp$subjectHits]$seq, 2, 4)
  otemp$start <- temp[otemp$subjectHits]$start_reads
  otemp$over <- temp[otemp$subjectHits]$overlap
  final <- rbind(final, otemp[,c(3:7)])
  
  start(cg.loc_temp) <- start(cg.loc_temp) - 1
  end(cg.loc_temp) <- start(cg.loc_temp)
  otemp <- as_tibble(as.data.frame(findOverlaps(cg.loc_temp, temp)))
  otemp$pos <- paste0(seqnames(cg.loc_temp[otemp$queryHits]), ":", start(cg.loc_temp[otemp$queryHits]) + 1)
  otemp$beta <- cg.loc_temp[otemp$queryHits]$beta
  otemp$seq <- substr(temp[otemp$subjectHits]$seq, 2, 5)
  otemp$start <- temp[c(otemp$subjectHits)]$start_reads
  otemp$over <- temp[c(otemp$subjectHits)]$overlap
  final <- rbind(final, otemp[,c(3:7)])
  
  end(cg.loc_temp) <- end(cg.loc_temp) + 2
  start(cg.loc_temp) <- end(cg.loc_temp)
  otemp <- as_tibble(as.data.frame(findOverlaps(cg.loc_temp, temp)))
  otemp$pos <- paste0(seqnames(cg.loc_temp[otemp$queryHits]), ":", start(cg.loc_temp[otemp$queryHits]) - 1 )
  otemp$beta <- cg.loc_temp[otemp$queryHits]$beta
  otemp$seq <- as.character(reverseComplement(DNAStringSet(substr(temp[otemp$subjectHits]$seq, 2, 4))))
  otemp$start <- temp[c(otemp$subjectHits)]$start_reads
  otemp$over <- temp[c(otemp$subjectHits)]$overlap
  final <- rbind(final, otemp[,c(3:7)])
  
  end(cg.loc_temp) <- end(cg.loc_temp) + 1
  start(cg.loc_temp) <- end(cg.loc_temp)
  otemp <- as_tibble(as.data.frame(findOverlaps(cg.loc_temp, temp)))
  otemp$pos <- paste0(seqnames(cg.loc_temp[otemp$queryHits]), ":", start(cg.loc_temp[otemp$queryHits]) - 2 )
  otemp$beta <- cg.loc_temp[otemp$queryHits]$beta
  otemp$seq <- as.character(reverseComplement(DNAStringSet(substr(temp[otemp$subjectHits]$seq, 1, 4))))
  otemp$start <- temp[c(otemp$subjectHits)]$start_reads
  otemp$over <- temp[c(otemp$subjectHits)]$overlap
  final <- rbind(final, otemp[,c(3:7)])
}
saveRDS(final, paste0("../../../cfepigenetics_data/Loyfer_combined/sub_", case, ".rds"))



