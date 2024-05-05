library(rtracklayer)
library(GenomicRanges)
library(tidyverse)
library(annotatr)
library(BSgenome.Hsapiens.UCSC.hg19)

args = commandArgs(trailingOnly=TRUE)
case <- as.integer(args[1])

he <- readRDS("../../../cfepigenetics_data/Moss/results/healthy.rds")
he <- he[which(he$chr %in% paste0("chr", c(1:22))),]
hegr <- GRanges(paste0(he$chr, ":", he$pos, "-", he$pos+1))
hese <- as.character(getSeq(Hsapiens, hegr))
hegr$beta <- rowMeans(he[,c(4:11)])
hegr <- hegr[which(hese == "CG")]

annots = c('hg19_cpgs')
annotations = build_annotations(genome = 'hg19', annotations = annots)
annotations <- annotations[which(seqnames(annotations) %in% paste0("chr", c(1:22)))]
annotations <- annotations[which(annotations$type == "hg19_cpg_islands")]
over <- as_tibble(data.frame(findOverlaps(hegr, annotations)))
over$beta <- hegr[over$queryHits]$beta
over <- over %>%
  group_by(subjectHits) %>%
  summarize(mean_beta = mean(beta))
annotations$beta <- 2
annotations[over$subjectHits,]$beta <- over$mean_beta
annotations <- annotations[which(annotations$beta != 2)]
annotations <- annotations[order(annotations$beta)]
tssgr <- annotations
start(tssgr) <- (c(end(tssgr) - start(tssgr))/2) + start(tssgr)
end(tssgr) <- start(tssgr)
tssgrspec <- tssgr

width <- 2500
start(tssgrspec) <- start(tssgrspec) - width
end(tssgrspec) <- end(tssgrspec) + width

loc <- "../../../cfepigenetics_data/reads_per_tilesmall"
ls <- list.files(loc)
g <- sub(".rds", "", ls)
g <- sub("_", ":", g)
g <- sub("_", "-", g)
gr <- GRanges(g)
ls <- ls[order(gr)]
gr <- gr[order(gr)]

seqlevels(tssgrspec) <- seqlevels(gr)
seqinfo(tssgrspec) <- seqinfo(gr)
o <- data.frame(findOverlaps(tssgrspec, gr))
o <- o[order(o$subjectHits),]
ou <- unique(o$subjectHits)
oui <- round((1:length(ou))/(length(ou)/150)+0.5)

div <- (width * 2) / 500

final <- tibble()

for (i in ou[which(oui == case)]) {
  print(paste0(match(i, unique(o$subjectHits)), "/", length(unique(o$subjectHits)), " - ", i))
  temp <- readRDS(file.path(loc, ls[i]))
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
      tempall <- c(tempall, tempplus)
    }
  }
  otemp <- as_tibble(data.frame(findOverlaps(temp, tssgrspec)))
  #otemp$temppos <- 0
  otemp$temppos <- start(temp)[otemp$queryHits]
  otemp$tssgrspecpos <- start(tssgrspec)[otemp$subjectHits]
  otemp$strand <- as.character(strand(tssgrspec))[otemp$subjectHits]
  otemp$relpos <- otemp$temppos - otemp$tssgrspecpos + 1
  otemp <- otemp[order(otemp$relpos),]
  otemp$relpos_500 <- 0
  otemp[which(otemp$relpos %in% c(1:width)),]$relpos_500 <- round((c(otemp[which(otemp$relpos %in% c(1:width)),]$relpos) + (div/2) - 0.00001)/div)
  otemp[which(otemp$relpos == (width + 1)),]$relpos_500 <- 251
  otemp[which(otemp$relpos %in% c((width + 2):(2*width+1))),]$relpos_500 <- round((c(otemp[which(otemp$relpos %in% c((width + 2):(2*width+1))),]$relpos) + (div/2) - 1.00001)/div) + 1
  #otemp[which(otemp$strand == "-"),]$relpos_500 <- - otemp[which(otemp$strand == "-"),]$relpos_500 + 502
  otemp$cov <- temp$cov[otemp$queryHits]
  otemp$size <- temp$size[otemp$queryHits]
  otemp$wps <- temp$wps_adj_norm[otemp$queryHits]
  prefinal <- otemp %>% 
    group_by(subjectHits, relpos_500) %>%
    summarize(mean_cov = mean(cov), mean_size = mean(size), mean_wps = mean(wps))
  final <- rbind(final, prefinal)
}

saveRDS(final, paste0("../../../cfepigenetics_data/cov_size_bins_cpg2500/subtab_", case, ".rds"))
