library(BSgenome.Hsapiens.UCSC.hg19)
library(rtracklayer)
library(GenomicRanges)
library(tidyverse)
library(annotatr)

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
over$n <- 1
over <- over %>%
  group_by(subjectHits) %>%
  summarize(mean_beta = mean(beta), n = sum(n))
annotations$beta <- 2
annotations[over$subjectHits,]$beta <- over$mean_beta
annotations$amount <- 0
annotations[over$subjectHits,]$amount <- over$n
annotations <- annotations[which(annotations$beta != 2)]
annotations <- annotations[which(annotations$amount > 2)]
annotations <- annotations[order(annotations$beta)]
annotations_low <- annotations[c(1:1000)]
annotations_low$group <- "cpg_lowmet"
annotations_high <- annotations[order(-annotations$beta)][c(1:1000)]
annotations_high$group <- "cpg_highmet"
tssgr <- c(annotations_high, annotations_low)
start(tssgr) <- (c(end(tssgr) - start(tssgr))/2) + start(tssgr)
end(tssgr) <- start(tssgr)
tssgrspec <- tssgr

width <- 500000
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

o <- data.frame(findOverlaps(tssgrspec, gr))
o <- o[order(o$subjectHits),]
ou <- unique(o$subjectHits)
oui <- round((1:length(ou))/(length(ou)/150)+0.5)

div <- (width * 2) / 500

final <- tibble()

for (i in ou[which(oui == case)]) {
  print(paste0(match(i, unique(o$subjectHits)), "/", length(unique(o$subjectHits)), " - ", i))
  temp <- readRDS(file.path(loc, ls[i]))
  
  otemp <- as_tibble(data.frame(findOverlaps(temp, tssgrspec)))
  otemp$temppos <- start(temp)[otemp$queryHits]
  otemp$tssgrspecpos <- start(tssgrspec)[otemp$subjectHits]
  otemp$beta <- tssgrspec$group[otemp$subjectHits]
  otemp$strand <- as.character(strand(tssgrspec))[otemp$subjectHits]
  otemp$relpos <- otemp$temppos - otemp$tssgrspecpos + 1
  otemp <- otemp[order(otemp$relpos),]
  otemp$relpos_500 <- 0
  otemp[which(otemp$relpos %in% c(1:width)),]$relpos_500 <- round((c(otemp[which(otemp$relpos %in% c(1:width)),]$relpos) + (div/2) - 0.00001)/div)
  otemp[which(otemp$relpos == (width + 1)),]$relpos_500 <- 251
  otemp[which(otemp$relpos %in% c((width + 2):(2*width+1))),]$relpos_500 <- round((c(otemp[which(otemp$relpos %in% c((width + 2):(2*width+1))),]$relpos) + (div/2) - 1.00001)/div) + 1
  otemp[which(otemp$strand == "-"),]$relpos_500 <- - otemp[which(otemp$strand == "-"),]$relpos_500 + 502
  otemp$cov <- temp$cov[otemp$queryHits]
  otemp$size <- temp$size[otemp$queryHits]
  prefinal <- otemp %>% 
    group_by(subjectHits, relpos_500, beta) %>%
    summarize(mean_cov = mean(cov), mean_size = mean(size))
  final <- rbind(final, prefinal)
}
saveRDS(final, paste0("../../../cfepigenetics_data/cpg_500000_topbot1000/subset_",case, ".rds"))
