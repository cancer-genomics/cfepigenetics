library(BSgenome.Hsapiens.UCSC.hg19)
library(rtracklayer)
library(GenomicRanges)
library(tidyverse)
library(annotatr)

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
div <- (width * 2) / 500

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

otemp <- as_tibble(data.frame(findOverlaps(cg.loc, tssgrspec)))
otemp$temppos <- start(cg.loc)[otemp$queryHits]
otemp$tssgrspecpos <- start(tssgrspec)[otemp$subjectHits]
otemp$beta_group <- tssgrspec$group[otemp$subjectHits]
otemp$beta <- cg.loc$beta[otemp$queryHits]
otemp$relpos <- otemp$temppos - otemp$tssgrspecpos + 1
otemp <- otemp[order(otemp$relpos),]
otemp$relpos_500 <- 0
otemp[which(otemp$relpos %in% c(1:width)),]$relpos_500 <- round((c(otemp[which(otemp$relpos %in% c(1:width)),]$relpos) + (div/2) - 0.00001)/div)
otemp[which(otemp$relpos == (width + 1)),]$relpos_500 <- 251
otemp[which(otemp$relpos %in% c((width + 2):(2*width+1))),]$relpos_500 <- round((c(otemp[which(otemp$relpos %in% c((width + 2):(2*width+1))),]$relpos) + (div/2) - 1.00001)/div) + 1
#otemp[which(otemp$strand == "-"),]$relpos_500 <- - otemp[which(otemp$strand == "-"),]$relpos_500 + 502
final <- otemp %>% 
  group_by(subjectHits, relpos_500, beta_group) %>%
  summarize(mean_beta = mean(beta))
saveRDS(final, "../../data/cpg_500000_topbot1000_beta.rds")





