library(BSgenome.Hsapiens.UCSC.hg19)
library(rtracklayer)
library(GenomicRanges)
library(tidyverse)

tss <- read.table(file = "../../data/transcriptAnno-GRCh37.75.tsv")
tss <- tss[which(tss$V2 %in% c(1:22)),]
tss$V2 <- paste0("chr", tss$V2)
tissue_key <- read.table("../../data/tissue_key.tsv",as.is=T,header=T,sep="\t",quote="\"")
conv <- read.table('../../data/labels.txt',as.is=T,header=T,sep="\t",quote="\"")
rna <- read.table('../../data/RNAtable.tsv.gz',as.is=T,header=T)
express <- conv[which(conv$Category %in% c("Myeloid")),]$RName
express_df <- rna[c("GeneID", express)]
express_df$mean <- rowMeans(express_df[,c(2:(length(express) + 1))], na.rm=T)
express_df <- express_df[order(-express_df$mean),]
geneid_ordered <- express_df$GeneID
tss$expression <- 0
tss <- tss[which(tss$V1 %in% express_df$GeneID),]
tss[match(express_df[which(express_df$GeneID %in% tss$V1),]$GeneID, tss$V1),]$expression <- express_df[which(express_df$GeneID %in% tss$V1),]$mean
tss <- tss[order(-tss$expression),]
colnames(tss) <- c("transcript", "seqnames", "start", "end", "strand", "expression")
tssgr <- GRanges(paste0(tss$seqnames, ":", tss$start, "-", tss$end, ":", tss$strand))
tssgr$expression <- 0
tssgr$expression <- tss$expression
tssgrspec <- tssgr
tssgrspec <- tssgrspec[c(1:1000, (length(tssgrspec)-999):(length(tssgrspec)))]
tssgrspec$group <- "NA"
tssgrspec[which(tssgrspec$expression >40)]$group <- "tss_highexp"
tssgrspec[which(tssgrspec$expression <40)]$group <- "tss_lowexp"

width <- 500000
end(tssgrspec[which(strand(tssgrspec) == "+")]) <- start(tssgrspec[which(strand(tssgrspec) == "+")])
start(tssgrspec[which(strand(tssgrspec) == "+")]) <- start(tssgrspec[which(strand(tssgrspec) == "+")]) - width
end(tssgrspec[which(strand(tssgrspec) == "+")]) <- end(tssgrspec[which(strand(tssgrspec) == "+")]) + width

start(tssgrspec[which(strand(tssgrspec) == "-")]) <- end(tssgrspec[which(strand(tssgrspec) == "-")])
start(tssgrspec[which(strand(tssgrspec) == "-")]) <- start(tssgrspec[which(strand(tssgrspec) == "-")]) - width
end(tssgrspec[which(strand(tssgrspec) == "-")]) <- end(tssgrspec[which(strand(tssgrspec) == "-")]) + width

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
otemp$strand <- as.character(strand(tssgrspec))[otemp$subjectHits]
otemp$exp_group <- tssgrspec$group[otemp$subjectHits]
otemp$beta <- cg.loc$beta[otemp$queryHits]
otemp$relpos <- otemp$temppos - otemp$tssgrspecpos + 1
otemp <- otemp[order(otemp$relpos),]
otemp$relpos_500 <- 0
otemp[which(otemp$relpos %in% c(1:width)),]$relpos_500 <- round((c(otemp[which(otemp$relpos %in% c(1:width)),]$relpos) + (div/2) - 0.00001)/div)
otemp[which(otemp$relpos == (width + 1)),]$relpos_500 <- 251
otemp[which(otemp$relpos %in% c((width + 2):(2*width+1))),]$relpos_500 <- round((c(otemp[which(otemp$relpos %in% c((width + 2):(2*width+1))),]$relpos) + (div/2) - 1.00001)/div) + 1
otemp[which(otemp$strand == "-"),]$relpos_500 <- - otemp[which(otemp$strand == "-"),]$relpos_500 + 502
final <- otemp %>% 
  group_by(subjectHits, relpos_500, exp_group) %>%
  summarize(mean_beta = mean(beta))
saveRDS(final, "../../data/tss_500000_topbot1000_beta.rds")





