library(BSgenome.Hsapiens.UCSC.hg19)
library(rtracklayer)
library(GenomicRanges)
library(tidyverse)


args = commandArgs(trailingOnly=TRUE)
case <- as.integer(args[1])

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


width <- 500000
end(tssgrspec[which(strand(tssgrspec) == "+")]) <- start(tssgrspec[which(strand(tssgrspec) == "+")])
start(tssgrspec[which(strand(tssgrspec) == "+")]) <- start(tssgrspec[which(strand(tssgrspec) == "+")]) - width
end(tssgrspec[which(strand(tssgrspec) == "+")]) <- end(tssgrspec[which(strand(tssgrspec) == "+")]) + width

start(tssgrspec[which(strand(tssgrspec) == "-")]) <- end(tssgrspec[which(strand(tssgrspec) == "-")])
start(tssgrspec[which(strand(tssgrspec) == "-")]) <- start(tssgrspec[which(strand(tssgrspec) == "-")]) - width
end(tssgrspec[which(strand(tssgrspec) == "-")]) <- end(tssgrspec[which(strand(tssgrspec) == "-")]) + width

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
  otemp$expression <- tssgrspec$expression[otemp$subjectHits]
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
    group_by(subjectHits, relpos_500) %>%
    summarize(mean_cov = mean(cov), mean_size = mean(size), expression = mean(expression))
  final <- rbind(final, prefinal)
}
saveRDS(final, paste0("../../../cfepigenetics_data/tss_500000_topbot1000/subset_",case, ".rds"))
