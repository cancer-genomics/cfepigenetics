library("BSgenome.Hsapiens.UCSC.hg19")
library("GenomicRanges")
library("annotatr")

annots <- c('hg19_cpgs')
annots <- annotatr::build_annotations(genome = 'hg19', annotations = annots)

args = commandArgs(trailingOnly=TRUE)
case <- as.integer(args[1])

motifs <- c("ACG", "CCG", "GCG", "TCG",
            "ACCG", "CCCG", "GCCG", "TCCG")
motif <- motifs[case]
print(motif)
chrs <- c(paste0("chr", c(1:22, "X")))
temp <- vmatchPattern(motif, Hsapiens)
temp <- temp[which(seqnames(temp) %in% chrs)]

over <- findOverlaps(temp, annots)
temp <- temp[queryHits(over)]
temp$group <- annots[subjectHits(over)]$type

saveRDS(temp, paste0("/dcs04/scharpf/data/mnoe/motif/", motif, ".rds"))
