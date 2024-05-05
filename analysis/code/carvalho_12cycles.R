library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg19)

args = commandArgs(trailingOnly=TRUE)
case <- as.integer(args[1])

ls <- list.files("/dcl01/scharpf1/data/cristiano/projects/DELFI_paper/granges", full.names = T)
pat <- read.table("/dcl02/leased/cglab/rscharpf/mnoe/PROJECTS/delfi_annotation/tidy_data/patient.csv", sep = ",", header = T)
sam <- read.table("/dcl02/leased/cglab/rscharpf/mnoe/PROJECTS/delfi_annotation/tidy_data/sample.csv", sep = ",", header = T)
met <- read.table("/dcl01/scharpf1/data/cristiano/projects/DELFI_paper/metadata/sample_reference.csv", sep = ",", header = T)
dmc <- readxl::read_xlsx("/dcl02/leased/cglab/rscharpf/mnoe/PROJECTS/cfdna_methylation_data/Carvalho/dmc/Supplementary_Table_3.xlsx")
dmc <- dmc[,c(1:6)]
dmc <- dmc[-1,]
colnames(dmc) <- c('chr', 'start', 'end', 'diff', 'pvalue', 'qvalue')
dmc_gr <- GRanges(paste0(dmc$chr, ":", dmc$start, "-", dmc$end))
dmc_gr$seq <- as.character(getSeq(x = Hsapiens, dmc_gr))
start(dmc_gr[which(dmc_gr$seq == "G"),]) <- start(dmc_gr[which(dmc_gr$seq == "G"),]) - 1
end(dmc_gr[which(dmc_gr$seq == "C"),]) <- end(dmc_gr[which(dmc_gr$seq == "C"),]) + 1
start(dmc_gr) <- start(dmc_gr) - 2
end(dmc_gr) <- end(dmc_gr) + 2
dmc_gr$seq <- as.character(getSeq(x = Hsapiens, dmc_gr))

file <- ls[case]
temp <- readRDS(file)
final <- GRanges()

sub <- dmc_gr[which((substr(dmc_gr$seq, 2,4) == "ACG"))]
sub$diff <- dmc[which((substr(dmc_gr$seq, 2,4) == "ACG")),]$diff
sub$pvalue <- dmc[which((substr(dmc_gr$seq, 2,4) == "ACG")),]$pvalue
sub$group <- "ACG"
sub$start_reads <- 0
sub$over_reads <- 0
sub$start_extra_reads <- 0
sub$over_extra_reads <- 0
start(sub) <- start(sub) + 2
end(sub) <- start(sub)
start_table <- table(queryHits(findOverlaps(sub, temp, type="start")))
sub_large <- sub
start(sub_large) <- start(sub_large) - 50
end(sub_large) <- end(sub_large) + 50
over_table <- table(queryHits(findOverlaps(sub_large, temp)))
start(sub) <- start(sub) - 1
end(sub) <- start(sub)
start1_table <- table(queryHits(findOverlaps(sub, temp, type="start")))
sub_large <- sub
start(sub_large) <- start(sub_large) - 50
end(sub_large) <- end(sub_large) + 50
over1_table <- table(queryHits(findOverlaps(sub_large, temp)))
sub[as.numeric(names(start_table))]$start_reads <- as.numeric(start_table)
sub[as.numeric(names(over_table))]$over_reads <- as.numeric(over_table)
sub[as.numeric(names(start1_table))]$start_extra_reads <- as.numeric(start1_table)
sub[as.numeric(names(over1_table))]$over_extra_reads <- as.numeric(over1_table)

final <- c(final, sub)

sub <- dmc_gr[which((substr(dmc_gr$seq, 3,5) == "CGT"))]
sub$diff <- dmc[which((substr(dmc_gr$seq, 3,5) == "CGT")),]$diff
sub$pvalue <- dmc[which((substr(dmc_gr$seq, 3,5) == "CGT")),]$pvalue
sub$group <- "ACG"
sub$start_reads <- 0
sub$over_reads <- 0
sub$start_extra_reads <- 0
sub$over_extra_reads <- 0
start(sub) <- start(sub) + 3
end(sub) <- start(sub)
start_table <- table(queryHits(findOverlaps(sub, temp, type="end")))
sub_large <- sub
start(sub_large) <- start(sub_large) - 50
end(sub_large) <- end(sub_large) + 50
over_table <- table(queryHits(findOverlaps(sub_large, temp)))
start(sub) <- start(sub) + 1
end(sub) <- start(sub)
start1_table <- table(queryHits(findOverlaps(sub, temp, type="end")))
sub_large <- sub
start(sub_large) <- start(sub_large) - 50
end(sub_large) <- end(sub_large) + 50
over1_table <- table(queryHits(findOverlaps(sub_large, temp)))
sub[as.numeric(names(start_table))]$start_reads <- as.numeric(start_table)
sub[as.numeric(names(over_table))]$over_reads <- as.numeric(over_table)
sub[as.numeric(names(start1_table))]$start_extra_reads <- as.numeric(start1_table)
sub[as.numeric(names(over1_table))]$over_extra_reads <- as.numeric(over1_table)

final <- c(final, sub)

sub <- dmc_gr[which((substr(dmc_gr$seq, 2,4) == "TCG"))]
sub$diff <- dmc[which((substr(dmc_gr$seq, 2,4) == "TCG")),]$diff
sub$pvalue <- dmc[which((substr(dmc_gr$seq, 2,4) == "TCG")),]$pvalue
sub$group <- "TCG"
sub$start_reads <- 0
sub$over_reads <- 0
sub$start_extra_reads <- 0
sub$over_extra_reads <- 0
start(sub) <- start(sub) + 2
end(sub) <- start(sub)
start_table <- table(queryHits(findOverlaps(sub, temp, type="start")))
sub_large <- sub
start(sub_large) <- start(sub_large) - 50
end(sub_large) <- end(sub_large) + 50
over_table <- table(queryHits(findOverlaps(sub_large, temp)))
start(sub) <- start(sub) - 1
end(sub) <- start(sub)
start1_table <- table(queryHits(findOverlaps(sub, temp, type="start")))
sub_large <- sub
start(sub_large) <- start(sub_large) - 50
end(sub_large) <- end(sub_large) + 50
over1_table <- table(queryHits(findOverlaps(sub_large, temp)))
sub[as.numeric(names(start_table))]$start_reads <- as.numeric(start_table)
sub[as.numeric(names(over_table))]$over_reads <- as.numeric(over_table)
sub[as.numeric(names(start1_table))]$start_extra_reads <- as.numeric(start1_table)
sub[as.numeric(names(over1_table))]$over_extra_reads <- as.numeric(over1_table)

final <- c(final, sub)

sub <- dmc_gr[which((substr(dmc_gr$seq, 3,5) == "CGA"))]
sub$diff <- dmc[which((substr(dmc_gr$seq, 3,5) == "CGA")),]$diff
sub$pvalue <- dmc[which((substr(dmc_gr$seq, 3,5) == "CGA")),]$pvalue
sub$group <- "TCG"
sub$start_reads <- 0
sub$over_reads <- 0
sub$start_extra_reads <- 0
sub$over_extra_reads <- 0
start(sub) <- start(sub) + 3
end(sub) <- start(sub)
start_table <- table(queryHits(findOverlaps(sub, temp, type="end")))
sub_large <- sub
start(sub_large) <- start(sub_large) - 50
end(sub_large) <- end(sub_large) + 50
over_table <- table(queryHits(findOverlaps(sub_large, temp)))
start(sub) <- start(sub) + 1
end(sub) <- start(sub)
start1_table <- table(queryHits(findOverlaps(sub, temp, type="end")))
sub_large <- sub
start(sub_large) <- start(sub_large) - 50
end(sub_large) <- end(sub_large) + 50
over1_table <- table(queryHits(findOverlaps(sub_large, temp)))
sub[as.numeric(names(start_table))]$start_reads <- as.numeric(start_table)
sub[as.numeric(names(over_table))]$over_reads <- as.numeric(over_table)
sub[as.numeric(names(start1_table))]$start_extra_reads <- as.numeric(start1_table)
sub[as.numeric(names(over1_table))]$over_extra_reads <- as.numeric(over1_table)

final <- c(final, sub)

sub <- dmc_gr[which((substr(dmc_gr$seq, 2,4) == "CCG"))]
sub$diff <- dmc[which((substr(dmc_gr$seq, 2,4) == "CCG")),]$diff
sub$pvalue <- dmc[which((substr(dmc_gr$seq, 2,4) == "CCG")),]$pvalue
sub$group <- "CCG"
sub$start_reads <- 0
sub$over_reads <- 0
sub$start_extra_reads <- 0
sub$over_extra_reads <- 0
start(sub) <- start(sub) + 2
end(sub) <- start(sub)
start_table <- table(queryHits(findOverlaps(sub, temp, type="start")))
sub_large <- sub
start(sub_large) <- start(sub_large) - 50
end(sub_large) <- end(sub_large) + 50
over_table <- table(queryHits(findOverlaps(sub_large, temp)))
start(sub) <- start(sub) - 1
end(sub) <- start(sub)
start1_table <- table(queryHits(findOverlaps(sub, temp, type="start")))
sub_large <- sub
start(sub_large) <- start(sub_large) - 50
end(sub_large) <- end(sub_large) + 50
over1_table <- table(queryHits(findOverlaps(sub_large, temp)))
sub[as.numeric(names(start_table))]$start_reads <- as.numeric(start_table)
sub[as.numeric(names(over_table))]$over_reads <- as.numeric(over_table)
sub[as.numeric(names(start1_table))]$start_extra_reads <- as.numeric(start1_table)
sub[as.numeric(names(over1_table))]$over_extra_reads <- as.numeric(over1_table)

final <- c(final, sub)

sub <- dmc_gr[which((substr(dmc_gr$seq, 3,5) == "CGG"))]
sub$diff <- dmc[which((substr(dmc_gr$seq, 3,5) == "CGG")),]$diff
sub$pvalue <- dmc[which((substr(dmc_gr$seq, 3,5) == "CGG")),]$pvalue
sub$group <- "CCG"
sub$start_reads <- 0
sub$over_reads <- 0
sub$start_extra_reads <- 0
sub$over_extra_reads <- 0
start(sub) <- start(sub) + 3
end(sub) <- start(sub)
start_table <- table(queryHits(findOverlaps(sub, temp, type="end")))
sub_large <- sub
start(sub_large) <- start(sub_large) - 50
end(sub_large) <- end(sub_large) + 50
over_table <- table(queryHits(findOverlaps(sub_large, temp)))
start(sub) <- start(sub) + 1
end(sub) <- start(sub)
start1_table <- table(queryHits(findOverlaps(sub, temp, type="end")))
sub_large <- sub
start(sub_large) <- start(sub_large) - 50
end(sub_large) <- end(sub_large) + 50
over1_table <- table(queryHits(findOverlaps(sub_large, temp)))
sub[as.numeric(names(start_table))]$start_reads <- as.numeric(start_table)
sub[as.numeric(names(over_table))]$over_reads <- as.numeric(over_table)
sub[as.numeric(names(start1_table))]$start_extra_reads <- as.numeric(start1_table)
sub[as.numeric(names(over1_table))]$over_extra_reads <- as.numeric(over1_table)

final <- c(final, sub)

sub <- dmc_gr[which((substr(dmc_gr$seq, 2,4) == "GCG"))]
sub$diff <- dmc[which((substr(dmc_gr$seq, 2,4) == "GCG")),]$diff
sub$pvalue <- dmc[which((substr(dmc_gr$seq, 2,4) == "GCG")),]$pvalue
sub$group <- "GCG"
sub$start_reads <- 0
sub$over_reads <- 0
sub$start_extra_reads <- 0
sub$over_extra_reads <- 0
start(sub) <- start(sub) + 2
end(sub) <- start(sub)
start_table <- table(queryHits(findOverlaps(sub, temp, type="start")))
sub_large <- sub
start(sub_large) <- start(sub_large) - 50
end(sub_large) <- end(sub_large) + 50
over_table <- table(queryHits(findOverlaps(sub_large, temp)))
start(sub) <- start(sub) - 1
end(sub) <- start(sub)
start1_table <- table(queryHits(findOverlaps(sub, temp, type="start")))
sub_large <- sub
start(sub_large) <- start(sub_large) - 50
end(sub_large) <- end(sub_large) + 50
over1_table <- table(queryHits(findOverlaps(sub_large, temp)))
sub[as.numeric(names(start_table))]$start_reads <- as.numeric(start_table)
sub[as.numeric(names(over_table))]$over_reads <- as.numeric(over_table)
sub[as.numeric(names(start1_table))]$start_extra_reads <- as.numeric(start1_table)
sub[as.numeric(names(over1_table))]$over_extra_reads <- as.numeric(over1_table)

final <- c(final, sub)

sub <- dmc_gr[which((substr(dmc_gr$seq, 3,5) == "CGC"))]
sub$diff <- dmc[which((substr(dmc_gr$seq, 3,5) == "CGC")),]$diff
sub$pvalue <- dmc[which((substr(dmc_gr$seq, 3,5) == "CGC")),]$pvalue
sub$group <- "GCG"
sub$start_reads <- 0
sub$over_reads <- 0
sub$start_extra_reads <- 0
sub$over_extra_reads <- 0
start(sub) <- start(sub) + 3
end(sub) <- start(sub)
start_table <- table(queryHits(findOverlaps(sub, temp, type="end")))
sub_large <- sub
start(sub_large) <- start(sub_large) - 50
end(sub_large) <- end(sub_large) + 50
over_table <- table(queryHits(findOverlaps(sub_large, temp)))
start(sub) <- start(sub) + 1
end(sub) <- start(sub)
start1_table <- table(queryHits(findOverlaps(sub, temp, type="end")))
sub_large <- sub
start(sub_large) <- start(sub_large) - 50
end(sub_large) <- end(sub_large) + 50
over1_table <- table(queryHits(findOverlaps(sub_large, temp)))
sub[as.numeric(names(start_table))]$start_reads <- as.numeric(start_table)
sub[as.numeric(names(over_table))]$over_reads <- as.numeric(over_table)
sub[as.numeric(names(start1_table))]$start_extra_reads <- as.numeric(start1_table)
sub[as.numeric(names(over1_table))]$over_extra_reads <- as.numeric(over1_table)

final <- c(final, sub)

sub <- dmc_gr[which(substr(dmc_gr$seq, 1,4) == "ACCG")]
sub$diff <- dmc[which(substr(dmc_gr$seq, 1,4) == "ACCG"),]$diff
sub$pvalue <- dmc[which(substr(dmc_gr$seq, 1,4) == "ACCG"),]$pvalue
sub$group <- "ACCG"
sub$start_reads <- 0
sub$over_reads <- 0
sub$start_extra_reads <- 0
sub$over_extra_reads <- 0
start(sub) <- start(sub) + 1
end(sub) <- start(sub)
start_table <- table(queryHits(findOverlaps(sub, temp, type="start")))
sub_large <- sub
start(sub_large) <- start(sub_large) - 50
end(sub_large) <- end(sub_large) + 50
over_table <- table(queryHits(findOverlaps(sub_large, temp)))
end(sub) <- end(sub) + 1
start(sub) <- end(sub)
start1_table <- table(queryHits(findOverlaps(sub, temp, type="start")))
sub_large <- sub
start(sub_large) <- start(sub_large) - 50
end(sub_large) <- end(sub_large) + 50
over1_table <- table(queryHits(findOverlaps(sub_large, temp)))
sub[as.numeric(names(start_table))]$start_reads <- as.numeric(start_table)
sub[as.numeric(names(over_table))]$over_reads <- as.numeric(over_table)
sub[as.numeric(names(start1_table))]$start_extra_reads <- as.numeric(start1_table)
sub[as.numeric(names(over1_table))]$over_extra_reads <- as.numeric(over1_table)

final <- c(final, sub)

sub <- dmc_gr[which(substr(dmc_gr$seq, 3,6) == "CGGT")]
sub$diff <- dmc[which(substr(dmc_gr$seq, 3,6) == "CGGT"),]$diff
sub$pvalue <- dmc[which(substr(dmc_gr$seq, 3,6) == "CGGT"),]$pvalue
sub$group <- "ACCG"
sub$start_reads <- 0
sub$over_reads <- 0
sub$start_extra_reads <- 0
sub$over_extra_reads <- 0
end(sub) <- end(sub) - 1
start(sub) <- end(sub)
start_table <- table(queryHits(findOverlaps(sub, temp, type="end")))
sub_large <- sub
start(sub_large) <- start(sub_large) - 50
end(sub_large) <- end(sub_large) + 50
over_table <- table(queryHits(findOverlaps(sub_large, temp)))
start(sub) <- start(sub) - 1
end(sub) <- start(sub)
start1_table <- table(queryHits(findOverlaps(sub, temp, type="end")))
sub_large <- sub
start(sub_large) <- start(sub_large) - 50
end(sub_large) <- end(sub_large) + 50
over1_table <- table(queryHits(findOverlaps(sub_large, temp)))
sub[as.numeric(names(start_table))]$start_reads <- as.numeric(start_table)
sub[as.numeric(names(over_table))]$over_reads <- as.numeric(over_table)
sub[as.numeric(names(start1_table))]$start_extra_reads <- as.numeric(start1_table)
sub[as.numeric(names(over1_table))]$over_extra_reads <- as.numeric(over1_table)

final <- c(final, sub)

sub <- dmc_gr[which(substr(dmc_gr$seq, 1,4) == "TCCG")]
sub$diff <- dmc[which(substr(dmc_gr$seq, 1,4) == "TCCG"),]$diff
sub$pvalue <- dmc[which(substr(dmc_gr$seq, 1,4) == "TCCG"),]$pvalue
sub$group <- "TCCG"
sub$start_reads <- 0
sub$over_reads <- 0
sub$start_extra_reads <- 0
sub$over_extra_reads <- 0
start(sub) <- start(sub) + 1
end(sub) <- start(sub)
start_table <- table(queryHits(findOverlaps(sub, temp, type="start")))
sub_large <- sub
start(sub_large) <- start(sub_large) - 50
end(sub_large) <- end(sub_large) + 50
over_table <- table(queryHits(findOverlaps(sub_large, temp)))
end(sub) <- end(sub) + 1
start(sub) <- end(sub)
start1_table <- table(queryHits(findOverlaps(sub, temp, type="start")))
sub_large <- sub
start(sub_large) <- start(sub_large) - 50
end(sub_large) <- end(sub_large) + 50
over1_table <- table(queryHits(findOverlaps(sub_large, temp)))
sub[as.numeric(names(start_table))]$start_reads <- as.numeric(start_table)
sub[as.numeric(names(over_table))]$over_reads <- as.numeric(over_table)
sub[as.numeric(names(start1_table))]$start_extra_reads <- as.numeric(start1_table)
sub[as.numeric(names(over1_table))]$over_extra_reads <- as.numeric(over1_table)

final <- c(final, sub)

sub <- dmc_gr[which(substr(dmc_gr$seq, 3,6) == "CGGA")]
sub$diff <- dmc[which(substr(dmc_gr$seq, 3,6) == "CGGA"),]$diff
sub$pvalue <- dmc[which(substr(dmc_gr$seq, 3,6) == "CGGA"),]$pvalue
sub$group <- "TCCG"
sub$start_reads <- 0
sub$over_reads <- 0
sub$start_extra_reads <- 0
sub$over_extra_reads <- 0
end(sub) <- end(sub) - 1
start(sub) <- end(sub)
start_table <- table(queryHits(findOverlaps(sub, temp, type="end")))
sub_large <- sub
start(sub_large) <- start(sub_large) - 50
end(sub_large) <- end(sub_large) + 50
over_table <- table(queryHits(findOverlaps(sub_large, temp)))
start(sub) <- start(sub) - 1
end(sub) <- start(sub)
start1_table <- table(queryHits(findOverlaps(sub, temp, type="end")))
sub_large <- sub
start(sub_large) <- start(sub_large) - 50
end(sub_large) <- end(sub_large) + 50
over1_table <- table(queryHits(findOverlaps(sub_large, temp)))
sub[as.numeric(names(start_table))]$start_reads <- as.numeric(start_table)
sub[as.numeric(names(over_table))]$over_reads <- as.numeric(over_table)
sub[as.numeric(names(start1_table))]$start_extra_reads <- as.numeric(start1_table)
sub[as.numeric(names(over1_table))]$over_extra_reads <- as.numeric(over1_table)

final <- c(final, sub)

sub <- dmc_gr[which(substr(dmc_gr$seq, 1,4) == "CCCG")]
sub$diff <- dmc[which(substr(dmc_gr$seq, 1,4) == "CCCG"),]$diff
sub$pvalue <- dmc[which(substr(dmc_gr$seq, 1,4) == "CCCG"),]$pvalue
sub$group <- "CCCG"
sub$start_reads <- 0
sub$over_reads <- 0
sub$start_extra_reads <- 0
sub$over_extra_reads <- 0
start(sub) <- start(sub) + 1
end(sub) <- start(sub)
start_table <- table(queryHits(findOverlaps(sub, temp, type="start")))
sub_large <- sub
start(sub_large) <- start(sub_large) - 50
end(sub_large) <- end(sub_large) + 50
over_table <- table(queryHits(findOverlaps(sub_large, temp)))
end(sub) <- end(sub) + 1
start(sub) <- end(sub)
start1_table <- table(queryHits(findOverlaps(sub, temp, type="start")))
sub_large <- sub
start(sub_large) <- start(sub_large) - 50
end(sub_large) <- end(sub_large) + 50
over1_table <- table(queryHits(findOverlaps(sub_large, temp)))
sub[as.numeric(names(start_table))]$start_reads <- as.numeric(start_table)
sub[as.numeric(names(over_table))]$over_reads <- as.numeric(over_table)
sub[as.numeric(names(start1_table))]$start_extra_reads <- as.numeric(start1_table)
sub[as.numeric(names(over1_table))]$over_extra_reads <- as.numeric(over1_table)

final <- c(final, sub)

sub <- dmc_gr[which(substr(dmc_gr$seq, 3,6) == "CGGG")]
sub$diff <- dmc[which(substr(dmc_gr$seq, 3,6) == "CGGG"),]$diff
sub$pvalue <- dmc[which(substr(dmc_gr$seq, 3,6) == "CGGG"),]$pvalue
sub$group <- "CCCG"
sub$start_reads <- 0
sub$over_reads <- 0
sub$start_extra_reads <- 0
sub$over_extra_reads <- 0
end(sub) <- end(sub) - 1
start(sub) <- end(sub)
start_table <- table(queryHits(findOverlaps(sub, temp, type="end")))
sub_large <- sub
start(sub_large) <- start(sub_large) - 50
end(sub_large) <- end(sub_large) + 50
over_table <- table(queryHits(findOverlaps(sub_large, temp)))
start(sub) <- start(sub) - 1
end(sub) <- start(sub)
start1_table <- table(queryHits(findOverlaps(sub, temp, type="end")))
sub_large <- sub
start(sub_large) <- start(sub_large) - 50
end(sub_large) <- end(sub_large) + 50
over1_table <- table(queryHits(findOverlaps(sub_large, temp)))
sub[as.numeric(names(start_table))]$start_reads <- as.numeric(start_table)
sub[as.numeric(names(over_table))]$over_reads <- as.numeric(over_table)
sub[as.numeric(names(start1_table))]$start_extra_reads <- as.numeric(start1_table)
sub[as.numeric(names(over1_table))]$over_extra_reads <- as.numeric(over1_table)

final <- c(final, sub)

sub <- dmc_gr[which(substr(dmc_gr$seq, 1,4) == "GCCG")]
sub$diff <- dmc[which(substr(dmc_gr$seq, 1,4) == "GCCG"),]$diff
sub$pvalue <- dmc[which(substr(dmc_gr$seq, 1,4) == "GCCG"),]$pvalue
sub$group <- "GCCG"
sub$start_reads <- 0
sub$over_reads <- 0
sub$start_extra_reads <- 0
sub$over_extra_reads <- 0
start(sub) <- start(sub) + 1
end(sub) <- start(sub)
start_table <- table(queryHits(findOverlaps(sub, temp, type="start")))
sub_large <- sub
start(sub_large) <- start(sub_large) - 50
end(sub_large) <- end(sub_large) + 50
over_table <- table(queryHits(findOverlaps(sub_large, temp)))
end(sub) <- end(sub) + 1
start(sub) <- end(sub)
start1_table <- table(queryHits(findOverlaps(sub, temp, type="start")))
sub_large <- sub
start(sub_large) <- start(sub_large) - 50
end(sub_large) <- end(sub_large) + 50
over1_table <- table(queryHits(findOverlaps(sub_large, temp)))
sub[as.numeric(names(start_table))]$start_reads <- as.numeric(start_table)
sub[as.numeric(names(over_table))]$over_reads <- as.numeric(over_table)
sub[as.numeric(names(start1_table))]$start_extra_reads <- as.numeric(start1_table)
sub[as.numeric(names(over1_table))]$over_extra_reads <- as.numeric(over1_table)

final <- c(final, sub)

sub <- dmc_gr[which(substr(dmc_gr$seq, 3,6) == "CGGC")]
sub$diff <- dmc[which(substr(dmc_gr$seq, 3,6) == "CGGC"),]$diff
sub$pvalue <- dmc[which(substr(dmc_gr$seq, 3,6) == "CGGC"),]$pvalue
sub$group <- "GCCG"
sub$start_reads <- 0
sub$over_reads <- 0
sub$start_extra_reads <- 0
sub$over_extra_reads <- 0
end(sub) <- end(sub) - 1
start(sub) <- end(sub)
start_table <- table(queryHits(findOverlaps(sub, temp, type="end")))
sub_large <- sub
start(sub_large) <- start(sub_large) - 50
end(sub_large) <- end(sub_large) + 50
over_table <- table(queryHits(findOverlaps(sub_large, temp)))
start(sub) <- start(sub) - 1
end(sub) <- start(sub)
start1_table <- table(queryHits(findOverlaps(sub, temp, type="end")))
sub_large <- sub
start(sub_large) <- start(sub_large) - 50
end(sub_large) <- end(sub_large) + 50
over1_table <- table(queryHits(findOverlaps(sub_large, temp)))
sub[as.numeric(names(start_table))]$start_reads <- as.numeric(start_table)
sub[as.numeric(names(over_table))]$over_reads <- as.numeric(over_table)
sub[as.numeric(names(start1_table))]$start_extra_reads <- as.numeric(start1_table)
sub[as.numeric(names(over1_table))]$over_extra_reads <- as.numeric(over1_table)

final <- c(final, sub)

saveRDS(final, paste0("../../../cfepigenetics_data/Carvalho/12cycles/",basename(file)))

