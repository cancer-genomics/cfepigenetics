library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg19)
library(tidyverse)

args = commandArgs(trailingOnly=TRUE)
case <- as.integer(args[1])

li <- list.files("/dcs04/scharpf/data/annapragada/Liver_Curated/granges", full.names = T)
va <- list.files("/dcl01/scharpf1/data/cristiano/projects/delfi-validation/granges/", full.names = T)
lu <- list.files("/dcl02/leased/cglab/rscharpf/cristiano/projects/lucas/granges/", full.names = T)
ci <- list.files("/dcl01/scharpf1/data/cristiano/projects/DELFI_paper/granges", full.names = T)

ls <- c(li, va, lu, ci)

dmc_gr <- readRDS("/dcl02/leased/cglab/rscharpf/mnoe/PROJECTS/cfepigenetics/data/DMR_cfdna_fallo.rds")

file <- ls[case]
temp <- readRDS(file)
final <- GRanges()

sub <- dmc_gr[which((substr(dmc_gr$seq, 2,4) == "ACG"))]
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

temp <- as_tibble(final)
temp <- temp %>%
  group_by(dir, group) %>%
  summarise(start_reads = sum(start_reads), over_reads = sum(over_reads), start_reads_extra = sum(start_extra_reads), over_reads_extra = sum(over_extra_reads))
temp$ratio <- temp$start_reads / temp$over_reads
temp$ratio_extra <- temp$start_reads_extra / temp$over_reads_extra
temp$id <- basename(file)

saveRDS(temp, paste0("/dcs04/scharpf/data/mnoe/DMR_Loyfer/fallo/",basename(file)))

