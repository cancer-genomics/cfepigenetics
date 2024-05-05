library("GenomicRanges")
library("BSgenome.Hsapiens.UCSC.hg19")
library("tidyverse")

args = commandArgs(trailingOnly=TRUE)
case <- as.integer(args[1])

val_files <- readRDS("../../data/selected_validation_files.rds")
luc_files <- readRDS("../../data/selected_lucas_files.rds")
val <- val_files[which(val_files$patient_type == "healthy"),]$pgdx_id
luc_h <- luc_files[which(luc_files$patient_type == 2),]$pgdx_id
h <- c(val, luc_h)
c <- luc_files[which(luc_files$patient_type == 3),]$pgdx_id
all <- c(h,c)

path_methylation <- "../../../cfepigenetics_data/Moss"
meth <- readRDS(file.path(path_methylation, "results", "healthy.rds"))
meth_gr <- GRanges(paste0(meth$chr, ":", meth$pos - 2, "-", meth$pos + 3))
meth$seq <- as.character(getSeq(Hsapiens, meth_gr))
meth$start_ncg <- 0
meth$start_ncg_over <- 0
meth$start_cg <- 0
meth$start_cg_over <- 0
meth$end_cg <- 0
meth$end_cg_over <- 0
meth$end_cgn <- 0
meth$end_cgn_over <- 0

temp <- readRDS(h[case])
meth$beta <- round(rowMeans(meth[,c(4:11)]), digits = 1)
meth <- meth[,c(1,2,3,12,13,14,15,16,17,18,19,20,21)]
meth <- meth[which(meth$chr %in% paste0("chr", c(1:22))),]

tcg <- tibble()
meth_cg <- meth[which(substr(meth$seq, 3, 4) == "CG"),]
t <- as_tibble(findOverlaps(temp, GRanges(paste0(meth_cg$chr, ":", meth_cg$pos, "-", meth_cg$pos + 1))))
t$read_pos <- start(temp[t$queryHits])
t$cg_pos <- meth_cg[t$subjectHits,]$pos
t$beta <- meth_cg[t$subjectHits,]$beta
t$diff <- t$cg_pos - t$read_pos + 1
t1 <- t %>% 
  count(beta, diff)
t$read_pos <- end(temp[t$queryHits])
t$diff <- t$read_pos - t$cg_pos
t2 <- t %>%
  count(beta, diff)
tcg <- rbind(tcg, t1, t2) %>%
  group_by(beta, diff) %>%
  summarise(n = sum(n))

tccg <- tibble()
meth_ccg <- meth[which(substr(meth$seq, 2, 4) == "CCG"),]
t <- as.tibble(findOverlaps(temp, GRanges(paste0(meth_ccg$chr, ":", meth_ccg$pos-1, "-", meth_ccg$pos + 1))))
t$read_pos <- start(temp[t$queryHits])
t$ccg_pos <- meth_ccg[t$subjectHits,]$pos - 1
t$beta <- meth_ccg[t$subjectHits,]$beta
t$diff <- t$ccg_pos - t$read_pos + 1
t1 <- t %>% 
  count(beta, diff)

meth_ccg <- meth[which(substr(meth$seq, 3, 5) == "CGG"),]
t <- as.tibble(findOverlaps(temp, GRanges(paste0(meth_ccg$chr, ":", meth_ccg$pos, "-", meth_ccg$pos + 2))))
t$read_pos <- end(temp[t$queryHits])
t$ccg_pos <- meth_ccg[t$subjectHits,]$pos + 2
t$beta <- meth_ccg[t$subjectHits,]$beta
t$diff <- t$read_pos - t$ccg_pos + 1
t2 <- t %>% 
  count(beta, diff)

tccg <- rbind(tccg, t1, t2) %>%
  group_by(beta, diff) %>%
  summarise(n = sum(n))

tcg$motif <- "CG"
tccg$motif <- "CCG"

t <- rbind(tcg, tccg)

saveRDS(t, paste0("../../../cfepigenetics_data/position_in_read/", basename(h[case])))
