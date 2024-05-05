library("GenomicRanges")
library("BSgenome.Hsapiens.UCSC.hg19")

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
meth$beta <- rowMeans(meth[,c(4:11)])
meth <- meth[,c(1,2,3,12,13,14,15,16,17,18,19,20,21)]

t <- table(subjectHits(findOverlaps(temp, GRanges(paste0(meth$chr, ":", meth$pos - 1)), type = "start")))
meth[as.integer(names(t)),]$start_ncg <- as.integer(t)

t <- table(subjectHits(findOverlaps(temp, GRanges(paste0(meth$chr, ":", meth$pos - 51, "-", meth$pos + 49)))))
meth[as.integer(names(t)),]$start_ncg_over <- as.integer(t)

t <- table(subjectHits(findOverlaps(temp, GRanges(paste0(meth$chr, ":", meth$pos)), type = "start")))
meth[as.integer(names(t)),]$start_cg <- as.integer(t)

t <- table(subjectHits(findOverlaps(temp, GRanges(paste0(meth$chr, ":", meth$pos - 50, "-", meth$pos + 50)))))
meth[as.integer(names(t)),]$start_cg_over <- as.integer(t)

t <- table(subjectHits(findOverlaps(temp, GRanges(paste0(meth$chr, ":", meth$pos + 1)), type = "end")))
meth[as.integer(names(t)),]$end_cg <- as.integer(t)

t <- table(subjectHits(findOverlaps(temp, GRanges(paste0(meth$chr, ":", meth$pos - 49, "-", meth$pos + 51)))))
meth[as.integer(names(t)),]$end_cg_over <- as.integer(t)

t <- table(subjectHits(findOverlaps(temp, GRanges(paste0(meth$chr, ":", meth$pos + 2)), type = "end")))
meth[as.integer(names(t)),]$end_cgn <- as.integer(t)

t <- table(subjectHits(findOverlaps(temp, GRanges(paste0(meth$chr, ":", meth$pos - 48, "-", meth$pos + 52)))))
meth[as.integer(names(t)),]$end_cgn_over <- as.integer(t)

saveRDS(meth, paste0("../../../cfepigenetics_data/Moss_cases/", basename(h[case])))
