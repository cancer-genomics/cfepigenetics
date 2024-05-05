library("GenomicRanges")
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

temp <- readRDS(h[case])
t <- tibble("total" = 0,
            "less_100" = 0,
            "more_100_less_220" = 0,
            "more_220" =0,
            "somatic_all" = 0,
            "somatic_small" = 0,
            "sex_all" = 0,
            "sex_small" = 0)
t$total <- length(temp)
t$less_100 <- length(temp[which(width(temp) < 100)])
t$more_100_less_220 <- length(temp[which((width(temp) >= 100) & (width(temp) <= 220))])
t$more_220 <- length(temp[which(width(temp) > 220)])
t$somatic_all <- length(temp[which(seqnames(temp) %in% paste0("chr", c(1:22)))])
t$somatic_small <- length(temp[which((seqnames(temp) %in% paste0("chr", c(1:22))) & (width(temp) >= 100) & (width(temp) <= 220))])
t$sex_all <- length(temp[which(seqnames(temp) == "chrX")])
t$sex_small <- length(temp[which((seqnames(temp) == "chrX") & (width(temp) >= 100) & (width(temp) <= 220))])


saveRDS(t, paste0("../../../cfepigenetics_data/read_counts/", basename(h[case])))
