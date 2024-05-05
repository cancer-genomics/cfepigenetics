library("tidyverse")

args = commandArgs(trailingOnly=TRUE)
case <- as.integer(args[1])

ls <- list.files("../../../cfepigenetics_data/reads_per_tilesmall", full.names=T)
done <- list.files("../../../cfepigenetics_data/cfdna_ratio_motif", full.names=F)
ls <- ls[which(!basename(ls) %in% done)]
file <- ls[case]

if ((length(grep("2023-04-06", file.info(file)$ctime)) == 1) |
    (length(grep("2023-04-07", file.info(file)$ctime)) == 1) |
    (length(grep("2023-04-10", file.info(file)$ctime)) == 1)) {
    print(file)
    temp <- readRDS(file)
    ratio1 <- round(temp$start_reads / temp$overlap, digits = 4)
    ratio1 <- ratio1[which(temp$overlap >= 200)]
    ratio2 <- round(temp$end_reads / temp$overlap, digits = 4)
    ratio2 <- ratio2[which(temp$overlap >= 200)]
    seq1 <- substr(as.character(temp$seq), 2, 5)
    seq1 <- seq1[which(temp$overlap >= 200)]
    seq2 <- substr(as.character(reverseComplement(temp$seq)), 2, 5)
    seq2 <- seq2[which(temp$overlap >= 200)]
    temp_4df <- tibble("motif" = c(seq1, seq2),
                       "ratio" = c(ratio1, ratio2),
                       "count" = 1)
    temp_4df <- temp_4df %>%
      group_by(motif, ratio) %>%
      summarize(num = sum(count))
    seq1 <- substr(as.character(temp$seq), 2, 4)
    seq1 <- seq1[which(temp$overlap >= 200)]
    seq2 <- substr(as.character(reverseComplement(temp$seq)), 2, 4)
    seq2 <- seq2[which(temp$overlap >= 200)]
    temp_3df <- tibble("motif" = c(seq1, seq2),
                       "ratio" = c(ratio1, ratio2),
                       "count" = 1)
    temp_3df <- temp_3df %>%
      group_by(motif, ratio) %>%
      summarize(num = sum(count))
    df <- rbind(temp_3df, temp_4df)
    saveRDS(df, paste0("../../../cfepigenetics_data/cfdna_ratio_motif/", basename(file)))
  }
