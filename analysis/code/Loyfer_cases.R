library("GenomicRanges")
library("BSgenome.Hsapiens.UCSC.hg19")
library("tidyverse")
library("Biostrings")

args = commandArgs(trailingOnly=TRUE)
case <- as.integer(args[1])

val_files <- readRDS("../../data/selected_validation_files.rds")
luc_files <- readRDS("../../data/selected_lucas_files.rds")
val <- val_files[which(val_files$patient_type == "healthy"),]$pgdx_id
luc_h <- luc_files[which(luc_files$patient_type == 2),]$pgdx_id
h <- c(val, luc_h)

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
m <- matrix(0, nrow=length(cg.loc), ncol=8)
temp <- readRDS(h[case])
temp <- temp[which((width(temp) <=220) & (width(temp) >=100))]

end(cg.loc) <- start(cg.loc)
otemp <- as_tibble(as.data.frame(findOverlaps(cg.loc, temp, type = "start")))
otemp$n <- 1
otemp <- otemp %>%
  group_by(queryHits) %>%
  summarise(start = sum(n))
m[otemp$queryHits,1] <- otemp$start

cg.loc_large <- cg.loc
start(cg.loc_large) <- start(cg.loc_large) - 50
end(cg.loc_large) <- end(cg.loc_large) + 50
otemp <- as_tibble(as.data.frame(findOverlaps(cg.loc_large, temp)))
otemp$n <- 1
otemp <- otemp %>%
  group_by(queryHits) %>%
  summarise(start = sum(n))
m[otemp$queryHits,2] <- otemp$start

start(cg.loc) <- start(cg.loc) - 1
end(cg.loc) <- start(cg.loc)
otemp <- as_tibble(as.data.frame(findOverlaps(cg.loc, temp, type = "start")))
otemp$n <- 1
otemp <- otemp %>%
  group_by(queryHits) %>%
  summarise(start = sum(n))
m[otemp$queryHits,3] <- otemp$start

cg.loc_large <- cg.loc
start(cg.loc_large) <- start(cg.loc_large) - 50
end(cg.loc_large) <- end(cg.loc_large) + 50
otemp <- as_tibble(as.data.frame(findOverlaps(cg.loc_large, temp)))
otemp$n <- 1
otemp <- otemp %>%
  group_by(queryHits) %>%
  summarise(start = sum(n))
m[otemp$queryHits,4] <- otemp$start

end(cg.loc) <- end(cg.loc) + 2
start(cg.loc) <- end(cg.loc)
otemp <- as_tibble(as.data.frame(findOverlaps(cg.loc, temp, type = "end")))
otemp$n <- 1
otemp <- otemp %>%
  group_by(queryHits) %>%
  summarise(start = sum(n))
m[otemp$queryHits,5] <- otemp$start

cg.loc_large <- cg.loc
start(cg.loc_large) <- start(cg.loc_large) - 50
end(cg.loc_large) <- end(cg.loc_large) + 50
otemp <- as_tibble(as.data.frame(findOverlaps(cg.loc_large, temp)))
otemp$n <- 1
otemp <- otemp %>%
  group_by(queryHits) %>%
  summarise(start = sum(n))
m[otemp$queryHits,6] <- otemp$start

end(cg.loc) <- end(cg.loc) + 1
end(cg.loc) <- start(cg.loc)
otemp <- as_tibble(as.data.frame(findOverlaps(cg.loc, temp, type = "end")))
otemp$n <- 1
otemp <- otemp %>%
  group_by(queryHits) %>%
  summarise(start = sum(n))
m[otemp$queryHits,7] <- otemp$start

cg.loc_large <- cg.loc
start(cg.loc_large) <- start(cg.loc_large) - 50
end(cg.loc_large) <- end(cg.loc_large) + 50
otemp <- as_tibble(as.data.frame(findOverlaps(cg.loc_large, temp)))
otemp$n <- 1
otemp <- otemp %>%
  group_by(queryHits) %>%
  summarise(start = sum(n))
m[otemp$queryHits,8] <- otemp$start

saveRDS(m, paste0("../../../cfepigenetics_data/Loyfer_cases/", basename(h[case])))



