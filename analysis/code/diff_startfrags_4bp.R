library("BSgenome.Hsapiens.UCSC.hg19")
library("GenomicRanges")
args = commandArgs(trailingOnly=TRUE)
case <- as.integer(args[1])
val_files <- readRDS("../../data/selected_validation_files.rds")
luc_files <- readRDS("../../data/selected_lucas_files.rds")
val <- val_files[which(val_files$patient_type == "healthy"),]$pgdx_id
luc_h <- luc_files[which(luc_files$patient_type == 2),]$pgdx_id
h <- c(val, luc_h)
c <- luc_files[which(luc_files$patient_type == 3),]$pgdx_id
he <- readRDS("../../../cfepigenetics_data/Moss/results/healthy.rds")
hegr <- GRanges(paste0(he$chr, ":", he$pos))
end(hegr) <- start(hegr) + 3
start(hegr) <- start(hegr) - 2
he$seq <- as.character(getSeq(Hsapiens, hegr))
he <- he[which(he$chr %in% paste0("chr", c(1:22))),]
hegr <- GRanges(paste0(he$chr, ":", he$pos))
end(hegr) <- start(hegr) + 3
start(hegr) <- start(hegr) - 2
hepos <- he[which(substr(he$seq, 1, 4) %in% c("ACCG", "CCCG", "GCCG", "TCCG")),]
hepos$strand <- "+"
hepos$seq <- substr(hepos$seq, 1, 4)
heneg <- he[which(substr(he$seq, 3, 6) %in% c("CGGT", "CGGG", "CGGC", "CGGA")),]
heneg$strand <- "-"
heneg$seq <- substr(heneg$seq, 3, 6)
heneg$seq <- reverseComplement(DNAStringSet(heneg$seq))
he <- rbind(hepos,heneg)
he$mean <- rowMeans(he[,c(4:11)])
he <- he[which((he$mean <= 0.3)|(he$mean >= 0.7)),]
temp <- readRDS(h[case])
temp <- temp[which((width(temp) >= 100) & (width(temp) <= 220))]
hepos <- he[which(he$strand == "+"),]
heneg <- he[which(he$strand == "-"),]
heposaccg <- hepos[which(hepos$seq == "ACCG"),]
heposcccg <- hepos[which(hepos$seq == "CCCG"),]
heposgccg <- hepos[which(hepos$seq == "GCCG"),]
hepostccg <- hepos[which(hepos$seq == "TCCG"),]
henegaccg <- heneg[which(heneg$seq == "ACCG"),]
henegcccg <- heneg[which(heneg$seq == "CCCG"),]
heneggccg <- heneg[which(heneg$seq == "GCCG"),]
henegtccg <- heneg[which(heneg$seq == "TCCG"),]
m <- matrix(0,16,51)
for (i in seq(-26,24,by=1)) {
  print(i)
  pos_accg_m <- length(unique(queryHits(findOverlaps(temp, GRanges(paste0(heposaccg[which(heposaccg$mean > 0.5),]$chr, ":", heposaccg[which(heposaccg$mean > 0.5),]$pos + i)), type="start"))))
  pos_cccg_m <- length(unique(queryHits(findOverlaps(temp, GRanges(paste0(heposcccg[which(heposcccg$mean > 0.5),]$chr, ":", heposcccg[which(heposcccg$mean > 0.5),]$pos + i)), type="start"))))
  pos_gccg_m <- length(unique(queryHits(findOverlaps(temp, GRanges(paste0(heposgccg[which(heposgccg$mean > 0.5),]$chr, ":", heposgccg[which(heposgccg$mean > 0.5),]$pos + i)), type="start"))))
  pos_tccg_m <- length(unique(queryHits(findOverlaps(temp, GRanges(paste0(hepostccg[which(hepostccg$mean > 0.5),]$chr, ":", hepostccg[which(hepostccg$mean > 0.5),]$pos + i)), type="start"))))
  neg_accg_m <- length(unique(queryHits(findOverlaps(temp, GRanges(paste0(henegaccg[which(henegaccg$mean > 0.5),]$chr, ":", henegaccg[which(henegaccg$mean > 0.5),]$pos - i + 1)), type="end"))))
  neg_cccg_m <- length(unique(queryHits(findOverlaps(temp, GRanges(paste0(henegcccg[which(henegcccg$mean > 0.5),]$chr, ":", henegcccg[which(henegcccg$mean > 0.5),]$pos - i + 1)), type="end"))))
  neg_gccg_m <- length(unique(queryHits(findOverlaps(temp, GRanges(paste0(heneggccg[which(heneggccg$mean > 0.5),]$chr, ":", heneggccg[which(heneggccg$mean > 0.5),]$pos - i + 1)), type="end"))))
  neg_tccg_m <- length(unique(queryHits(findOverlaps(temp, GRanges(paste0(henegtccg[which(henegtccg$mean > 0.5),]$chr, ":", henegtccg[which(henegtccg$mean > 0.5),]$pos - i + 1)), type="end"))))
  
  pos_accg_u <- length(unique(queryHits(findOverlaps(temp, GRanges(paste0(heposaccg[which(heposaccg$mean < 0.5),]$chr, ":", heposaccg[which(heposaccg$mean < 0.5),]$pos + i)), type="start"))))
  pos_cccg_u <- length(unique(queryHits(findOverlaps(temp, GRanges(paste0(heposcccg[which(heposcccg$mean < 0.5),]$chr, ":", heposcccg[which(heposcccg$mean < 0.5),]$pos + i)), type="start"))))
  pos_gccg_u <- length(unique(queryHits(findOverlaps(temp, GRanges(paste0(heposgccg[which(heposgccg$mean < 0.5),]$chr, ":", heposgccg[which(heposgccg$mean < 0.5),]$pos + i)), type="start"))))
  pos_tccg_u <- length(unique(queryHits(findOverlaps(temp, GRanges(paste0(hepostccg[which(hepostccg$mean < 0.5),]$chr, ":", hepostccg[which(hepostccg$mean < 0.5),]$pos + i)), type="start"))))
  neg_accg_u <- length(unique(queryHits(findOverlaps(temp, GRanges(paste0(henegaccg[which(henegaccg$mean < 0.5),]$chr, ":", henegaccg[which(henegaccg$mean < 0.5),]$pos - i + 1)), type="end"))))
  neg_cccg_u <- length(unique(queryHits(findOverlaps(temp, GRanges(paste0(henegcccg[which(henegcccg$mean < 0.5),]$chr, ":", henegcccg[which(henegcccg$mean < 0.5),]$pos - i + 1)), type="end"))))
  neg_gccg_u <- length(unique(queryHits(findOverlaps(temp, GRanges(paste0(heneggccg[which(heneggccg$mean < 0.5),]$chr, ":", heneggccg[which(heneggccg$mean < 0.5),]$pos - i + 1)), type="end"))))
  neg_tccg_u <- length(unique(queryHits(findOverlaps(temp, GRanges(paste0(henegtccg[which(henegtccg$mean < 0.5),]$chr, ":", henegtccg[which(henegtccg$mean < 0.5),]$pos - i + 1)), type="end"))))
  
  m[,(i+27)] <- c(pos_accg_m, neg_accg_m, pos_cccg_m, neg_cccg_m, pos_gccg_m, neg_gccg_m, pos_tccg_m, neg_tccg_m, 
                  pos_accg_u, neg_accg_u, pos_cccg_u, neg_cccg_u, pos_gccg_u, neg_gccg_u, pos_tccg_u, neg_tccg_u)
}
saveRDS(m, paste0("../../../cfepigenetics_data/diff_startfrags_4bp/", basename(h[case])))
