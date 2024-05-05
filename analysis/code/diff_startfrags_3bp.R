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
hepos <- he[which(substr(he$seq, 2, 4) %in% c("ACG", "CCG", "GCG", "TCG")),]
hepos$strand <- "+"
hepos$seq <- substr(hepos$seq, 2, 4)
heneg <- he[which(substr(he$seq, 3, 5) %in% c("CGT", "CGG", "CGC", "CGA")),]
heneg$strand <- "-"
heneg$seq <- substr(heneg$seq, 3, 5)
heneg$seq <- reverseComplement(DNAStringSet(heneg$seq))
he <- rbind(hepos,heneg)
he$mean <- rowMeans(he[,c(4:11)])
he <- he[which((he$mean <= 0.3)|(he$mean >= 0.7)),]
temp <- readRDS(h[case])
temp <- temp[which((width(temp) >= 100) & (width(temp) <= 220))]
hepos <- he[which(he$strand == "+"),]
heneg <- he[which(he$strand == "-"),]
heposacg <- hepos[which(hepos$seq == "ACG"),]
heposccg <- hepos[which(hepos$seq == "CCG"),]
heposgcg <- hepos[which(hepos$seq == "GCG"),]
hepostcg <- hepos[which(hepos$seq == "TCG"),]
henegacg <- heneg[which(heneg$seq == "ACG"),]
henegccg <- heneg[which(heneg$seq == "CCG"),]
heneggcg <- heneg[which(heneg$seq == "GCG"),]
henegtcg <- heneg[which(heneg$seq == "TCG"),]
m <- matrix(0,16,51)
for (i in seq(-26,24,by=1)) {
  print(i)
  pos_acg_m <- length(unique(queryHits(findOverlaps(temp, GRanges(paste0(heposacg[which(heposacg$mean > 0.5),]$chr, ":", heposacg[which(heposacg$mean > 0.5),]$pos + i)), type="start"))))
  pos_ccg_m <- length(unique(queryHits(findOverlaps(temp, GRanges(paste0(heposccg[which(heposccg$mean > 0.5),]$chr, ":", heposccg[which(heposccg$mean > 0.5),]$pos + i)), type="start"))))
  pos_gcg_m <- length(unique(queryHits(findOverlaps(temp, GRanges(paste0(heposgcg[which(heposgcg$mean > 0.5),]$chr, ":", heposgcg[which(heposgcg$mean > 0.5),]$pos + i)), type="start"))))
  pos_tcg_m <- length(unique(queryHits(findOverlaps(temp, GRanges(paste0(hepostcg[which(hepostcg$mean > 0.5),]$chr, ":", hepostcg[which(hepostcg$mean > 0.5),]$pos + i)), type="start"))))
  neg_acg_m <- length(unique(queryHits(findOverlaps(temp, GRanges(paste0(henegacg[which(henegacg$mean > 0.5),]$chr, ":", henegacg[which(henegacg$mean > 0.5),]$pos - i + 1)), type="end"))))
  neg_ccg_m <- length(unique(queryHits(findOverlaps(temp, GRanges(paste0(henegccg[which(henegccg$mean > 0.5),]$chr, ":", henegccg[which(henegccg$mean > 0.5),]$pos - i + 1)), type="end"))))
  neg_gcg_m <- length(unique(queryHits(findOverlaps(temp, GRanges(paste0(heneggcg[which(heneggcg$mean > 0.5),]$chr, ":", heneggcg[which(heneggcg$mean > 0.5),]$pos - i + 1)), type="end"))))
  neg_tcg_m <- length(unique(queryHits(findOverlaps(temp, GRanges(paste0(henegtcg[which(henegtcg$mean > 0.5),]$chr, ":", henegtcg[which(henegtcg$mean > 0.5),]$pos - i + 1)), type="end"))))
  
  pos_acg_u <- length(unique(queryHits(findOverlaps(temp, GRanges(paste0(heposacg[which(heposacg$mean < 0.5),]$chr, ":", heposacg[which(heposacg$mean < 0.5),]$pos + i)), type="start"))))
  pos_ccg_u <- length(unique(queryHits(findOverlaps(temp, GRanges(paste0(heposccg[which(heposccg$mean < 0.5),]$chr, ":", heposccg[which(heposccg$mean < 0.5),]$pos + i)), type="start"))))
  pos_gcg_u <- length(unique(queryHits(findOverlaps(temp, GRanges(paste0(heposgcg[which(heposgcg$mean < 0.5),]$chr, ":", heposgcg[which(heposgcg$mean < 0.5),]$pos + i)), type="start"))))
  pos_tcg_u <- length(unique(queryHits(findOverlaps(temp, GRanges(paste0(hepostcg[which(hepostcg$mean < 0.5),]$chr, ":", hepostcg[which(hepostcg$mean < 0.5),]$pos + i)), type="start"))))
  neg_acg_u <- length(unique(queryHits(findOverlaps(temp, GRanges(paste0(henegacg[which(henegacg$mean < 0.5),]$chr, ":", henegacg[which(henegacg$mean < 0.5),]$pos - i + 1)), type="end"))))
  neg_ccg_u <- length(unique(queryHits(findOverlaps(temp, GRanges(paste0(henegccg[which(henegccg$mean < 0.5),]$chr, ":", henegccg[which(henegccg$mean < 0.5),]$pos - i + 1)), type="end"))))
  neg_gcg_u <- length(unique(queryHits(findOverlaps(temp, GRanges(paste0(heneggcg[which(heneggcg$mean < 0.5),]$chr, ":", heneggcg[which(heneggcg$mean < 0.5),]$pos - i + 1)), type="end"))))
  neg_tcg_u <- length(unique(queryHits(findOverlaps(temp, GRanges(paste0(henegtcg[which(henegtcg$mean < 0.5),]$chr, ":", henegtcg[which(henegtcg$mean < 0.5),]$pos - i + 1)), type="end"))))
  
  m[,(i+27)] <- c(pos_acg_m, neg_acg_m, pos_ccg_m, neg_ccg_m, pos_gcg_m, neg_gcg_m, pos_tcg_m, neg_tcg_m, 
                  pos_acg_u, neg_acg_u, pos_ccg_u, neg_ccg_u, pos_gcg_u, neg_gcg_u, pos_tcg_u, neg_tcg_u)
}
saveRDS(m, paste0("../../../cfepigenetics_data/diff_startfrags_3bp/", basename(h[case])))
