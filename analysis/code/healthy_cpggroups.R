library("tidyverse")
library("BSgenome.Hsapiens.UCSC.hg19")
library("GenomicRanges")
library("annotatr")

args = commandArgs(trailingOnly=TRUE)
case <- as.integer(args[1])

val_files <- readRDS("../../data/selected_validation_files.rds")
luc_files <- readRDS("../../data/selected_lucas_files.rds")
val <- val_files[which(val_files$patient_type == "healthy"),]$pgdx_id
luc_h <- luc_files[which(luc_files$patient_type == 2),]$pgdx_id
h <- c(val, luc_h)
c <- luc_files[which(luc_files$patient_type == 3),]$pgdx_id
normal <- readRDS("../../../cfepigenetics_data/Moss/results/healthy.rds")
ngr <- GRanges(paste0(normal$chr, ":", normal$pos, "-", normal$pos))
files <- list.files("../../../cfepigenetics_data/Moss_cases", full.names=T)

annots = c('hg19_cpgs')
annotations = build_annotations(genome = 'hg19', annotations = annots)
normal$annot <- "NAp"
over <- findOverlaps(ngr, annotations)
normal[queryHits(over),]$annot <- annotations[subjectHits(over)]$type
start(ngr) <- start(ngr) - 2
end(ngr) <- start(ngr) + 5
normal$seq <- as.character(getSeq(Hsapiens, ngr))
normal$beta <- round(rowMeans(normal[,c(4:11)]), digits = 1)
motifs3_fw <- c("ACG", "CCG", "GCG", "TCG")
motifs3_re <- c("CGT", "CGG", "CGC", "CGA")
motifs4_fw <- c("ACCG", "CCCG", "GCCG", "TCCG")
motifs4_re <- c("CGGT", "CGGG", "CGGC", "CGGA")
t <- tibble()
i <- h[case]
print(paste0(basename(i), " : ", match(i, h), " / ", length(h)))
temp <- readRDS(files[match(basename(i), basename(files))]) 
for (chr in unique(normal$annot)) {
  for (b in seq(0,1,0.1)) {
    print(b)
    for (motif in motifs3_fw) {
      start <- sum(temp[which((normal$annot == chr) & (normal$beta == as.numeric(as.character(b))) & (substr(normal$seq, 2, 4) == motif) & (normal$chr != "chrX") & (normal$chr != "chrY")),]$start_cg)
      end <- sum(temp[which((normal$annot == chr) & (normal$beta == as.numeric(as.character(b))) & (substr(normal$seq, 3, 5) == motifs3_re[match(motif, motifs3_fw)])  & (normal$chr != "chrX") & (normal$chr != "chrY")),]$end_cg)
      start_over <- sum(temp[which((normal$annot == chr) & (normal$beta == as.numeric(as.character(b))) & (substr(normal$seq, 2, 4) == motif)  & (normal$chr != "chrX") & (normal$chr != "chrY")),]$start_cg_over)
      end_over <- sum(temp[which((normal$annot == chr) & (normal$beta == as.numeric(as.character(b))) & (substr(normal$seq, 3, 5) == motifs3_re[match(motif, motifs3_fw)])  & (normal$chr != "chrX") & (normal$chr != "chrY")),]$end_cg_over)
      t_temp <- tibble("case" = sub(".rds", "", basename(i)),
                       "group" = chr,
                       "beta" = b,
                       "motif" = motif,
                       "start" = start + end,
                       "over" = start_over + end_over,
                       "ratio" = (start + end) / (start_over + end_over))
      t <- rbind(t, t_temp)
    }
    for (motif in motifs4_fw) {
      start <- sum(temp[which((normal$annot == chr) & (normal$beta == as.numeric(as.character(b))) & (substr(normal$seq, 1, 4) == motif) & (normal$chr != "chrX") & (normal$chr != "chrY")),]$start_ncg)
      end <- sum(temp[which((normal$annot == chr) & (normal$beta == as.numeric(as.character(b))) & (substr(normal$seq, 3, 6) == motifs4_re[match(motif, motifs4_fw)])  & (normal$chr != "chrX") & (normal$chr != "chrY")),]$end_cgn)
      start_over <- sum(temp[which((normal$annot == chr) & (normal$beta == as.numeric(as.character(b))) & (substr(normal$seq, 1, 4) == motif)  & (normal$chr != "chrX") & (normal$chr != "chrY")),]$start_ncg_over)
      end_over <- sum(temp[which((normal$annot == chr) & (normal$beta == as.numeric(as.character(b))) & (substr(normal$seq, 3, 6) == motifs4_re[match(motif, motifs4_fw)])  & (normal$chr != "chrX") & (normal$chr != "chrY")),]$end_cgn_over)
      t_temp <- tibble("case" = sub(".rds", "", basename(i)),
                       "group" = chr,
                       "beta" = b,
                       "motif" = motif,
                       "start" = start + end,
                       "over" = start_over + end_over,
                       "ratio" = (start + end) / (start_over + end_over))
      t <- rbind(t, t_temp)
    }
  }
}
saveRDS(t, paste0("../../../cfepigenetics_data/healthy_cpggroups/", basename(i)))
