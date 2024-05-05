library("tidyverse")
library("readr")
library("GenomicRanges")
library("BSgenome.Hsapiens.UCSC.hg19")
library("biomaRt")
library("annotatr")
library("readxl")
library("cowplot")
library("ggrepel")

loc_n <- '/dcs04/scharpf/data/nniknafs/delfi-brain/data/xenograft'
meta <- read_tsv(file.path(loc_n, "meta", "metadata.tsv"))
meta_all <- meta[!is.na(meta$idh),]$pgdx.id
meta_fl <- meta[which(meta$loc == "flank"),]$pgdx.id

loc_c <- '/dcs05/scharpf/data/ccherry/2023.12_mouse_multiomics'
#met <- read_csv(file.path(loc_c, "methylation", "genotype_all_samples.csv"))
met <- read_csv(file.path(loc_c, "methylation", "genotype_flank_samples.csv"))
rna <- read_csv(file.path(loc_c, "rna", "flank_genoptye_de.csv"))

loc_l <- '/dcl02/leased/cglab/rscharpf/mnoe/PROJECTS/cfepigenetics/data'
tra <- read_tsv(file.path(loc_l, "transcriptAnno-GRCh37.75.tsv"), col_names=F)

#get CpG-islands, -shores, -shelves and open sea
annots = c('hg19_cpgs')
annotations = build_annotations(genome = 'hg19', annotations = annots)
#prep differential methylation data
met <- met[order(met$stat),]
met_gr <- GRanges(seqnames=met$seqnames, ranges=IRanges(start=met$start, end=met$end + 1))
met$seq <- as.character(getSeq(Hsapiens, met_gr))
met <- met[which(met$seq == "CG"),]
met <- met[which(met$seqnames %in% paste0("chr", c(1:22))),]
met_gr <- GRanges(seqnames=met$seqnames, ranges=IRanges(start=met$start, end=met$end))
o <- findOverlaps(met_gr, annotations)
met$group <- "NA"
met[queryHits(o),]$group <- annotations[subjectHits(o),]$type
met$group_id <- 0
met[queryHits(o),]$group_id <- subjectHits(o)
met_cpg <- met %>%
  group_by(group_id) %>%
  summarise(stat_mean = mean(stat))
ann_cpg <- annotations[met_cpg$group_id]
ann_cpg$stat_mean <- met_cpg$stat_mean
ann_cpg <- ann_cpg[order(ann_cpg$stat_mean)]
ann_cpg$middle <- round((start(ann_cpg) + end(ann_cpg))/2)
start(ann_cpg) <- ann_cpg$middle - 250
end(ann_cpg) <- ann_cpg$middle + 250
ann_cpg_is <- ann_cpg[which(ann_cpg$type == "hg19_cpg_islands")]

ensembl = useEnsembl(biomart="ensembl")
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
searchDatasets(mart = mart, pattern = "hsapiens")
x <- as_tibble(getBM(attributes = c('ensembl_gene_id', 'hgnc_symbol'), mart = mart))
tra$gene <- x[match(tra$X1, x$ensembl_gene_id),]$hgnc_symbol
tra$tss <- 0
tra <- tra[which(tra$X2 %in% c(1:22)),]
tra[which(tra$X5 == "+"),]$tss <- tra[which(tra$X5 == "+"),]$X3
tra[which(tra$X5 == "-"),]$tss <- tra[which(tra$X5 == "-"),]$X4
rna <- rna[!is.na(rna$stat),]
rna <- rna[order(rna$stat),]
rna$ensembl <- "NA"
rna$ensembl <- x[match(rna$...1, x$hgnc_symbol),]$ensembl_gene_id
rna$seqnames <- "NA"
rna$seqnames <- paste0("chr", tra[match(rna$ensembl, tra$X1),]$X2)
rna[which(rna$seqnames == "chrNA"),]$seqnames <- "NA"
rna$tss <- "NA"
rna$tss <- tra[match(rna$ensembl, tra$X1),]$tss
rna <- rna[which(!is.na(rna$tss)),]
rna_gr <- GRanges(seqnames=rna$seqnames, ranges = IRanges(start = rna$tss - 250, end = rna$tss + 250))

loc_m <- '/dcl02/leased/cglab/rscharpf/mnoe/PROJECTS/mouse/data'
idh <- read_excel(file.path(loc_m, "PMID31292202_SuppTable2.xlsx"), col_names = T, skip = 1)
idh$KO_to_Par.logFC <- as.numeric(idh$KO_to_Par.logFC)
idh$KO_to_Par.P.Value <- as.numeric(idh$KO_to_Par.P.Value)
idh$KO_to_Par.t <- as.numeric(idh$KO_to_Par.t)
idh$KO_to_Par.P.Value <- idh$KO_to_Par.P.Value * sign(idh$KO_to_Par.logFC)
idh <- idh[order(-idh$KO_to_Par.t),]
idh$seqnames <- "NA"
q <- readRDS("/dcl02/leased/cglab/rscharpf/mnoe/PROJECTS/cfepigenetics_data/Moss/results/healthy.rds")
idh$seqnames <- q[match(idh$CpGsite, q$Row.names),]$chr
idh$pos <- 0
idh$pos <- q[match(idh$CpGsite, q$Row.names),]$pos
idh_gr <- GRanges(seqnames = idh$seqnames, ranges = IRanges(start = idh$pos, end = idh$pos))


n <- c(1000)
t <- tibble()
for (m in n) {
  print(m)
  rna_gr_min <- rna_gr[c(1:m)]
  rna_gr_plus <- rna_gr[c((length(rna_gr)-(m-1)):length(rna_gr))]
  met_gr_min <- met_gr[c(1:m)]
  met_gr_plus <- met_gr[c((length(met_gr)-(m-1)):length(met_gr))]
  met_gr_is_min <- met_gr[which(met$group == "hg19_cpg_islands")][c(1:m)]
  met_gr_is_plus <- met_gr[which(met$group == "hg19_cpg_islands")][c((length(met_gr[which(met$group == "hg19_cpg_islands")])-(m-1)):length(met_gr[which(met$group == "hg19_cpg_islands")]))]
  ann_cpg_is_min <- ann_cpg_is[c(1:m)]
  ann_cpg_is_plus <- ann_cpg_is[c((length(ann_cpg_is)-(m-1)):length(ann_cpg_is))]
  idh_gr_min <- idh_gr[which(idh$KO_to_Par.t > 0),][c(1:m)]
  if (m > length(idh_gr[which(idh$KO_to_Par.t < 0),])) {
    idh_gr_plus <- idh_gr[which(idh$KO_to_Par.t < 0),]
  } else {
    idh_gr_plus <- idh_gr[which(idh$KO_to_Par.t < 0),][c((length(idh_gr[which(idh$KO_to_Par.t < 0),])-(m-1)):length(idh_gr[which(idh$KO_to_Par.t < 0),]))]
  }
    for (i in meta_all) {
    print(i)
    temp <- readRDS(file.path(loc_n, "granges", paste0(i,".graft.rds")))
    temp <- temp[which((width(temp) >=100) & (width(temp) <= 220))]
    ann_cpg_is_min_length <- sum(width(ann_cpg_is_min))
    ann_cpg_is_plus_length = sum(width(ann_cpg_is_plus))
    temp_tb <- tibble("id" = i,
                      "mut" = meta[match(i, meta$pgdx.id),]$idh,
                      "inplant" = meta[match(i, meta$pgdx.id),]$loc,
                      "top" = m,
                      "rna_min" = length(queryHits(findOverlaps(temp, rna_gr_min))),
                      "rna_min_length" = sum(width(rna_gr_min)),
                      "rna_plus" = length(queryHits(findOverlaps(temp, rna_gr_plus))),
                      "rna_plus_length" = sum(width(rna_gr_plus)),
                      "met_min" = length(queryHits(findOverlaps(temp, met_gr_min))),
                      "met_min_length" = sum(width(met_gr_min)),
                      "met_plus" = length(queryHits(findOverlaps(temp, met_gr_plus))),
                      "met_plus_length" =  sum(width(met_gr_plus)),
                      "met_is_min" = length(queryHits(findOverlaps(temp, met_gr_is_min))),
                      "met_is_min_length" = sum(width(met_gr_is_min)),
                      "met_is_plus" = length(queryHits(findOverlaps(temp, met_gr_is_plus))),
                      "met_is_plus_length" = sum(width(met_gr_is_plus)),
                      "ann_cpg_is_min" = length(queryHits(findOverlaps(temp, ann_cpg_is_min))),
                      "ann_cpg_is_min_length" = ann_cpg_is_min_length,
                      "ann_cpg_is_plus" = length(queryHits(findOverlaps(temp, ann_cpg_is_plus))),
                      "ann_cpg_is_plus_length" = ann_cpg_is_plus_length,
                      "idh_min" = length(queryHits(findOverlaps(temp, idh_gr_min))),
                      "idh_min_length" = sum(width(idh_gr_min)),
                      "idh_plus" = length(queryHits(findOverlaps(temp, idh_gr_plus))),
                      "idh_plus_length" = sum(width(idh_gr_plus)),
                      "idh_all_min" = length(queryHits(findOverlaps(temp, idh_gr[which(idh$KO_to_Par.t > 0),]))),
                      "idh_all_min_length" = sum(width(idh_gr[which(idh$KO_to_Par.t > 0),])),
                      "idh_all_plus" = length(queryHits(findOverlaps(temp, idh_gr[which(idh$KO_to_Par.t < 0),]))),
                      "idh_all_plus_length" = sum(width(idh_gr[which(idh$KO_to_Par.t < 0),])),
                      "idh_all_is_min" = length(queryHits(findOverlaps(temp, idh_gr[which((idh$KO_to_Par.t > 0) & (idh$KO_to_Par.cgi == "island"))]))),
                      "idh_all_is_min_length" = sum(width(idh_gr[which((idh$KO_to_Par.t > 0) & (idh$KO_to_Par.cgi == "island"))])),
                      "idh_all_is_plus" = length(queryHits(findOverlaps(temp, idh_gr[which((idh$KO_to_Par.t < 0) & (idh$KO_to_Par.cgi == "island"))]))),
                      "idh_all_is_plus_length" = sum(width(idh_gr[which((idh$KO_to_Par.t < 0) & (idh$KO_to_Par.cgi == "island"))])),
                      "all" = length(temp)
                      )
    t <- rbind(t, temp_tb)
  }
}
saveRDS(t, "/dcl02/leased/cglab/rscharpf/mnoe/PROJECTS/cfepigenetics/data/mouse_frag.rds")

all <- tibble()
for (i in meta_all) {
  if (file.exists(file.path(loc_n, "granges", paste0(i,".graft.rds")))) {
    print(i)
    temp <- readRDS(file.path(loc_n, "granges", paste0(i,".graft.rds")))
    table(width(temp))
    temp_tb <- tibble("case" = i,
                      "inplant" = meta[match(i, meta$pgdx.id),]$loc,
                      "amount" = length(temp),
                      "size" = as.numeric(names(table(width(temp)))),
                      "size_amount" = as.numeric(table(width(temp))),
                      "origin" = "graft")
    all <- rbind(all, temp_tb)
  }
}
for (i in meta_all) {
  if (file.exists(file.path(loc_n, "granges", paste0(i,".host.rds")))){
    print(i)
    temp <- readRDS(file.path(loc_n, "granges", paste0(i,".host.rds")))
    temp_tb <- tibble("case" = i,
                      "inplant" = meta[match(i, meta$pgdx.id),]$loc,
                      "amount" = length(temp),
                      "size" = as.numeric(names(table(width(temp)))),
                      "size_amount" = as.numeric(table(width(temp))),
                      "origin" = "host")
    all <- rbind(all, temp_tb)
  }
}
saveRDS(all, "/dcl02/leased/cglab/rscharpf/mnoe/PROJECTS/cfepigenetics/data/mouse_all_frag.rds")

       