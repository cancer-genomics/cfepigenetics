t <- readRDS("../../mouse_frag.rds")
#t <- t[which(t$top == 1000),]
t <- t[-which(t$id == "PGDX27350P"),] ### remove a file with problems
t <- t[-which(t$id == "PGDX27354P"),] ### remove a file with problems
t <- t[-which(t$inplant == "brain"),]
t[which(t$mut == "MT"),]$mut <- "mutant"
t[which(t$mut == "WT"),]$mut <- "wild-type"
t2 <- t %>%
  group_by(mut,id,top) %>%
  summarise(rna_min = (sum(rna_plus) / sum(all)) / rna_plus_length,
            rna_plus = (sum(rna_min) / sum(all)) / rna_min_length,
            met_min = (sum(met_min) / sum(all)) / met_min_length,
            met_plus = (sum(met_plus) / sum(all)) / met_plus_length,
            met_is_min = (sum(met_is_min) / sum(all)) / met_is_min_length,
            met_is_plus = (sum(met_is_plus) / sum(all)) / met_is_plus_length,
            ann_is_min = (sum(ann_cpg_is_min) / sum(all)) / ann_cpg_is_min_length,
            ann_is_plus = (sum(ann_cpg_is_plus) / sum(all)) / ann_cpg_is_plus_length,
            idh_all_min = (sum(idh_all_min) / sum(all)) / idh_all_min_length,
            idh_all_plus = (sum(idh_all_plus) / sum(all)) / idh_all_plus_length,
            idh_all_is_min = (sum(idh_all_is_plus) / sum(all)) / idh_all_is_plus_length,
            idh_all_is_plus = (sum(idh_all_is_min) / sum(all)) / idh_all_is_min_length,
            all = sum(all))
t2l <- reshape2::melt(t2,id.vars = c("mut", "id", "top"))

t2l_sum_met <- t2l[which(t2l$variable %in% c('met_is_min', 'met_is_plus')),]
t2l_sum_met$group <- "Methylation" 
t2l_sum_rna <- t2l[which((t2l$variable %in% c('rna_min', 'rna_plus'))),]
t2l_sum_rna$group <- "Expression"
t2l_sum <- rbind (t2l_sum_met, t2l_sum_rna)
t2l_sum$group <- factor(t2l_sum$group, levels=c("Methylation", "Expression"))
t2l_sum_max <- t2l_sum %>%
  group_by(top, variable) %>%
  summarise(max = max(value))
t2l_sum$max <- t2l_sum_max[match(paste0(t2l_sum$top, t2l_sum$variable), paste0(t2l_sum_max$top, t2l_sum_max$variable)),]$max
t2l_sum$norm_value <- t2l_sum$value / t2l_sum$max

set.seed(1234)
labels <- rep(c(0, 1), each=3)
observed <- paste(labels, collapse="")
tmp <- replicate(10000, sample(labels, size=6, replace=FALSE))
chars <- apply(tmp, 2, paste, collapse="")
perms <- unique(chars)
keep <- tmp[, chars %in% perms & !duplicated(chars)]

#permutation using medians
all <- tibble()
set.seed(1234)
for (j in c(1:10000)) {
  print(j)
  P <- j
  for (i in c(500, 1000, 2000, 3000, 4000, 5000)) {
    t2l_sum_temp <- t2l_sum[which(t2l_sum$top == i),]
    for (q in unique(t2l_sum_temp$variable)) {
      t2l_sum_temp_temp <-  t2l_sum_temp[which( t2l_sum_temp$variable == q),]
      pindex <- keep[, sample(c(1:20), size = 1, replace = T)] + 1
      permutation <- c("mutant", "wild-type")[pindex]
      t2l_sum_temp_temp$mut <- permutation
      #t2l_sum_temp_temp[which(t2l_sum_temp_temp$mut == "mutant"),]$value <- sample(t2l_sum_temp_temp[which(t2l_sum_temp_temp$mut == "mutant"),]$value, size = 3, replace = T)
      #t2l_sum_temp_temp[which(t2l_sum_temp_temp$mut == "wild-type"),]$value <- sample(t2l_sum_temp_temp[which(t2l_sum_temp_temp$mut == "wild-type"),]$value, size = 3, replace = T)
      sum <- t2l_sum_temp_temp %>%
        group_by(mut) %>%
        summarise(med = median(norm_value))
      if (q %in% c("met_is_min", "rna_min")) {
        sum_temp <- tibble("P" = P,
                           "top" = i,
                           "group" = q,
                           "direction" = sign(sum[which(sum$mut == "wild-type"),]$med - sum[which(sum$mut == "mutant"),]$med))
      }
      else {
        sum_temp <- tibble("P" = P,
                           "top" = i,
                           "group" = q,
                           "direction" = sign(sum[which(sum$mut == "mutant"),]$med - sum[which(sum$mut == "wild-type"),]$med))
      }
      all <- rbind(all, sum_temp)
    }
  }
}
all[which(all$direction == -1),]$direction <- 0
all_sum <- all %>%
  group_by(P) %>%
  summarise(count_positive = sum(direction))
hist(all_sum$count_positive, breaks = 30)
all_sum_median <- all_sum
hist(all_sum_median$count_positive, breaks = 30)
max(all_sum_median$count_positive)
saveRDS(all_sum_median,"../data/mouse_permutations_median.rds")

#permutation using means
all <- tibble()
for (j in c(1:10000)) {
  print(j)
  P <- j
  for (i in c(500, 1000, 2000, 3000, 4000, 5000)) {
    t2l_sum_temp <- t2l_sum[which(t2l_sum$top == i),]
    for (q in unique(t2l_sum_temp$variable)) {
      t2l_sum_temp_temp <-  t2l_sum_temp[which( t2l_sum_temp$variable == q),]
      pindex <- keep[, sample(c(1:20), size = 1, replace = T)] + 1
      permutation <- c("mutant", "wild-type")[pindex]
      t2l_sum_temp_temp$mut <- permutation
      #t2l_sum_temp_temp[which(t2l_sum_temp_temp$mut == "mutant"),]$value <- sample(t2l_sum_temp_temp[which(t2l_sum_temp_temp$mut == "mutant"),]$value, size = 3, replace = T)
      #t2l_sum_temp_temp[which(t2l_sum_temp_temp$mut == "wild-type"),]$value <- sample(t2l_sum_temp_temp[which(t2l_sum_temp_temp$mut == "wild-type"),]$value, size = 3, replace = T)
      sum <- t2l_sum_temp_temp %>%
        group_by(mut) %>%
        summarise(med = mean(norm_value))
      if (q %in% c("met_is_min", "rna_min")) {
        sum_temp <- tibble("P" = P,
                           "top" = i,
                           "group" = q,
                           "direction" = sign(sum[which(sum$mut == "wild-type"),]$med - sum[which(sum$mut == "mutant"),]$med))
      }
      else {
        sum_temp <- tibble("P" = P,
                           "top" = i,
                           "group" = q,
                           "direction" = sign(sum[which(sum$mut == "mutant"),]$med - sum[which(sum$mut == "wild-type"),]$med))
      }
      all <- rbind(all, sum_temp)
    }
  }
}
all[which(all$direction == -1),]$direction <- 0
all_sum <- all %>%
  group_by(P) %>%
  summarise(count_positive = sum(direction))
all_sum_mean <- all_sum
hist(all_sum_mean$count_positive, breaks = 30)
max(all_sum_mean$count_positive)
saveRDS(all_sum_mean,"../data/mouse_permutations_mean.rds")

#calculating the amount of times a metric is in the correct direction, when using correct labels
all_norm <- tibble()
for (i in c(500, 1000, 2000, 3000, 4000, 5000)){
  t2l_sum_temp <- t2l_sum[which(t2l_sum$top == i),]
  for (q in unique(t2l_sum$variable)) {
    t2l_sum_temp_temp <-  t2l_sum_temp[which( t2l_sum_temp$variable == q),]
    sum <- t2l_sum_temp_temp %>%
      group_by(mut) %>%
      summarise(med = mean(norm_value))
    if (q %in% c("met_is_min", "rna_min")) {
      sum_temp <- tibble("P" = P,
                         "top" = i,
                         "group" = q,
                         "direction" = sign(sum[which(sum$mut == "wild-type"),]$med - sum[which(sum$mut == "mutant"),]$med))
    }
    else {
      sum_temp <- tibble("P" = P,
                         "top" = i,
                         "group" = q,
                         "direction" = sign(sum[which(sum$mut == "mutant"),]$med - sum[which(sum$mut == "wild-type"),]$med))
    }
    all_norm <- rbind(all_norm, sum_temp)
  }
}
all_norm[which(all_norm$direction == -1),]$direction <- 0
all_norm_sum <- all_norm %>%
  group_by(P) %>%
  summarise(count_positive = sum(direction))
