# Figure 4
# load libraries
library(tidyverse)
# library(limma)
library(fgsea)
# library(org.Hs.eg.db)
library(cowplot)

# # load data
# ## overexpression in MDA
# df <- read_tsv('data/GSE128983_Normalized_counts.txt')
# mat <- as.matrix(df[, -1])
# rownames(mat) <- as.character(df[,1, drop = TRUE])
# 
# group <- str_split(colnames(mat), '_', simplify = TRUE)[, 1]
# group <- ifelse(group == 'CONTROL', 'WT', 'OE')
# 
# oe <- list(mat = mat, group = group)
# 
# ## knockdown in mcf7
# kd <- read_rds('data/msg_kd_eset.rds')
# kd <- kd[, kd$target %in% c('Scramble', 'PEBP1')]
# group <- ifelse(kd$target == 'Scramble', 'WT', 'KD')
# kd <- list(mat = exprs(kd), group = group)
# 
# plot(kd$mat['PEBP1',] ~ as.factor(group))
# 
# dat <- list(MDA = oe, MCF7 = kd) %>%
#   map_df(~{
#     tibble(expr = as.numeric(scale(.x$mat['PEBP1',])),
#            group = .x$group)
#   }, .id = 'cell')
# 
# write_csv(dat, 'data/kd_oe.csv')
dat <- read_csv('data/kd_oe.csv')
dat <- mutate(dat, cell = ifelse(cell == 'MDA', 'MDA_MB_231', cell))
set.seed(123)
pp1 <- dat %>% 
  group_by(cell, group) %>%
  mutate(ave = mean(expr), sd = sd(expr)) %>%
  ggplot(aes(x = group, y = expr)) +
  geom_jitter(width = .2) +
  geom_point(aes(x = group, y = ave), color = 'red') +
  geom_hline(yintercept = 0, lty = 2, color = 'red') +
  geom_errorbar(aes(x = group, ymin = ave -sd, ymax = ave + sd), color = 'red', width = .1) +
  facet_grid(~cell, scales = 'free_x') +
  labs(x = 'Experimental group',
       y = 'RKIP Expression (standardized)') + 
  theme(panel.grid = element_blank(),
        axis.ticks.length=unit(.2, "cm"),
        axis.ticks = element_line(color= 'gray'),
        axis.title = element_text(size = 14))

pp1

# # differential expression
# fac <- relevel(as.factor(kd$group), ref = 'WT')
# mod <- model.matrix(~fac)
# mat <- kd$mat[rowSums(kd$mat) > 10,]
# fit <- lmFit(log2(mat), mod)
# fit <- eBayes(fit)
# tt1 <- topTable(fit, genelist = rownames(kd$mat), number = Inf)
# 
# table(tt1$adj.P.Val < .2)
# 
# fac <- relevel(as.factor(oe$group), ref = 'WT')
# mod <- model.matrix(~fac)
# mat <- oe$mat[rowSums(oe$mat) > 10,]
# fit <- lmFit(log2(mat + 1), mod)
# fit <- eBayes(fit)
# tt2 <- topTable(fit, genelist = rownames(oe$mat), number = Inf)
# 
# tt <- bind_rows(list(MDA = tt2, MCF7 = tt1), .id = 'cell') %>% as_tibble()
# table(tt2$adj.P.Val < .2)
# 
# # enrichment
# # GO:0001837
# pathway <- select(org.Hs.eg.db,
#                   'GO:0010719',
#                   'SYMBOL',
#                   'GO')
# 
# emt <- unique(pathway$SYMBOL)
# 
# fc <- split(tt$t, tt$cell)
# nms <- split(tt$ID, tt$cell)
# stats <- map2(fc, nms, ~{
#   vec <- .x
#   names(vec) <- .y
#   vec
# })
# 
# enrichment <- map(stats, ~{
#   fgsea(list(EMT = emt),
#         stats = .x)
# })
# 
# write_rds(stats, 'data/stats.rds')
# write_rds(emt, 'data/emt.rds')

stats <- read_rds('data/stats.rds')
emt <- read_rds('data/emt.rds')

pp2 <- plotEnrichment(emt, stats$MCF7) +
  facet_grid(~'MCF7') +
  theme_gray()+
  theme(panel.grid = element_blank(),
        axis.ticks.length=unit(.2, "cm"),
        axis.ticks = element_line(color= 'gray'),
        axis.title = element_text(size = 14)) +
  labs(x = 'Gene rank (KD vs WT)',
       y = 'Enrichment Scores (ES)')

pp2 <- pp2$data %>%
  ggplot() +
  geom_line(aes(x = x, y = y), color = 'blue', size = 2) +
  geom_segment(aes(x = x, xend = x, y = -.3, yend = -.22)) +
  theme(panel.grid = element_blank(),
        axis.ticks.length=unit(.2, "cm"),
        axis.ticks = element_line(color= 'gray'),
        axis.title = element_text(size = 14)) +
  geom_hline(yintercept = -.26) +
  geom_hline(yintercept = 0, color = 'red', lty = 2) +
  labs(x = 'Gene rank (KD vs WT)',
       y = 'Enrichment Scores (ES)') +
  facet_grid(~'MCF7')
pp2

pp3 <- plotEnrichment(emt, stats$MDA) +
  facet_grid(~'MDA_MB_231') +
  theme_gray()+
  theme(panel.grid = element_blank(),
        axis.ticks.length=unit(.2, "cm"),
        axis.ticks = element_line(color= 'gray'),
        axis.title = element_text(size = 14)) +
  labs(x = 'Gene rank (OE vs WT)',
       y = 'Enrichment Scores (ES)')

pp3 <- pp3$data %>%
  ggplot() +
  geom_line(aes(x = x, y = y), color = 'blue', size = 2) +
  geom_segment(aes(x = x, xend = x, y = -.05, yend = -.11)) +
  theme(panel.grid = element_blank(),
        axis.ticks.length=unit(.2, "cm"),
        axis.ticks = element_line(color= 'gray'),
        axis.title = element_text(size = 14)) +
  geom_hline(yintercept = -.08) +
  scale_y_continuous(breaks = seq(0, .3, .1)) +
  geom_hline(yintercept = 0, color = 'red', lty = 2) +
  labs(x = 'Gene rank (OE vs WT)',
       y = 'Enrichment Scores (ES)') +
  facet_grid(~'MDA_MB_231')
pp3

plot_grid(pp1, pp2, pp3,
          scale = .9,
          nrow = 1,
          labels = 'AUTO',
          label_fontface = 'plain') %>%
  ggsave(plot = .,
         filename = 'output/04_geo.png',
         height = 3.5, width = 10)
