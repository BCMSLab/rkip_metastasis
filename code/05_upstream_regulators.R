# Figure 5
# load libraries
library(tidyverse)
library(readxl)
library(cowplot)
library(pcr)
# library(ggbreak)

# load data
snai1 <- read_excel('data/WB_RKIP_SNAI_and_NME1.xlsx', 1)
nme1 <- read_excel('data/WB_RKIP_SNAI_and_NME1.xlsx', 2)

# make figure
d <- snai1
d$RKIP <- d$rkip/d$beta_actin
d$RKIP <- d$RKIP/mean(d$RKIP[d$cell == 'WT'])

# test significance across replicates
TukeyHSD(aov(lm(RKIP ~ cell, d)))

set.seed(123)
pp1 <- d %>%
  mutate(cell = factor(cell, levels = c('WT', 'OE', 'KD'))) %>%
  group_by(cell) %>%
  mutate(ave = mean(RKIP), sd = sd(RKIP)) %>%
  ggplot(aes(x = cell, y = RKIP)) +
  geom_jitter(width = .2) +
  geom_errorbar(aes(ymin = ave - sd, ymax = ave + sd), width = .2, color = 'red') +
  geom_point(aes(y = ave), color = 'red') +
  lims(y = c(0, 5)) +
  labs(x = '', y = 'Relative RKIP level') +
  theme(panel.grid = element_blank(),
        axis.ticks.length=unit(.2, "cm"),
        axis.ticks = element_line(color= 'gray'),
        panel.spacing = unit('.05', 'cm'),
        axis.title = element_text(size = 10),
        plot.margin = unit(c(.2,.2,-.5,.2), 'cm'))
pp1

d <- nme1
d$RKIP <- d$rkip/d$beta_actin
d$RKIP <- d$RKIP/mean(d$RKIP[d$cell == 'WT'])

# test significance across replicates
TukeyHSD(aov(lm(RKIP ~ cell, d)))

set.seed(123)
pp3 <- d %>%
  mutate(cell = relevel(as.factor(cell), ref = 'WT')) %>%
  group_by(cell) %>%
  mutate(ave = mean(RKIP), sd = sd(RKIP)) %>%
  ggplot(aes(x = cell, y = RKIP)) +
  geom_jitter(width = .2) +
  geom_errorbar(aes(ymin = ave - sd, ymax = ave + sd), width = .2, color = 'red') +
  geom_point(aes(y = ave), color = 'red') +
  labs(x = '', y = 'Relative RKIP level') +
  scale_y_continuous(breaks = seq(0, 10, 2), limits = c(0, 11)) +
  theme(panel.grid = element_blank(),
        axis.ticks.length=unit(.2, "cm"),
        axis.ticks = element_line(color= 'gray'),
        panel.spacing = unit('.05', 'cm'),
        axis.title = element_text(size = 10),
        plot.margin = unit(c(.2,.2,-.5,.2), 'cm'))
pp3

# mrna
snai1 <- read_excel('data/WB_RKIP_SNAI_and_NME1.xlsx', 3)
ct <- select(snai1, -cell)
group <- snai1$cell

# test significance
pcr_test(ct,
         group_var = group,
         reference_gene = 'GAPDH',
         reference_group = 'WT',
         test = 'lm')

res <- pcr_analyze(ct,
                   group_var = group,
                   reference_gene = 'GAPDH',
                   reference_group = 'WT')
pp2_a <- res %>%
  filter(gene == 'SNAI1') %>%
  mutate(group = relevel(as.factor(group), ref = 'WT')) %>%
  ggplot(aes(x = group, y = relative_expression)) +
  geom_col() +
  facet_wrap(~gene, scales = 'free_y') +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = .2) +
  labs(x = '', y = 'Relative gene expression') +
  lims(y = c(0, 11.8)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid = element_blank(),
        axis.ticks.length=unit(.2, "cm"),
        axis.ticks = element_line(color= 'gray'),
        panel.spacing = unit('.05', 'cm'),
        axis.title = element_text(size = 10),
        plot.margin = unit(c(.2,.2,-.5,.2), 'cm'))

pp2_b <- res %>%
  filter(gene == 'RKIP') %>%
  mutate(group = relevel(as.factor(group), ref = 'WT')) %>%
  ggplot(aes(x = group, y = relative_expression)) +
  geom_col() +
  facet_wrap(~gene, scales = 'free_y') +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = .2) +
  labs(x = '', y = '') +
  lims(y = c(0, 1.8)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid = element_blank(),
        axis.ticks.length=unit(.2, "cm"),
        axis.ticks = element_line(color= 'gray'),
        panel.spacing = unit('.05', 'cm'),
        axis.title = element_text(size = 10),
        plot.margin = unit(c(.2,.2,-.5,-.2), 'cm'))

nme1 <- read_excel('data/WB_RKIP_SNAI_and_NME1.xlsx', 4)

ct <- select(nme1, -cell)
group <- nme1$cell

res <- pcr_analyze(ct,
                   group_var = group,
                   reference_gene = 'GAPDH',
                   reference_group = 'WT')

# test significance
pcr_test(ct,
         group_var = group,
         reference_gene = 'GAPDH',
         reference_group = 'WT',
         test = 'lm')

pp4 <- res %>%
  mutate(group = relevel(as.factor(group), ref = 'WT')) %>%
  ggplot(aes(x = group, y = relative_expression)) +
  geom_col() +
  lims(y = c(0, 2.7)) +
  facet_wrap(~gene) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = .2) +
  labs(x = '', y = 'Relative gene expression') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid = element_blank(),
        axis.ticks.length=unit(.2, "cm"),
        axis.ticks = element_line(color= 'gray'),
        panel.spacing = unit('.05', 'cm'),
        axis.title = element_text(size = 10),
        plot.margin = unit(c(.2,.2,-.5,.2), 'cm'))
pp4

# save
plot_grid(
  plot_grid(pp2_a, pp2_b, rel_widths = c(1, .9)),
  plot_grid(NULL, pp1, ncol = 1, rel_heights = c(1, 1.5)),
  pp4,
  plot_grid(NULL, pp3, ncol = 1, rel_heights = c(1, 1.5)),
  nrow = 2,
  scale = .9,
  labels = 'AUTO',
  label_fontface = 'plain'
) %>%
  ggsave(plot = .,
         filename = 'output/05_upstream_regulators.png',
         width = 7, height = 7)
