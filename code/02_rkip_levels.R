# Figure 2
# load libraries
library(tidyverse)
library(readxl)
library(cowplot)
library(pcr)

# load data
rkip_protein <- read_excel('data/RKIP_level_in_cells.xlsx', 1)
rkip_mrna <- read_excel('data/RKIP_level_in_cells.xlsx', 2)

# RKIP levels in different cancer cell lines
# visualize the band intensities
set.seed(123)
pp1 <- rkip_protein %>%
  pivot_wider(values_from = value, names_from = protein, values_fn = list) %>%
  unnest(cols = c(RKIP, b_actin)) %>%
  mutate(relative = RKIP/b_actin,
         cell = factor(cell, levels = unique(rkip_protein$cell))) %>%
  group_by(cell) %>%
  mutate(ave = mean(relative), sd = sd(relative)) %>%
  ggplot(aes(x = cell, y = relative)) +
  geom_jitter(width = .2) +
  geom_errorbar(aes(ymin = ave - sd, ymax = ave + sd), width = .2, color = 'red') +
  geom_point(aes(y = ave), color = 'red') +
  lims(y = c(0, 2.5)) +
  labs(x = '', y = 'Relative RKIP level') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid = element_blank(),
        axis.ticks.length=unit(.2, "cm"),
        axis.ticks = element_line(color= 'gray'),
        panel.spacing = unit('.05', 'cm'),
        axis.title = element_text(size = 10),
        plot.margin = unit(c(.2,.2,-.5,.2), 'cm'))
pp1

# test significance across replicates
rkip_protein_cells <- rkip_protein %>%
  pivot_wider(values_from = value, names_from = protein, values_fn = list) %>%
  unnest(cols = c(RKIP, b_actin)) %>%
  mutate(relative = RKIP/b_actin,
         cell = factor(cell, levels = unique(rkip_protein$cell))) %>%
  with(split(relative, cell))

map_df(rkip_protein_cells, ~{
  broom::tidy(t.test(rkip_protein_cells$MCF7, .x))
}, .id = 'cell') 

# visualize the relative expression of RKIP in the two cell lines
rkip_mrna2 <- rkip_mrna %>%
  pivot_wider(values_from = value, names_from = gene, values_fn = list) %>%
  unnest(cols = c(RKIP, GAPDH))
ct <- select(rkip_mrna2, -cell)
group <- rkip_mrna2$cell

res <- pcr_analyze(ct, 
            reference_gene = 'GAPDH',
            reference_group = 'MCF7',
            group_var = group)

pp2 <- res %>%
  ggplot(aes(x = group, y = relative_expression)) +
  geom_col() +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = .3) +
  labs(x = '', y = 'Relative RKIP expression') +
  lims(y = c(0, 2)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid = element_blank(),
        axis.ticks.length=unit(.2, "cm"),
        axis.ticks = element_line(color= 'gray'),
        panel.spacing = unit('.05', 'cm'),
        axis.title = element_text(size = 10),
        plot.margin = unit(c(.2,.2,-.5,.2), 'cm'))

# test the significance across replicates
pcr_test(ct, 
            reference_gene = 'GAPDH',
            reference_group = 'MCF7',
            group_var = group)

# save the figures
plot_grid(
  plot_grid(NULL, pp1, ncol = 1, rel_heights = c(1, 3)),
  pp2,
  rel_widths = c(1, .8),
  nrow = 1,
  scale = .9,
  labels = 'AUTO',
  label_fontface = 'plain'
) %>%
  ggsave(plot = .,
         filename = 'output/02_rkip_levels.png',
         width = 7, height = 3)