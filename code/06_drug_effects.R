# Figure 7
# load libraries
library(tidyverse)
library(pcr)
library(cowplot)
# library(ggbreak)

# load data
all_drugs <- readxl::read_excel('data/ct_values.xlsx', 'alldrugs')
pathways <- readxl::read_excel('data/ct_values.xlsx', 'pathways')

drug_groups <- select(all_drugs, type = Group, group = Treatment) %>%
  unique() %>%
  mutate(type = factor(type, levels = c('CTR', 'Activator', 'Repressor')))

# all drugs
all_drugs <- all_drugs %>%
  # filter(`Type of assay` == 1) %>%
  filter(`Date of running` == '2021-11-24 UTC') %>%
  na.omit()
ct1 <- select(all_drugs, GAPDH, RKIP)

group <- pull(all_drugs, Treatment)
group <- relevel(factor(group), ref = 'DMSO')

res1 <- pcr_analyze(ct1,
                    group_var = group,
                    reference_gene = 'GAPDH',
                    reference_group = 'DMSO')

res1_test <- map_df(levels(group)[-1], ~{
  com <- c('DMSO', .x)
  ct <- filter(ct1, group %in% com)
  g <- factor(group[group %in% com], levels = com)
  
  pcr_test(ct,
           group_var = g,
           reference_gene = 'GAPDH',
           reference_group = 'DMSO') %>%
    mutate(group = .x)
}) %>%
  select(group, p_value)

# pcr_test(ct1,
#             group_var = group,
#             reference_gene = 'GAPDH',
#             reference_group = 'DMSO',
#          test = 'lm')

pp1 <- res1 %>%
    left_join(drug_groups) %>%
    left_join(res1_test) %>%
    mutate(sig = ifelse(p_value < .05, '*', '')) %>%
    ggplot(aes(x = group, y = relative_expression,
               ymin = lower, ymax = upper)) +
    geom_col() +
    geom_errorbar(width = .3) +
    # geom_text(aes(y = upper + .1, label = sig)) +
    facet_grid(~type, scales = 'free_x', space = 'free') +
    labs(x = '', y = 'Relative RKIP expression') +
    lims(y = c(0, 1.7)) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid = element_blank(),
        axis.ticks.length=unit(.2, "cm"),
        axis.ticks = element_line(color= 'gray'),
        panel.spacing = unit('.05', 'cm'),
        axis.title = element_text(size = 10),
        plot.margin = unit(c(.2,.2,-.5,.2), 'cm'))

pp1
# pathways
res2 <- pathways %>%
  filter(`Type of assay` == 2,
         Experiment != 7) %>%
  na.omit() %>%
  group_split(Experiment) %>%
  map_df(~{
    group <- .x$Treatment
    group <- relevel(factor(group), ref = 'DMSO')
    ct <- select(.x, GAPDH, RKIP, NME1, ESR1, RELA, SNAI1)
    
    pcr_analyze(ct,
                group_var = group,
                reference_group = 'DMSO',
                reference_gene = 'GAPDH') %>%
      as_tibble()
  })

# test significance
res2_test <- pathways %>%
  filter(`Type of assay` == 2,
         Experiment != 7) %>%
  na.omit() %>%
  group_split(Experiment) %>%
  map_df(~{
    group <- .x$Treatment
    group <- relevel(factor(group), ref = 'DMSO')
    ct <- select(.x, GAPDH, RKIP, NME1, ESR1, RELA, SNAI1)
    
    pcr_test(ct,
             group_var = group,
             reference_group = 'DMSO',
             reference_gene = 'GAPDH') %>%
      as_tibble() %>%
      mutate(group = unique(.x$Treatment)[-1])
  }) %>%
  select(group, gene, p_value) %>%
  unique()

pp2 <- res2 %>%
  filter(group != 'DMSO') %>%
  mutate(relative_expression = ifelse(relative_expression > 3, 3, relative_expression)) %>%
  ggplot(aes(x = group, y = gene, fill = relative_expression)) +
  geom_tile() +
  labs(x = '', y = '') +
  scale_fill_gradient2(low = 'blue', high = 'red', midpoint = 1) +
  theme(legend.position = 'none') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid = element_blank(),
        axis.ticks.length=unit(.2, "cm"),
        axis.ticks = element_line(color= 'gray'),
        panel.spacing = unit('.05', 'cm'),
        axis.title = element_text(size = 10),
        plot.margin = unit(c(.2,.2,-.5,-.2), 'cm'))
pp2

# save
plot_grid(NULL, pp1, pp2, nrow = 1, scale = .9,
          rel_widths = c(.8, 1, 1),
          labels = 'AUTO',
          label_fontface = 'plain') %>%
  ggsave(plot = .,
         filename = 'output/06_drug_effects.png',
         width = 7, height = 3)
