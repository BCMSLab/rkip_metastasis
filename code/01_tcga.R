# Figure 1
# load libraries
library(tidyverse)
library(reshape2)
library(cowplot)

# # load data
# # METABRIC
# metabric <- list(cna = 'data/brca_metabric/data_cna.txt',
#      mrna = 'data/brca_metabric/data_mrna_agilent_microarray_zscores_ref_diploid_samples.txt') %>%
#   map(~read_tsv(.x) %>% filter(Hugo_Symbol == 'PEBP1') %>% dplyr::select(starts_with('MB'))) %>%
#   bind_rows() %>%
#   t() %>%
#   as.data.frame() %>%
#   set_names(c('cna', 'mrna')) %>%
#   rownames_to_column('SAMPLE_ID') %>%
#   right_join(read_tsv('data/brca_metabric/data_clinical_sample.txt', skip = 4)) %>%
#   mutate(grade = GRADE)
# 
# # MBC
# mbc <- list(cna = 'data/brca_mbcproject_wagle_2017/data_cna.txt',
#      mrna = 'data/brca_mbcproject_wagle_2017/data_mrna_seq_v2_rsem_zscores_ref_diploid_samples.txt') %>%
#   map(~read_tsv(.x) %>% filter(Hugo_Symbol == 'PEBP1') %>% dplyr::select(starts_with('MB'))) %>%
#   bind_rows() %>%
#   t() %>%
#   as.data.frame() %>%
#   set_names(c('cna', 'mrna')) %>%
#   rownames_to_column('SAMPLE_ID') %>%
#   right_join(read_tsv('data/brca_mbcproject_wagle_2017/data_clinical_sample.txt', skip = 4)) %>%
#   dplyr::select(ER_STATUS = BX_ER, PR_STATUS = BX_PR, HER2_STATUS = BX_HER2OVERALL, everything()) %>%
#   mutate(grade = as.numeric(as.factor(BX_GRADE))) 
# 
# # combine
# dat <- bind_rows(metabric = metabric, mbc = mbc, .id = 'study') %>% 
#   as_tibble() %>%
#   mutate_at(vars(ends_with('STATUS')), tolower) %>%
#   mutate(study = toupper(study))
# 
# write_rds(dat, 'data/tcga.rds')

dat <- read_rds('data/tcga.rds')

# plots
# percentages of alterations
pp1 <- as_tibble(table(dat$study, dat$cna), .name_repair = 'unique') %>%
  setNames(c('study', 'alteration', 'number')) %>%
  group_split(study) %>%
  map(~{
    ggplot(data = .x, aes(x = '', y = number, fill = alteration)) +
      geom_bar(stat = 'identity', width = 1) +
      coord_polar('y', start = 1) +
      facet_wrap(~study) +
      theme_void() +
      theme(legend.position = 'none',
            plot.margin = unit(c(0, 3, 0, 0), 'cm'))
  })
pp1
round(table(dat$study, dat$cna) / rowSums(table(dat$study, dat$cna)), 2)

1# expression vs cna
pp2 <- dat %>%
  dplyr::select(study, cna, mrna) %>%
  na.omit() %>%
  ggplot(aes(x = as.factor(cna), y = mrna)) +
  geom_boxplot(position = 'dodge', size = 1) +
  labs(x = 'Copy Number Alterations (CNA)',
       y = 'RKIP Expression (standardized)') +
  geom_hline(yintercept = 0, lty = 2, color = 'red') + 
  scale_y_continuous(limits = c(-3.5, 3.5),
                     breaks = seq(-3, 3)) +
  facet_grid(~study) + 
  theme(legend.position = 'top',
        panel.grid = element_blank(),
        axis.ticks.length=unit(.2, "cm"),
        axis.ticks = element_line(color= 'gray'),
        axis.title = element_text(size = 14))
pp1

# expression vs receptor status
set.seed(123)
pp3 <- dat %>%
  dplyr::select(study, mrna, ER_STATUS, HER2_STATUS, PR_STATUS) %>%
  gather(receptor, status, ER_STATUS, HER2_STATUS, PR_STATUS) %>%
  mutate(receptor = str_replace_all(receptor, '_STATUS', '')) %>%
  filter(status %in% c('positive', 'negative')) %>%
  mutate(status = ifelse(status == 'positive', '+ve', '-ve')) %>%
  group_by(study, receptor, status) %>%
  mutate(ave = mean(mrna, na.rm = TRUE),
         sd = sd(mrna, na.rm = TRUE)) %>%
  ggplot(aes(x = status, y = mrna)) +
  facet_grid(receptor ~ study) +
  geom_jitter(width = .3, alpha = .5) +
  geom_point(aes(x = status, y = ave), color = 'red') +
  geom_errorbar(aes(x = status, ymin = ave-sd, ymax = ave+sd), color = 'red', width = .2) +
  geom_hline(yintercept = 0, lty = 2, color = 'red') + 
  # scale_y_continuous(limits = c(-3.5, 3.5),
  #                    breaks = seq(-3, 3)) +
  labs(x = 'Receptor Status',
       y = 'RKIP Expression (standardized)',
       fill = '') +
  theme(legend.position = 'top',
        panel.grid = element_blank(),
        axis.ticks.length=unit(.2, "cm"),
        axis.ticks = element_line(color= 'gray'),
        panel.spacing = unit('.05', 'cm'),
        axis.title = element_text(size = 14))
pp2

plot_grid(plot_grid(pp1[[1]], pp1[[2]], ncol = 1),
          pp2, pp3,
          scale = .9,
          nrow = 1,
          labels = 'AUTO',
          label_fontface = 'plain') %>%
  ggsave(plot = .,
         filename = 'output/01_tcga.png',
         height = 4, width = 12)
