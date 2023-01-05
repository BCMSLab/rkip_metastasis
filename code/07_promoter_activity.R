# Figure 7
# load libraries
library(tidyverse)
library(readxl)
library(cowplot)

# load data
regions <- read_excel('data/RKIP promoter activity (Luciferase).xlsx', 1)
drugs <- read_excel('data/RKIP promoter activity (Luciferase).xlsx', 2)

# test for significance
TukeyHSD(aov(lm(firefly/renilla ~ vector, regions)))

group_split(drugs, vector) %>%
  map(~{
    # test for significance
    TukeyHSD(aov(lm(firefly/renilla ~ treatment, .x)))
  })

# make figure
set.seed(1234)

mean_ctr <- (regions$firefly/regions$renilla)
mean_ctr <- mean(mean_ctr[regions$vector == 'CTR'])

pp1 <- regions %>%
  mutate(fold = (firefly/renilla)/mean_ctr,
         vector = factor(vector, unique(vector))) %>%
  group_by(vector) %>%
  mutate(ave = mean(fold), sd = sd(fold)) %>%
  ggplot(aes(x = vector, y = fold)) +
  geom_hline(yintercept = 1, lty = 2, color = 'darkgray') +
  geom_jitter(width = .3) +
  geom_errorbar(aes(ymin = ave - sd, ymax = ave + sd), width = .2, color = 'red') +
  geom_point(aes(y = ave), color = 'red') +
  labs(x = '', y = 'Relative Luciferase activity') +
  lims(y = c(0, 5)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid = element_blank(),
        axis.ticks.length=unit(.2, "cm"),
        axis.ticks = element_line(color= 'gray'),
        panel.spacing = unit('.05', 'cm'),
        axis.title = element_text(size = 10),
        plot.margin = unit(c(.2,.2,-.5,.2), 'cm'))

pp1

set.seed(123)

mean_veh <- drugs %>%
  mutate(fold = firefly/renilla) %>%
  filter(treatment == 'Vehicle') %>%
  group_by(vector) %>%
  summarise(mean_veh = mean(fold))

trts <- c("Vehicle", "Epirubicin", "Methotrexate", "Vorinostat", "Cisplatin", "Imatinib", "Sorafenib") 

pp2 <- drugs %>%
  mutate(treatment = factor(treatment, trts),
         fold = firefly/renilla) %>%
  left_join(mean_veh) %>%
  mutate(fold = fold/mean_veh) %>%
  group_by(vector, treatment) %>%
  mutate(ave = mean(fold), sd = sd(fold)) %>%
  ggplot(aes(x = treatment, y = fold)) +
  geom_hline(yintercept = 1, lty = 2, color = 'darkgray') +
  geom_jitter(width = .2) +
  geom_errorbar(aes(ymin = ave - sd, ymax = ave + sd), width = .2, color = 'red') +
  geom_point(aes(y = ave), color = 'red') +
  facet_grid(~vector) +
  lims(y = c(0, 3)) +
  labs(x = '', y = 'Relative Luciferase activity') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid = element_blank(),
        axis.ticks.length=unit(.2, "cm"),
        axis.ticks = element_line(color= 'gray'),
        panel.spacing = unit('.05', 'cm'),
        axis.title = element_text(size = 10),
        plot.margin = unit(c(.2,.2,-.5,.2), 'cm'))
pp2

# save
plot_grid(
  plot_grid(NULL, pp1, ncol = 1, rel_heights = c(1, 3)),
  pp2,
  rel_widths = c(1, 2),
  nrow = 1,
  scale = .9,
  labels = 'AUTO',
  label_fontface = 'plain'
) %>%
  ggsave(plot = .,
         filename = 'output/07_promoter_activity.png',
         width = 7, height = 3)


# # load libraries
# library(tidyverse)
# library(rstatix)
# 
# # load data
# drug_luc <- read_csv('drug_luc_raw.csv')
# glimpse(drug_luc)
# 
# # process data
# drug_luc_processed <- drug_luc %>%
#   mutate(
#     treatment_type = case_when(
#       treatment == 'vehicle' ~ 'control',
#       treatment %in% c('cisplatin', 'imatinib', 'sorafenib') ~ 'repressor',
#       treatment %in% c('epirubicin', 'vorinostat', 'methotrexate') ~ 'activator'
#     ),
#     treatment = fct_relevel(treatment, 'vehicle'),
#     treatment_type = fct_relevel(treatment_type, 'control')) %>%
#   group_by(vector) %>%
#   mutate(ratio = firefly / renilla,             # divide firefly by renilla
#          ratio = as.numeric(scale(ratio))) %>%  # scale to make values comparable
#   ungroup()
# glimpse(drug_luc_processed)  
# 
# # summarize data
# drug_luc_summary <- drug_luc_processed %>%
#   group_by(vector, treatment) %>%
#   get_summary_stats(ratio,
#                     type = 'full',                  # calculate summary stats
#                     show = c('mean', 'sd', 'max'))  # show only mean, sd and max
# glimpse(drug_luc_summary)
# 
# # test data
# drug_luc_test <- drug_luc_processed %>%
#   group_by(vector) %>%
#   pairwise_t_test(ratio ~ treatment) %>%     # compare all groups to each other (t.test)
#   filter(group1 == 'vehicle') %>%            # keep only comparisons treatment - vehicle
#   mutate(signif = ifelse(p < .05, '*', ''))  # signif if p < .05
# glimpse(drug_luc_test)
# 
# # visualize data
# drug_luc_processed %>%
#   left_join(drug_luc_summary) %>%
#   left_join(drug_luc_test, by = c('treatment' = 'group2', 'vector')) %>% 
#   mutate(treatment = fct_relevel(treatment, 'vehicle')) %>%
#   ggplot(aes(x = treatment, y = ratio, color = treatment_type)) +
#   geom_point(position = position_jitter(width = .2, seed = 1234)) +
#   geom_point(aes(y = mean), color = 'black') +
#   geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = .2, color = 'black') +
#   geom_text(aes(label = signif, y = max + .3), size = 10, color = 'black') +
#   facet_wrap(~vector) +
#   labs(x = '', y = 'Relative Fluorescence\n (Firefly / Renilla)') +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1),
#         legend.position = 'none',
#         axis.title.x = element_blank()) 

