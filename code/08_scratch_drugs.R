# Figure 3
# load libraries
library(tidyverse)
library(readxl)

# load data
scratch <- read_excel('data/Scratch assay on MCF-7 and MDA-MB-231 cells.xlsx', 2)

# test significance across replicates
scratch %>%
  group_split(treatment) %>%
  map(~{TukeyHSD(aov(lm(.x$value ~ .x$time)))})

# visualize the migration area over time
d <- scratch
ctr <- mean(d$value[d$treatment == 'DMSO' & d$time == '0h'] )

d$relative <- (d$value/ctr) * 100

set.seed(123)
pp1 <- d %>%
  mutate(treatment = paste0('MCF7 + ', treatment)) %>%
  group_by(time, treatment) %>%
  mutate(ave = mean(relative), sd = sd(relative)) %>%
  ggplot(aes(x = time, y = relative)) +
  geom_jitter(width = .2) +
  geom_errorbar(aes(ymin = ave - sd, ymax = ave + sd), width = .2, color = 'red') +
  geom_point(aes(y = ave), color = 'red') +
  facet_wrap(~treatment, ncol = 1) +
  lims(y = c(0, 125)) +
  labs(x = '', y = 'Migration area (%)') +
  theme(panel.grid = element_blank(),
        axis.ticks.length=unit(.2, "cm"),
        axis.ticks = element_line(color= 'gray'),
        panel.spacing = unit('.05', 'cm'),
        axis.title = element_text(size = 10),
        plot.margin = unit(c(.2,.2,-.2,.2), 'cm'))
pp1


scratch2 <- read_excel('data/Scratch assay on MCF-7 and MDA-MB-231 cells.xlsx', 3)

# test significance across replicates
scratch %>%
  group_split(treatment) %>%
  map(~{TukeyHSD(aov(lm(.x$value ~ .x$time)))})

# visualize the migration area over time
d <- scratch2
ctr <- mean(d$value[d$treatment == 'Epirubicin' & d$time == '0h'] )

d %>%
  ggplot(aes(x = time, y = value/ctr)) +
  geom_point()

d$relative <- (d$value/ctr) * 100

set.seed(123)
pp2 <- d %>%
  mutate(treatment = paste0('MDA_MB_231 + ', treatment)) %>%
  group_by(time, treatment) %>%
  mutate(ave = mean(relative), sd = sd(relative)) %>%
  ggplot(aes(x = time, y = relative)) +
  geom_jitter(width = .2) +
  geom_errorbar(aes(ymin = ave - sd, ymax = ave + sd), width = .2, color = 'red') +
  geom_point(aes(y = ave), color = 'red') +
  facet_wrap(~treatment, ncol = 1) +
  lims(y = c(0, 130)) +
  labs(x = '', y = 'Migration area (%)') +
  theme(panel.grid = element_blank(),
        axis.ticks.length=unit(.2, "cm"),
        axis.ticks = element_line(color= 'gray'),
        panel.spacing = unit('.05', 'cm'),
        axis.title = element_text(size = 10),
        plot.margin = unit(c(.2,.2,-.2,.2), 'cm'))
pp2

# save 
plot_grid(NULL, pp1, NULL, pp2, rel_widths = c(2, 1),
          nrow = 2,
          scale = .9,
          labels = 'AUTO',
          label_fontface = 'plain') %>%
  
  ggsave(plot = .,
         filename = 'output/08_scratch_drugs.png',
         width = 8, height = 8)

