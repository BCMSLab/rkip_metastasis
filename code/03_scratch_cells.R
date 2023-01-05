# Figure 3
# load libraries
library(tidyverse)
library(readxl)

# load data
scratch <- read_excel('data/Scratch assay on MCF-7 and MDA-MB-231 cells.xlsx', 1)

# test significance of the time points as groups
scratch %>%
  group_split(cell) %>%
  map(~{TukeyHSD(aov(lm(.x$value ~ .x$time)))})

# visualize the migration area over time
pp <- map(unique(scratch$cell), ~{
  ave_0h <- mean(scratch$value[scratch$time == '0h'& scratch$ cell == .x])
  set.seed(123)
  scratch %>%
    filter(cell == .x) %>%
    mutate(normalized = value / ave_0h,
           normalized = normalized * 100) %>%
    group_by(time) %>%
    mutate(ave = mean(normalized), sd = sd(normalized)) %>%
    ggplot(aes(x = time, y = normalized)) +
    geom_jitter(width = .2) +
    geom_errorbar(aes(ymin = ave - sd, ymax = ave + sd), width = .2, color = 'red') +
    geom_point(aes(y = ave), color = 'red') +
    lims(y = c(0, 130)) +
    labs(x = '', y = 'Migration area (%)') +
    theme(panel.grid = element_blank(),
          axis.ticks.length=unit(.2, "cm"),
          axis.ticks = element_line(color= 'gray'),
          panel.spacing = unit('.05', 'cm'),
          axis.title = element_text(size = 10),
          plot.margin = unit(c(1,.2,-.2,.2), 'cm'))
})

# save
plot_grid(
  NULL, pp[[1]],
  NULL, pp[[2]],
  ncol = 2,
  rel_widths = c(1, .5),
  scale = .95,
  labels = c('A', '', 'B', ''),
  label_fontface = 'plain'
) %>%
  ggsave(plot = .,
         filename = 'output/03_scratch_cells.png',
         width = 8, height = 5)
