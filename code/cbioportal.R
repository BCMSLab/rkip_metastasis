# load required libraries
library(tidyverse)
library(readxl)

# load data
type <- read_excel('data/cbioportal/cbioportal_breast.xlsx', 1)
type %>%
  select(dataset, cancer_type = `Cancer Type Detailed`,
         alteration_type = `Alteration Type`,
         count = `Alteration Count`) %>%
  ggplot(aes(x = cancer_type, y = count)) +
  geom_col() +
  facet_wrap(~dataset, scales = 'free')

# replace absolute count with percent in the first study

cna <- read_excel('data/cbioportal/cbioportal_breast.xlsx', 2)


cna %>%
  select(cna = `PEBP1: Putative copy-number alterations from DNAcopy.`,
         expr = `PEBP1: mRNA expression (microarray)`,
         everything()) %>%
  group_by(dataset, cna) %>%
  mutate(expr = scale(expr),
         ave = median(expr),
            sd = sd(expr)) %>%
  ggplot(aes(x = cna, y = expr)) +
  geom_boxplot() +
  facet_wrap(~dataset) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
  

 receptors <- read_excel('data/cbioportal/cbioportal_breast.xlsx', 3)
head(receptors)

receptors %>%
  mutate(Status = tolower(Status)) %>%
  select(expr = `PEBP1: mRNA expression (microarray)`, everything()) %>%
  group_by(dataset, receptor, Status) %>%
  summarise(ave = mean(expr)) %>%
  ggplot(aes(x = receptor, y = log2(ave))) +
  geom_col() +
  facet_wrap(~dataset, scales = 'free_y')

alterations <- read_excel('data/cbioportal/cbioportal_breast.xlsx', 4)
alterations %>%
  group_by(dataset, Group) %>%
  summarise(ave = mean(Expression)) %>%
  ggplot(aes(x = Group, y = log2(ave))) +
  geom_col() +
  facet_wrap(~dataset)
