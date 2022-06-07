#--------------------------------------------------------------#
# look at alpha diversity and pielous evenness of the 16S data #
#--------------------------------------------------------------#

# what this script does
# recreates Figure 4

# load packages ####
library(phyloseq)
library(tidyverse)
library(patchwork)
library(palettetown)
library(flextable)
library(officer)

# load rarefied data for the diversity analysis
ps <- readRDS('data/phyloseq_16s_rarefied_nocyano_nomito.rds')

# show available ranks in the dataset
rank_names(ps)

# look at the number of reads per sample
sample_sums(ps)

# alpha diversity estimates ####

# prune OTUs that are not present in any of the samples
ps_sub <- prune_taxa(taxa_sums(ps) > 0, ps)

# metadata #
m <- sample_data(ps_sub) %>%
  data.frame() %>%
  janitor::clean_names()

# calculate diversity measures of each sample
a_div <- estimate_richness(ps_sub, measures = c('Shannon', 'Observed')) %>%
  mutate(., sample_id = row.names(.)) %>%
  # calculate pielou's evenness
  mutate(., pielou = Shannon / log(Observed)) %>%
  janitor::clean_names() %>%
  # join metadata in
  left_join(., m) %>%
  # rename week 4 week 5
  mutate(., id2 = paste(plastic, rep, sep = '_'),
         week = ifelse(week == '4', '5', week))

# plot these measures
head(a_div)

# calculate mean for each week by plastic combination
a_div_summary <- group_by(a_div, week, plastic) %>%
  summarise(observed = mean(observed),
            pielou = mean(pielou),
            .groups = 'drop')

#-------------------------------------------#
# look at differences in alpha diversity ####
#-------------------------------------------#

mod_diversity <- lm(log10(observed) ~ week*plastic, a_div)
mod_diversity2 <- lm(log10(observed) ~ week + plastic, a_div)
mod_diversity3 <- lm(log10(observed) ~ plastic, a_div)
mod_diversity4 <- lm(log10(observed) ~ week, a_div)
mod_diversity5 <- lm(log10(observed) ~ 1, a_div)
anova(mod_diversity, mod_diversity2)
anova(mod_diversity2, mod_diversity3)
anova(mod_diversity2, mod_diversity4)
anova(mod_diversity4, mod_diversity5)
anova(mod_diversity3, mod_diversity5)
# mod_diversity 4 is the best

# create pairwise contrast table
table <- emmeans::emmeans(mod_diversity4, pairwise ~ week)$contrasts %>%
  data.frame() %>%
  mutate(across(estimate:t.ratio, round, 3),
         p.value = round(p.value, 4),
         contrast = gsub('week', '', contrast)) %>%
  flextable() %>%
  align(align = 'center', part = 'all') %>%
  align(align = 'left', part = 'all', j = 2) %>%
  set_header_labels(df = 'd.f',
                    p.value = 'p value',
                    t.ratio = 't ratio') %>%
  italic(., italic = TRUE, part = "header", j = c(4)) %>%
  font(fontname = 'Times', part = 'all') %>%
  fontsize(size = 14, part = 'all') %>%
  bold(i = ~ p.value < 0.05, j = ~ p.value) %>%
  fix_border_issues() %>%
  autofit()

# get a png from the html file with webshot
save_as_image(table, 'tables/Table_S2.png', zoom = 3, webshot = 'webshot2')

# plot data and model predictions
preds_diversity <- emmeans::emmeans(mod_diversity4, pairwise~week, type = 'response') %>%
  .$emmeans %>%
  data.frame() %>%
  mutate(week = as.character(week))

p_div <- ggplot(a_div) +
  geom_linerange(aes(as.numeric(week), ymin = lower.CL, ymax = upper.CL), preds_diversity) +
  geom_point(aes(as.numeric(week), observed, col = week), shape = 21, fill = 'white', position = position_jitter(width = 0.05, height = 0)) +
  geom_line(aes(as.numeric(week), response), preds_diversity) +
  geom_point(aes(as.numeric(week), response, col = week), preds_diversity, size = 4) +
  theme_bw(base_size = 14) +
  theme(strip.background = element_blank(),
        strip.text = element_text(hjust = 0),
        legend.position = 'none') +
  labs(y = "Number of unique ASVs",
       x = 'Week',
       title = '(a)') +
  scale_color_poke(pokemon = 'charizard') +
  ylim(c(0, 2000))

#------------------------------------#
# look at differences in evenness ####
#------------------------------------#

mod_even <- lm(pielou ~ plastic*week, a_div)
mod_even2 <- lm(pielou ~ plastic + week, a_div)
mod_even3 <- lm(pielou ~ plastic, a_div)
mod_even4 <- lm(pielou ~ week, a_div)
mod_even5 <- lm(pielou ~ 1, a_div)
anova(mod_even, mod_even2)
anova(mod_even2, mod_even3)
anova(mod_even2, mod_even4)
anova(mod_even4, mod_even5)
anova(mod_even3, mod_even5)
# mod_even4 is best

table <- emmeans::emmeans(mod_even4, pairwise ~ week)$contrasts %>%
  data.frame() %>%
  mutate(across(estimate:t.ratio, round, 3),
         p.value = round(p.value, 4),
         contrast = gsub('week', '', contrast)) %>%
  flextable() %>%
  align(align = 'center', part = 'all') %>%
  align(align = 'left', part = 'all', j = 2) %>%
  set_header_labels(df = 'd.f',
                    p.value = 'p value',
                    t.ratio = 't ratio') %>%
  italic(., italic = TRUE, part = "header", j = c(4)) %>%
  font(fontname = 'Times', part = 'all') %>%
  fontsize(size = 14, part = 'all') %>%
  bold(i = ~ p.value < 0.05, j = ~ p.value) %>%
  fix_border_issues() %>%
  autofit()

# get a png from the html file with webshot
save_as_image(table, 'tables/Table_S3.png', zoom = 3, webshot = 'webshot2')

# plot data and model predictions
preds_even <- emmeans::emmeans(mod_even4, pairwise~week) %>%
  .$emmeans %>%
  data.frame() %>%
  mutate(week = as.character(week))

p_even <- ggplot(a_div) +
  geom_linerange(aes(as.numeric(week), ymin = lower.CL, ymax = upper.CL), preds_even) +
  geom_point(aes(as.numeric(week), pielou, col = week), shape = 21, fill = 'white', position = position_jitter(width = 0.05, height = 0)) +
  geom_line(aes(as.numeric(week), emmean), preds_even) +
  geom_point(aes(as.numeric(week), emmean, col = week), preds_even, size = 4)  +
  theme_bw(base_size = 14) +
  theme(strip.background = element_blank(),
        strip.text = element_text(hjust = 0),
        legend.position = 'none') +
  labs(y = "Pielou's evenness",
       x = 'Week',
       title = '(b)') +
  ylim(c(0.5,1)) +
  scale_color_poke(pokemon = 'charizard')

#--------------------#
# create Figure 4 ####
#--------------------#

p_div + p_even

ggsave('figures/Figure_4.pdf', last_plot(), height = 3.5, width = 8)
ggsave('figures/Figure_4.png', last_plot(), height = 3.5, width = 8)
