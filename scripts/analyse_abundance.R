#---------------------------------------------------------------------------#
# Analysis of abundance data from the location and colonisation experiments #
#---------------------------------------------------------------------------#

# what the script does
# recreates Figure 2 and Figure 3

# load packages
library(tidyverse)
library(lme4)
library(flextable)
library(officer)

#---------------------------------------------------------------------------------#
# LOCATION EXPERIMENT: does abundance differ across sites and types of plastic ####
#---------------------------------------------------------------------------------#

# load in data

# 7-week data across three sites
d <- read.csv('data/location_experiment_abundance.csv')

# calculate means
d_sum <- group_by(d, site, plastic, type, type2) %>%
  summarise(abundance_cor = mean(abundance_cor), .groups = 'drop')

# make Figure 2
ggplot(d, aes(site, abundance_cor)) +
  #MicrobioUoE::geom_pretty_boxplot(fill = 'black', col = 'black') +
  facet_wrap(~type2, labeller = label_parsed) +
  geom_point(aes(col = plastic), position = position_dodge(width = 0.7), size = 4.5, data = d_sum) +
  geom_point(aes(col = plastic), position = position_jitterdodge(dodge.width = 0.7, jitter.width = 0.05), fill = 'white', shape = 21, size = 2) +
  theme_bw(base_size = 14) +
  scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x, n = 4),
                labels = scales::trans_format("log10", scales::math_format(10^.x)),
                limits = c(10^0, 10^4),
                minor_breaks = NULL) +
  labs(y = expression(log[10]~Abundance~(cm^2)),
       x = 'Location') +
  theme(strip.background = element_blank(),
        strip.text = element_text(hjust = 0)) +
  scale_color_brewer('Polymer\ntype', type = 'qual', palette = 6)

ggsave('figures/Figure_2.pdf', last_plot(), width = 13, height = 4)
ggsave('figures/Figure_2.png', last_plot(), width = 13, height = 4)

# actually do the stats

# analysis on all culturable bacteria
model_culture <- lm(log_abundance ~ site * plastic, filter(d, type == 'culturable'), na.action = 'na.fail')
model_culture2 <- lm(log_abundance ~ site + plastic, filter(d, type == 'culturable'), na.action = 'na.fail')
model_culture3 <- lm(log_abundance ~ site, filter(d, type == 'culturable'), na.action = 'na.fail')
model_culture4 <- lm(log_abundance ~ plastic, filter(d, type == 'culturable'), na.action = 'na.fail')
anova(model_culture, model_culture2)
anova(model_culture2, model_culture3)
anova(model_culture2, model_culture4)
# model_culture2 is the best model

# get contrasts out for different sites and different plastics
culture_site <- emmeans::emmeans(model_culture2, pairwise~site) %>%
  .$contrasts %>%
  data.frame() %>%
  mutate(across(where(is.numeric), round, 3)) %>%
  mutate_all(as.character)
culture_plastic <- emmeans::emmeans(model_culture2, pairwise~plastic) %>%
  .$contrasts %>%
  data.frame() %>%
  mutate(across(where(is.numeric), round, 3)) %>%
  mutate_all(as.character)
culture <- bind_rows(culture_site, culture_plastic) %>%
  mutate(type = 'LB media')

# analysis of just the coliforms 
model_coli <- lm(log_abundance ~ site * plastic, filter(d, type == 'coliform'), na.action = 'na.fail')
model_coli2 <- lm(log_abundance ~ site + plastic, filter(d, type == 'coliform'), na.action = 'na.fail')
model_coli3 <- lm(log_abundance ~ site, filter(d, type == 'coliform'), na.action = 'na.fail')
model_coli4 <- lm(log_abundance ~ plastic, filter(d, type == 'coliform'), na.action = 'na.fail')
anova(model_coli, model_coli2)
anova(model_coli2, model_coli3)
anova(model_coli2, model_coli4)
# model_coli2 is best

# get contrasts out for different sites and different plastics
coli_site <- emmeans::emmeans(model_coli2, pairwise~site) %>%
  .$contrasts %>%
  data.frame() %>%
  mutate(across(where(is.numeric), round, 3)) %>%
  mutate_all(as.character)
coli_plastic <- emmeans::emmeans(model_coli2, pairwise~plastic) %>%
  .$contrasts %>%
  data.frame() %>%
  mutate(across(where(is.numeric), round, 3)) %>%
  mutate_all(as.character)
coli <- bind_rows(coli_site, coli_plastic) %>%
  mutate(type = 'coliform')

# analysis of just E coli 
model_ecoli <- lm(log_abundance ~ site * plastic, filter(d, type == 'ecoli'), na.action = 'na.fail')
model_ecoli2 <- lm(log_abundance ~ site + plastic, filter(d, type == 'ecoli'), na.action = 'na.fail')
model_ecoli3 <- lm(log_abundance ~ site, filter(d, type == 'ecoli'), na.action = 'na.fail')
model_ecoli4 <- lm(log_abundance ~ plastic, filter(d, type == 'ecoli'), na.action = 'na.fail')
anova(model_ecoli, model_ecoli2)
anova(model_ecoli2, model_ecoli3)
anova(model_ecoli2, model_ecoli4)
# model_ecoli2 is best

# get contrasts out for different sites and different plastics
ecoli_site <- emmeans::emmeans(model_ecoli2, pairwise~site) %>%
  .$contrasts %>%
  data.frame() %>%
  mutate(across(where(is.numeric), round, 3)) %>%
  mutate_all(as.character)
ecoli_plastic <- emmeans::emmeans(model_ecoli2, pairwise~plastic) %>%
  .$contrasts %>%
  data.frame() %>%
  mutate(across(where(is.numeric), round, 3)) %>%
  mutate_all(as.character)
ecoli <- bind_rows(coli_site, coli_plastic) %>%
  mutate(type = 'E.coli')

# make reproducible table
d_table <- bind_rows(select(culture, type, contrast, everything()),
                     select(coli, type, contrast, everything()),
                     select(ecoli, type, contrast, everything())) %>%
  mutate(p.value = ifelse(p.value == '0', '<0.001', p.value)) %>%
  mutate(contrast = gsub('\n', ' ', contrast))

table <- flextable(d_table) %>%
  align(align = 'center', part = 'all') %>%
  align(align = 'left', j = 'contrast', part = 'all') %>%
  set_header_labels(t.ratio = "t ratio",
                    p.value = "p value",
                    df = 'd.f.') %>%
  italic(j = 'df', part = 'header') %>%
  font(fontname = 'Times', part = 'all') %>%
  fontsize(size = 14, part = 'all') %>%
  bold(i = ~ p.value < 0.05, j = ~ p.value) %>%
  merge_v(., j = ~type) %>%
  valign(valign = 'top', j = 1, part = 'body') %>%
  hline(i = c(13, 26), border = fp_border_default(), part = 'body') %>%
  hline(i = c(3, 16, 29), j = 2:7, border = fp_border_default(), part = 'body') %>%
  fix_border_issues() %>%
  italic(i = 27, j =1) %>%
  vline(j = 1, part = 'body') %>%
  autofit()

# get a png from the html file with webshot
save_as_image(table, 'tables/Table_S1.png', zoom = 3, webshot = 'webshot2')

#----------------------------------------------------------------------------------------#
# Colonisation EXPERIMENT: does abundance change through time with different plastics ####
#----------------------------------------------------------------------------------------#

d2 <- read.csv('data/colonisation_experiment_abundance.csv')

d2_sum <- group_by(d2, plastic, type, week, type2) %>%
  summarise(abundance_cor = mean(abundance_cor), .groups = 'drop')

# do quick plot to look at relationship between abundance x week x plastic
p1 <- ggplot(d2, aes(week, abundance_cor)) +
  facet_wrap(~type2, labeller = label_parsed) +
  geom_point(aes(col = plastic), size = 5, data = d2_sum) +
  geom_point(aes(col = plastic), position = position_jitter(width = 0.1, height = 0), fill = 'white', shape = 21, size = 2) +
  theme_bw(base_size = 14) +
  scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x, n = 4),
                labels = scales::trans_format("log10", scales::math_format(10^.x)),
                limits = c(10^-0.25, 10^4),
                minor_breaks = NULL) +
  labs(y = 'log10 Abundance',
       x = 'Week') +
  theme(strip.background = element_blank(),
        strip.text = element_text(hjust = 0))
p1

# do the stats separately for each organism type

# only culturable bacteria
model_culture <- lm(log_abundance ~ week * plastic, filter(d2, type == 'culturable'), na.action = 'na.fail')
model_culture2 <- lm(log_abundance ~ week + plastic, filter(d2, type == 'culturable'), na.action = 'na.fail')
model_culture3 <- lm(log_abundance ~ week, filter(d2, type == 'culturable'), na.action = 'na.fail')
model_culture4 <- lm(log_abundance ~ plastic, filter(d2, type == 'culturable'), na.action = 'na.fail')
model_culture5 <- lm(log_abundance ~ 1, filter(d2, type == 'culturable'), na.action = 'na.fail')

anova(model_culture, model_culture2)
anova(model_culture2, model_culture3)
anova(model_culture2, model_culture4)
anova(model_culture4, model_culture5)
anova(model_culture3, model_culture5)
# model_culture 3 is the best - only includes time

# first only coliforms
model_coli <- lm(log_abundance ~ week * plastic, filter(d2, type == 'coliform'), na.action = 'na.fail')
model_coli2 <- lm(log_abundance ~ week + plastic, filter(d2, type == 'coliform'), na.action = 'na.fail')
model_coli3 <- lm(log_abundance ~ week, filter(d2, type == 'coliform'), na.action = 'na.fail')
model_coli4 <- lm(log_abundance ~ plastic, filter(d2, type == 'coliform'), na.action = 'na.fail')
model_coli5 <- lm(log_abundance ~ 1, filter(d2, type == 'coliform'), na.action = 'na.fail')

anova(model_coli, model_coli2)
anova(model_coli2, model_coli3)
anova(model_coli2, model_coli4)
anova(model_coli4, model_coli5)
anova(model_coli3, model_coli5)
# model_coli3 is the best - only includes time

# only E.coli
model_ecoli <- lm(log_abundance ~ week * plastic, filter(d2, type == 'ecoli'), na.action = 'na.fail')
model_ecoli2 <- lm(log_abundance ~ week + plastic, filter(d2, type == 'ecoli'), na.action = 'na.fail')

anova(model_ecoli, model_ecoli2)
# mode_ecoli is the best

# look at emtrends to see what is going on with the slopes between plastics
emmeans::emtrends(model_ecoli, pairwise~1, var = 'week')
emmeans::emtrends(model_ecoli, pairwise~plastic, var = 'week')

# make Figure 3
preds <- bind_rows(broom::augment(model_ecoli) %>% mutate(type = 'ecoli'),
                   broom::augment(model_coli3) %>% mutate(type = 'coliform'),
                   broom::augment(model_culture3) %>% mutate(type = 'culturable')) %>%
  mutate(type = forcats::fct_relevel(type, 'culturable', 'coliform', 'ecoli'),
         type2 = forcats::fct_recode(type, '(a)~LB~medium' = 'culturable', '(b)~coliform' = 'coliform', '(c)~italic("E.coli")' = 'ecoli'))
  

ggplot(d2, aes(week, abundance_cor)) +
  facet_wrap(~type2, labeller = label_parsed) +
  geom_line(aes(week, 10^.fitted, col = plastic), filter(preds, type == 'ecoli')) +
  geom_line(aes(week, 10^.fitted), col = 'black', filter(preds, type != 'ecoli')) +
  geom_point(aes(col = plastic), size = 5, data = d2_sum) +
  geom_point(aes(col = plastic), position = position_jitter(width = 0.1, height = 0), fill = 'white', shape = 21, size = 2) +
  theme_bw(base_size = 14) +
  scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x, n = 4),
                labels = scales::trans_format("log10", scales::math_format(10^.x)),
                limits = c(10^-0.25, 10^4.1),
                minor_breaks = NULL) +
  labs(y = expression(log[10]~Abundance~(cm^2)),
       x = 'Week') +
  theme(strip.background = element_blank(),
        strip.text = element_text(hjust = 0)) +
  scale_color_brewer('Polymer\ntype', type = 'qual', palette = 6)

ggsave('figures/Figure_3.pdf', last_plot(), width = 13, height = 4)
ggsave('figures/Figure_3.png', last_plot(), width = 13, height = 4)

