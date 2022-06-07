#-------------------------------------------------------------------------------------------#
# look at differences in community composition of biofilms through time and across plastics #
#-------------------------------------------------------------------------------------------#

# what this script does
# recreates Figure 5, Table S4 and Figure S2, Figure S3 and Figure S4.

# clean workspace
rm(list = ls())

# load packages ####
library(phyloseq)
library(vegan)
library(patchwork)
library(flextable)
library(tidyverse)
library(palettetown)
library(officer)

# source extra functions
source('scripts/extra_functions.R')

# load data - do not use the rarefied data
ps <- readRDS('data/phyloseq_16s_nocyano_nomito.rds')

sort(sample_sums(ps))
rank_names(ps) 

# transform counts to relative abundances for ordination
ps_prop <- transform_sample_counts(ps, function(x){x / sum(x)})

# wrangle the metadata
d_samp <- data.frame(sample_data(ps_prop))
d_samp <- mutate(d_samp, week_fac = as.factor(paste('w', week, sep = '_')),
                 plastic_fac = as.factor(plastic),
                 rep_fac = as.factor(rep)) %>%
  unite(., 'id', week_fac, plastic_fac, sep = ':', remove = FALSE)

# calculate distance matrix - use weighted Unifrac
ps_wunifrac <- phyloseq::distance(ps_prop, method = 'wunifrac')

#--------------------#
# create Figure 5 ####
#--------------------#

# run a betadisper
mod_betadisper <- betadisper(ps_wunifrac, d_samp$id)

# grab centroids and other data
d_fig <- get_betadisper_data(mod_betadisper)

# combine centroid and eigenvector dataframes for plotting
betadisper_lines <- merge(select(d_fig$centroids, group, PCoA1, PCoA2), select(d_fig$eigenvector, group, PCoA1, PCoA2), by = c('group'))

# add distances to eigenvector and lines data
betadisper_lines <- mutate(betadisper_lines, distances = dist_between_points(PCoA1.x, PCoA2.x, PCoA1.y, PCoA2.y))
d_fig$eigenvector$distances <- d_fig$distances$distances

# split up group into week and plastic context
betadisper_lines <- separate(betadisper_lines, group, c('week', 'plastic'), sep =':')
d_fig$centroids <- separate(d_fig$centroids, group, c('week', 'plastic'), sep =':')
d_fig$eigenvector <- separate(d_fig$eigenvector, group, c('week', 'plastic'), sep =':') 

d_fig$eigenvalue <- mutate(d_fig$eigenvalue, percent = eig/sum(eig)*100)

# get correct eigenvalues by applying a correction
correct_eigenvalues <- ape::pcoa(ps_wunifrac, correction = 'cailliez', d_samp$id) %>%
  .$values %>%
  pull(Rel_corr_eig)
correct_eigenvalues[1:2]

# add replicate into eigenvector
d_fig$eigenvector <- bind_cols(d_fig$eigenvector, select(d_samp, rep)) %>%
  mutate(id = paste(plastic, rep, sep = '_'))

# calculate limits to apply across each plot
min(d_fig$eigenvector$PCoA1)
min(d_fig$eigenvector$PCoA2)
max(d_fig$eigenvector$PCoA1)
max(d_fig$eigenvector$PCoA22)

# HDPE
p1 <- ggplot() +
  geom_point(aes(PCoA1, PCoA2, col = week), filter(d_fig$eigenvector, plastic == 'HDPE'), show.legend = FALSE) +
  geom_line(aes(PCoA1, PCoA2), filter(d_fig$centroids, plastic == 'HDPE')) +
  geom_point(aes(PCoA1, PCoA2, col = week), filter(d_fig$centroids, plastic == 'HDPE'), size = 4, show.legend = FALSE) +
  scale_color_poke(pokemon = 'charizard') +
  theme_bw(base_size = 14) +
  theme(strip.background = element_blank(),
        strip.text = element_text(hjust = 0, size = 12),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank()) +
  labs(x = 'Axis 1 (20.3%)',
       y = 'Axis 2 (12.9%)',
       title = '(a) HDPE') +
  xlim(c(-0.09, 0.3)) +
  ylim(c(-0.09, 0.135))

# LDPE
p2 <- ggplot() +
  geom_point(aes(PCoA1, PCoA2, col = week), filter(d_fig$eigenvector, plastic == 'LDPE'), show.legend = FALSE) +
  #ggrepel::geom_label_repel(aes(PCoA1, PCoA2, col = week, label = id), data = to_label) +
  geom_line(aes(PCoA1, PCoA2), filter(d_fig$centroids, plastic == 'LDPE')) +
  geom_point(aes(PCoA1, PCoA2, col = week), filter(d_fig$centroids, plastic == 'LDPE'), size = 4, show.legend = FALSE) +
  scale_color_poke(pokemon = 'charizard') +
  theme_bw(base_size = 14) +
  theme(strip.background = element_blank(),
        strip.text = element_text(hjust = 0, size = 12),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank()) +
  labs(x = 'Axis 1 (20.3%)',
       y = 'Axis 2 (12.9%)',
       title = '(b) LDPE') +
  xlim(c(-0.09, 0.3)) +
  ylim(c(-0.09, 0.135))

# PET
p3 <- ggplot() +
  geom_point(aes(PCoA1, PCoA2, col = week), filter(d_fig$eigenvector, plastic == 'PET'), show.legend = FALSE) +
  #ggrepel::geom_label_repel(aes(PCoA1, PCoA2, col = week, label = id), data = to_label) +
  geom_line(aes(PCoA1, PCoA2), filter(d_fig$centroids, plastic == 'PET')) +
  geom_point(aes(PCoA1, PCoA2, col = week), filter(d_fig$centroids, plastic == 'PET'), size = 4, show.legend = FALSE) +
  scale_color_poke(pokemon = 'charizard') +
  theme_bw(base_size = 14) +
  theme(strip.background = element_blank(),
        strip.text = element_text(hjust = 0, size = 12),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank()) +
  labs(x = 'Axis 1 (20.3%)',
       y = 'Axis 2 (12.9%)',
       title = '(c) PET') +
  xlim(c(-0.09, 0.3)) +
  ylim(c(-0.09, 0.135))

# PP
p4 <- ggplot() +
  geom_point(aes(PCoA1, PCoA2, col = week), filter(d_fig$eigenvector, plastic == 'PP'), show.legend = FALSE) +
  #ggrepel::geom_label_repel(aes(PCoA1, PCoA2, col = week, label = id), data = to_label) +
  geom_line(aes(PCoA1, PCoA2), filter(d_fig$centroids, plastic == 'PP')) +
  geom_point(aes(PCoA1, PCoA2, col = week), filter(d_fig$centroids, plastic == 'PP'), size = 4, show.legend = FALSE) +
  scale_color_poke(pokemon = 'charizard') +
  theme_bw(base_size = 14) +
  theme(strip.background = element_blank(),
        strip.text = element_text(hjust = 0, size = 12),
        axis.title.x = element_blank()) +
  labs(x = 'Axis 1 (20.3%)',
       y = 'Axis 2 (12.9%)',
       title = '(d) PP') +
  xlim(c(-0.09, 0.3)) +
  ylim(c(-0.09, 0.135))

# PVC
p5 <- ggplot() +
  geom_point(aes(PCoA1, PCoA2, col = week), filter(d_fig$eigenvector, plastic == 'PVC'), show.legend = FALSE) +
  #ggrepel::geom_label_repel(aes(PCoA1, PCoA2, col = week, label = id), data = to_label) +
  geom_line(aes(PCoA1, PCoA2), filter(d_fig$centroids, plastic == 'PVC')) +
  geom_point(aes(PCoA1, PCoA2, col = week), filter(d_fig$centroids, plastic == 'PVC'), size = 4, show.legend = FALSE) +
  scale_color_poke(pokemon = 'charizard') +
  theme_bw(base_size = 14) +
  theme(strip.background = element_blank(),
        strip.text = element_text(hjust = 0, size = 12),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank()) +
  labs(x = 'Axis 1 (20.3%)',
       y = 'Axis 2 (12.9%)',
       title = '(e) PVC') +
  xlim(c(-0.09, 0.3)) +
  ylim(c(-0.09, 0.135))

# plot PCoA axis 1
p6 <- ggplot(d_fig$eigenvector, aes(week, PCoA1, col = week, fill = week), show.legend = FALSE) +
  geom_hline(aes(yintercept = 0)) +
  MicrobioUoE::geom_pretty_boxplot() +
  geom_point(shape = 21, fill = 'white', size = 3, position = position_jitter(width = 0.1))  +
  theme_bw(base_size = 14, base_family = 'Helvetica') +
  ylab('PCoA Axis 1 (20.3%)') +
  xlab('') +
  theme(legend.position = 'none') +
  scale_x_discrete(labels = c('Week 1', 'Week 2', 'Week 3', 'Week 5')) +
  ggtitle('(f)') +
  scale_color_poke(pokemon = 'charizard') +
  scale_fill_poke(pokemon = 'charizard')

{p1 + p2 + p3} / {p4 + p5 + p6} + plot_layout(nrow = 2)

ggsave('figures/Figure_5.pdf', last_plot(), height = 7, width = 12)
ggsave('figures/Figure_5.png', last_plot(), height = 7, width = 12)

#-------------------------------------#
# make plot for graphical abstract ####
#-------------------------------------#

graph_abstract_1 <- ggplot() +
  geom_segment(aes(x = PCoA1.x, y = PCoA2.x, yend = PCoA2.y, xend = PCoA1.y, group = row.names(filter(betadisper_lines, week == 'w_4')), col = plastic), filter(betadisper_lines, week == 'w_4'), alpha  = 0.5) +
  geom_point(aes(PCoA1, PCoA2, col = plastic), fill = 'white', filter(d_fig$eigenvector, week == 'w_4'), show.legend = FALSE, shape = 21) +
  geom_point(aes(PCoA1, PCoA2, col = plastic), filter(d_fig$centroids, week == 'w_4'), size = 4) +
  theme_bw(base_size = 12) +
  theme(strip.background = element_blank(),
        strip.text = element_text(hjust = 0, size = 12),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank()) +
  labs(x = 'Axis 1',
       y = 'Axis 2') +
  scale_color_brewer('Plastic', type = 'qual', palette = 6)

# grab legend
legend <- cowplot::get_legend(graph_abstract_1)
cowplot::plot_grid(legend)
ggsave('figures/graphical_abstract_legend.pdf', height = 6, width = 6)

# make plot without legend
graph_abstract_1 +
  theme(legend.position = 'none')

ggsave('figures/pcoa_graphical_abstract.pdf', height = 2, width = 2.5)

#--------------------#
# formal analyses ####
#--------------------#

# first look for differences in position of centroids
# Do this using permanovas
mod_perm <- vegan::adonis2(ps_wunifrac ~ week_fac*plastic_fac, data = d_samp, n_perm = 9999)
mod_perm
mod_perm2 <- vegan::adonis2(ps_wunifrac ~ week_fac + plastic_fac, data = d_samp, n_perm = 9999)
# week makes a difference and plastic makes a very small difference, but only when it week is present
mod_perm4 <- vegan::adonis2(ps_wunifrac ~ plastic_fac, data = d_samp, n_perm = 9999)

# compare plastic composition within each week
dist_all = list(dm_loaded = ps_wunifrac, map_loaded = d_samp)
dist_week1 <- mctoolsr::filter_dm(dist_all, 'week_fac', keep_vals = 'w_1')
mult_comp_plastic_w1 <- calc_pairwise_permanovas(dist_week1$dm_loaded, dist_week1$map_loaded, 'plastic_fac', n_perm = 9999)
dist_week2 <- mctoolsr::filter_dm(dist_all, 'week_fac', keep_vals = 'w_2')
mult_comp_plastic_w2 <- calc_pairwise_permanovas(dist_week2$dm_loaded, dist_week2$map_loaded, 'plastic_fac', n_perm = 9999)
dist_week3 <- mctoolsr::filter_dm(dist_all, 'week_fac', keep_vals = 'w_3')
mult_comp_plastic_w3 <- calc_pairwise_permanovas(dist_week3$dm_loaded, dist_week3$map_loaded, 'plastic_fac', n_perm = 9999)
dist_week4 <- mctoolsr::filter_dm(dist_all, 'week_fac', keep_vals = 'w_4')
mult_comp_plastic_w4 <- calc_pairwise_permanovas(dist_week4$dm_loaded, dist_week4$map_loaded, 'plastic_fac', n_perm = 9999)
# none of these are significant

pvals <- c(mult_comp_plastic_w1$pval, mult_comp_plastic_w2$pval, mult_comp_plastic_w3$pval, mult_comp_plastic_w4$pval)
padj <- p.adjust(pvals, method = 'fdr')
padj[padj <= 0.05]

# check for differences in homogeneity of variance - beta diversity
mod1_dispers <- betadisper(ps_wunifrac, d_samp$id)

# grab the data out for distance from centroid
d_betadiversity <- get_betadisper_data(mod1_dispers)$distances %>%
  separate(group, c('week', 'plastic'), sep = ':')

# test for differences between week and plastic
mod_beta <- lm(distances ~ week*plastic, d_betadiversity)
mod_beta2 <- lm(distances ~ week + plastic, d_betadiversity)
mod_beta3 <- lm(distances ~ week, d_betadiversity)
mod_beta4 <- lm(distances ~ plastic, d_betadiversity)
mod_beta5 <- lm(distances ~ 1, d_betadiversity)
anova(mod_beta, mod_beta2)
anova(mod_beta2, mod_beta3)
anova(mod_beta2, mod_beta4)
anova(mod_beta3, mod_beta5)
anova(mod_beta4, mod_beta5)

#--------------------#
# create Table S4 ####
#--------------------#

data_for_table <- tibble(model_id = c(1,2,3,4,5,5),
                         model = c('distance to centroid ~ week * plastic',
                                   'distance to centroid ~ week + plastic',
                                   'distance to centroid ~ week',
                                   'distance to centroid ~ plastic',
                                   'distance to centroid ~ 1',
                                   'distance to centroid ~ 1'),
                         model_compared = c('', '1,2', '2,3', '2,4', '3,5', '4,5'))
data_for_table2 <- rbind(broom::glance(mod_beta),
                         broom::glance(mod_beta2),
                         broom::glance(mod_beta3),
                         broom::glance(mod_beta4),
                         broom::glance(mod_beta5),
                         broom::glance(mod_beta5)) %>%
  select(df, AIC, logLik)

data_for_table3 <- rbind(broom::tidy(anova(mod_beta, mod_beta2)),
                         broom::tidy(anova(mod_beta2, mod_beta3)),
                         broom::tidy(anova(mod_beta2, mod_beta4)),
                         broom::tidy(anova(mod_beta3, mod_beta5)),
                         broom::tidy(anova(mod_beta4, mod_beta5))) %>%
  filter(!is.na(p.value))

data_for_table <- bind_cols(data_for_table, data_for_table2) %>%
  mutate(., f.value = c(NA, data_for_table3$statistic),
         p.value = c(NA, data_for_table3$p.value)) %>%
  mutate(across(c(AIC,logLik, f.value), round, 2),
         p.value = round(p.value, 3)) %>%
  select(model_id, model, df, AIC, logLik, model_compared, f.value, p.value)
  
# make table using flextable
table <- flextable(data_for_table) %>%
  align(align = 'center', part = 'all') %>%
  align(align = 'left', part = 'all', j = 2) %>%
  set_header_labels(model_id = "",
                    model_compared = "ANOVA\ncomparison",
                    df = 'd.f.',
                    logLik = 'Log Lik',
                    f.value = 'F',
                    p.value = 'p value') %>%
  italic(., italic = TRUE, part = "header", j = c(3,7)) %>%
  font(fontname = 'Times', part = 'all') %>%
  fontsize(size = 14, part = 'all') %>%
  bold(i = ~ p.value < 0.05, j = ~ p.value) %>%
  fix_border_issues() %>%
  autofit()

# get a png from the html file with webshot
save_as_image(table, 'tables/Table_S4.png', zoom = 3, webshot = 'webshot2')

#-----------------------------------------------------------------#
# look at compositional differences are happening through time ####
#-----------------------------------------------------------------#

# first for Class as thats whats generally done in plastics research

# agglomerate at the Class level
ps_class <- speedyseq::tax_glom(ps_prop, 'Class')

# get top classes per sample
d_class <- speedyseq::psmelt(ps_class) %>%
  janitor::clean_names() %>%
  mutate(week = ifelse(week == 4, 5, week),
         week = as.numeric(week))

# calculate mean at each class by week by plastic combination
d_class_mean <- group_by(d_class, week, plastic, class) %>%
  summarise(mean_prop = mean(abundance),
            sd_prop = sd(abundance),
            .groups = 'drop')

# find dominant classes - average proportion > 0.1
d_class_mean_dominant <- group_by(d_class_mean, class) %>%
  filter(mean_prop > 0.1) %>%
  ungroup()

d_class_dominant <- filter(d_class, class %in% unique(d_class_mean_dominant$class))

# how many of the total abundance do these dominant classes cover
group_by(d_class_dominant, sample_id) %>%
  summarise(sum(abundance)) %>%
  pull(2) %>%
  mean()
# 92%

d_class <- filter(d_class, !class %in% unique(d_class_mean_dominant$class))
d_class_mean <- filter(d_class_mean, !class %in% unique(d_class_mean_dominant$class))

# plot
ggplot(d_class) +
  geom_point(aes(week, abundance, col = class), d_class_dominant, position = position_jitter(width = 0.1, height = 0)) +
  geom_line(aes(week, mean_prop, group = class, col = class), d_class_mean_dominant) +
  geom_point(aes(week, mean_prop, col = class), size = 3, d_class_mean_dominant) +
  geom_point(aes(week, abundance, col = class), d_class, col = 'grey', alpha = 0.05, position = position_jitter(width = 0.1, height = 0)) +
  facet_wrap(~plastic, labeller = labeller(plastic = MicrobioUoE::letter_facets)) +
  theme_bw(base_size = 14, base_family = 'Helvetica') +
  ylim(c(0,1)) +
  theme(strip.background = element_blank(),
        strip.text = element_text(hjust = 0),
        legend.position = c(0.85, 0.2),
        legend.text = element_text(face = 'italic')) +
  labs(y = 'Relative abundance',
       x = 'Week') +
  palettetown::scale_color_poke('Class', pokemon = 'bulbasaur')

ggsave('figures/Figure_S2.pdf', last_plot(), height = 7, width = 12)
ggsave('figures/Figure_S2.png', last_plot(), height = 7, width = 12)


# agglomerate at the Order level
ps_order <- speedyseq::tax_glom(ps_prop, 'Order')

# get top orders per sample
d_order <- speedyseq::psmelt(ps_order) %>%
  janitor::clean_names()  %>%
  mutate(week = ifelse(week == 4, 5, week),
         week = as.numeric(week))

# calculate mean for each order by week by plastic combination
d_order_mean <- group_by(d_order, week, plastic, order) %>%
  summarise(mean_prop = mean(abundance),
            sd_prop = sd(abundance),
            .groups = 'drop')

# find the dominant families - defined as any proportion over 0.075
d_order_mean_dominant <- group_by(d_order_mean, order) %>%
  filter(max(mean_prop) > 0.075) %>%
  ungroup()
d_order_dominant <- filter(d_order, order %in% unique(d_order_mean_dominant$order))

# calculate how much of the total abundance our dominant orders encompass
group_by(d_order_dominant, sample_id) %>%
  summarise(sum(abundance)) %>%
  pull(2) %>%
  mean()

d_order <- filter(d_order, !order %in% unique(d_order_mean_dominant$order))
d_order_mean <- filter(d_order_mean, !order %in% unique(d_order_mean_dominant$order))

ggplot(d_order) +
  geom_point(aes(week, abundance, col = order), d_order_dominant, position = position_jitter(width = 0.1, height = 0)) +
  geom_line(aes(week, mean_prop, group = order, col = order), d_order_mean_dominant) +
  geom_point(aes(week, mean_prop, col = order), size = 3, d_order_mean_dominant) +
  geom_point(aes(week, abundance, col = order), d_order, col = 'grey', alpha = 0.05, position = position_jitter(width = 0.1, height = 0)) +
  facet_wrap(~plastic, labeller = labeller(plastic = MicrobioUoE::letter_facets)) +
  theme_bw(base_size = 14, base_family = 'Helvetica') +
  ylim(c(0,0.5)) +
  theme(strip.background = element_blank(),
        strip.text = element_text(hjust = 0),
        legend.position = c(0.85, 0.2),
        legend.text = element_text(face = 'italic')) +
  labs(y = 'Relative abundance',
       x = 'Week') +
  palettetown::scale_color_poke('Order', pokemon = 'squirtle')

ggsave('figures/Figure_S3.pdf', last_plot(), height = 7, width = 12)
ggsave('figures/Figure_S3.png', last_plot(), height = 7, width = 12)

# agglomerate at the Family level
ps_family <- speedyseq::tax_glom(ps_prop, 'Family')

# get top families per sample
d_family <- speedyseq::psmelt(ps_family) %>%
  janitor::clean_names()  %>%
  mutate(week = ifelse(week == 4, 5, week),
         week = as.numeric(week))

# calculate means for each family by week by plastic combination
d_family_mean <- group_by(d_family, week, plastic, family) %>%
  summarise(mean_prop = mean(abundance),
            sd_prop = sd(abundance),
            .groups = 'drop')

# find the dominant families - defined as any proportion over 0.05
d_family_mean_dominant <- group_by(d_family_mean, family) %>%
  filter(max(mean_prop) > 0.05) %>%
  ungroup()
d_family_dominant <- filter(d_family, family %in% unique(d_family_mean_dominant$family))

d_family <- filter(d_family, !family %in% unique(d_family_mean_dominant$family))
d_family_mean <- filter(d_family_mean, !family %in% unique(d_family_mean_dominant$family))

# plot
ggplot(d_family) +
  geom_point(aes(week, abundance, col = family), d_family_dominant, position = position_jitter(width = 0.1, height = 0)) +
  geom_line(aes(week, mean_prop, group = family, col = family), d_family_mean_dominant) +
  geom_point(aes(week, mean_prop, col = family), size = 3, d_family_mean_dominant) +
  geom_point(aes(week, abundance, col = family), d_family, col = 'grey', alpha = 0.05, position = position_jitter(width = 0.1, height = 0)) +
  facet_wrap(~plastic, labeller = labeller(plastic = MicrobioUoE::letter_facets)) +
  theme_bw(base_size = 14, base_family = 'Helvetica') +
  ylim(c(0,0.4)) +
  theme(strip.background = element_blank(),
        strip.text = element_text(hjust = 0),
        legend.position = c(0.85, 0.2),
        legend.text = element_text(face = 'italic')) +
  labs(y = 'Relative abundance',
       x = 'Week') +
  palettetown::scale_color_poke('Family', pokemon = 'spearow')

ggsave('figures/Figure_S4.pdf', last_plot(), height = 7, width = 12)
ggsave('figures/Figure_S4.png', last_plot(), height = 7, width = 12)

