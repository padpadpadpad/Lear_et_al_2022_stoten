#---------------------------------------------------------------------------#
# Analysis of virulence data from the location and colonisation experiments #
#---------------------------------------------------------------------------#

# what the script does
# recreates Figure 6

# load in packages
library(survival)
library(rstanarm) # had to fork the rstanarm package on GitHub and download the survival branch. Opened the R project and clicked Install and Rebuild.
library(bayesplot)
library(tidyverse)
library(tidybayes)
library(patchwork)
library(palettetown)
library(flextable)
library(officer)

#--------------------------------------------------------#
# 1. Analyse survival curves from location experiment ####
#--------------------------------------------------------#

# load in data
d <- read.csv('data/location_experiment_virulence.csv')

# status: 1 = dead, 0 = alive

# run a parametric proportional hazards model
# follows instructions from "Bayesian Survival Analysis Using the rstanarm R Package, 2020, Brilleman et al. arXiv"

# model site and plastic differences
mod_location <- stan_surv(Surv(time_of_death, status) ~ site*plastic + (1|sample),
                        data = d,
                        chains = 3,
                        cores = 3,
                        seed = 42,
                        iter = 3000)

mod_location
summary(mod_location)

# what are the things we want to know...
# 1. does virulence (as measured by hazard ratio) differ overall between sites
# 2. does virulence (as measured by hazard ratio) differ overall between plastics
# 3. Are there wacky differences between plastic * sites

# get a list of the variables in the model
tidybayes::get_variables(mod_location)

# grab just the fixed effects we are interested in
to_plot <- tidybayes::get_variables(mod_location)[1:15]

# check key mcmc trace plots
mcmc_trace(mod_location, pars = to_plot)
# should be fuzzy caterpillars

# extract draws and calculate hazard for each site by plastic combination
params_location <- spread_draws(mod_location, !!!syms(to_plot)) %>%
        janitor::clean_names() %>%
        mutate(a_H = intercept,
               b_H = intercept + site_b,
               c_H = intercept + site_c,
               a_L = intercept + plastic_l,
               b_L = intercept + plastic_l + site_b + site_b_plastic_l,
               c_L = intercept + plastic_l + site_c + site_c_plastic_l,
               a_P = intercept + plastic_p,
               b_P = intercept + plastic_p + site_b + site_b_plastic_p,
               c_P = intercept + plastic_p + site_c + site_c_plastic_p,
               `a_T` = intercept + plastic_t,
               `b_T` = intercept + plastic_t + site_b + site_b_plastic_t,
               `c_T` = intercept + plastic_t + site_c + site_c_plastic_t,
               a_V = intercept + plastic_v,
               b_V = intercept + plastic_v + site_b + site_b_plastic_v,
               c_V = intercept + plastic_v + site_c + site_c_plastic_p,
               # calculate average hazard for each location/site
               a = (a_H + a_L + a_P + a_V + a_T)/5,
               b = (b_H+ b_L+ b_P+ b_V+ b_T)/5,
               c = (c_H + c_L + c_P + c_V + c_T)/5,
               # calculate average hazard for each plastic
               H = (a_H + b_H + c_H)/3,
               L = (a_L + b_L + c_L)/3,
               P = (a_P + b_P + c_P)/3,
               V = (a_V + b_V + c_V)/3,
               `T` = (a_T + b_T + c_T)/3,
               # calculate average hazard across sites and across plastics
               ave_site = (a + b + c)/3,
               ave_plastic = (H + L + P + V + `T`)/5)
# ave_site and ave_plastic are reassuringly the same!

# calculate some hazard ratios
# compare each site and plastic to the average site/plastic. This allows us to see which plastics/sites are the most/least virulent
hazard_ratios = mutate(params_location, hazard_ratio_a = exp(a - ave_site),
                       hazard_ratio_b = exp(b - ave_site),
                       hazard_ratio_c = exp(c - ave_site),
                       hazard_ratio_H = exp(H - ave_plastic),
                       hazard_ratio_T = exp(`T` - ave_plastic),
                       hazard_ratio_L = exp(L - ave_plastic),
                       hazard_ratio_P = exp(P - ave_plastic),
                       hazard_ratio_V = exp(V - ave_plastic))

# grab location hazard ratios
hazard_ratios_location <- select(hazard_ratios, starts_with('hazard')) %>%
        pivot_longer(cols = everything(), names_to = 'variable', values_to = 'estimate') %>%
        group_by(variable) %>%
        median_qi() %>%
        ungroup() %>%
        mutate(var = gsub('hazard_ratio_', '', variable))

# check whether any of them do not cross 0, no significant difference in hazard rate
ggplot(hazard_ratios_location, aes(y = var, x = estimate, xmin = .lower, xmax = .upper)) +
        geom_point() +
        geom_linerange() +
        theme_bw(base_size = 14) +
        geom_vline(aes(xintercept = 1), linetype = 2) +
        labs(y = 'Groups used',
             x = 'Hazard ratio') 

# calculate survival curves #

# for location only

# 1. predict over all random effects
d_preds_site <- select(d, site, plastic, sample) %>%
        distinct() %>%
        mutate(id = 1:n(),
               plastic2 = plastic,
               site2 = site) %>%
        nest_legacy(-c(site2)) %>%
        mutate(., preds = map(data, ~posterior_survfit(mod_location, newdata = .x, times = 0, standardise = TRUE, extrapolate = TRUE))) %>%
        unnest(preds) %>%
        select(site = site2, everything(), -data)

# 2. standardise over all ids - this gives a single survival curve for the whole dataset
d_preds_ave <- posterior_survfit(mod_location, times = 0, standardise = TRUE, extrapolate = TRUE) 

# plot these
p1 <- ggplot() +
        geom_line(aes(time, col = site, y= median), show.legend = FALSE,data = d_preds_site, alpha = 0.6) +
        geom_ribbon(aes(time, ymin = ci_lb, ymax = ci_ub, fill = site), col = NA, alpha = 0.05, show.legend = FALSE, data = d_preds_site) +
        geom_ribbon(aes(time, ymin = ci_lb, ymax = ci_ub), col = NA, alpha = 0.3, show.legend = FALSE, data = d_preds_ave) +
        geom_line(aes(time, y= median), show.legend = FALSE, data = d_preds_ave) +
        ylim(c(0,1)) +
        theme_bw(base_size = 14) +
        labs(y = 'Survival probability',
             x = 'Time (hours)',
             title = '(a) virulence across sites',
             subtitle = 'location experiment') +
        theme(strip.background = element_blank(),
              strip.text = element_text(hjust = 0),
              axis.title.x = element_blank()) +
        scale_color_poke(pokemon = 'oddish', spread = 3) +
        scale_fill_poke(pokemon = 'oddish', spread = 3)

# plot hazard ratios for this plot
p2 <- ggplot(filter(hazard_ratios_location, var %in% c('a', 'b', 'c')), aes(y = var, x = estimate, xmin = .lower, xmax = .upper, col = var)) +
        geom_vline(aes(xintercept = 1), linetype = 2) +
        geom_point(size = 3, show.legend = FALSE) +
        geom_linerange(show.legend = FALSE) +
        theme_bw(base_size = 10) +
        theme(axis.title.y = element_blank()) +
        labs(y = 'Site',
             x = 'Hazard ratio') +
        scale_y_discrete(labels = c(a = 'dockyard', b = 'falmouth wharf', c = 'mylor bank')) +
        xlim(c(0,3.5)) +
        scale_color_poke(pokemon = 'oddish', spread = 3)

panel_1 <- p1 + inset_element(p2, left = 0.01, bottom = 0.01, right = 0.45, top = 0.45)

# for the different plastics

# 1. predict over all random effects
d_preds_plastic <- select(d, site, plastic, sample) %>%
        distinct() %>%
        mutate(id = 1:n(),
               plastic2 = plastic,
               site2 = site) %>%
        nest_legacy(-c(plastic2)) %>%
        mutate(., preds = map(data, ~posterior_survfit(mod_location, newdata = .x, times = 0, standardise = TRUE, extrapolate = TRUE))) %>%
        unnest(preds) %>%
        select(plastic = plastic2, everything(), -data)

# plot these
p3 <- ggplot() +
        geom_line(aes(time, col = plastic, y= median), show.legend = FALSE,data = d_preds_plastic, alpha = 0.6) +
        geom_ribbon(aes(time, ymin = ci_lb, ymax = ci_ub, fill = plastic), col = NA, alpha = 0.05, show.legend = FALSE, data = d_preds_plastic) +
        geom_ribbon(aes(time, ymin = ci_lb, ymax = ci_ub), col = NA, alpha = 0.3, show.legend = FALSE, data = d_preds_ave) +
        geom_line(aes(time, y= median), show.legend = FALSE, data = d_preds_ave) +
        ylim(c(0,1)) +
        theme_bw(base_size = 14) +
        labs(y = 'Survival probability',
             x = 'Time (hours)',
             title = '(b) virulence across plastics',
             subtitle = 'location experiment') +
        theme(strip.background = element_blank(),
              strip.text = element_text(hjust = 0),
              axis.title.x = element_blank()) +
        scale_color_brewer('Plastic', type = 'qual', palette = 6) +
        scale_fill_brewer('Plastic', type = 'qual', palette = 6)

# plot hazard ratios for this plot
p4 <- ggplot(filter(hazard_ratios_location, !var %in% c('a', 'b', 'c')), aes(y = var, x = estimate, xmin = .lower, xmax = .upper, col = var)) +
        geom_vline(aes(xintercept = 1), linetype = 2) +
        geom_point(size = 3, show.legend = FALSE) +
        geom_linerange(show.legend = FALSE) +
        theme_bw(base_size = 10) +
        theme(axis.title.y = element_blank()) +
        labs(y = 'Plastic',
             x = 'Hazard ratio') +
        scale_y_discrete(labels = c(H = 'HDPE', L = 'LDPE', P = 'PP', `T` = 'PET', V = 'PVC')) +
        xlim(c(0,6)) +
        scale_color_brewer('Plastic', type = 'qual', palette = 6)

panel_2 <- p3 + inset_element(p4, left = 0.01, bottom = 0.01, right = 0.45, top = 0.45)

# make the plot used in the graphical abstract
ggplot() +
        geom_line(aes(time, col = plastic, y= median), show.legend = FALSE,data = d_preds_plastic, alpha = 0.6) +
        geom_ribbon(aes(time, ymin = ci_lb, ymax = ci_ub, fill = plastic), col = NA, alpha = 0.05, show.legend = FALSE, data = d_preds_plastic) +
        geom_ribbon(aes(time, ymin = ci_lb, ymax = ci_ub), col = NA, alpha = 0.3, show.legend = FALSE, data = d_preds_ave) +
        geom_line(aes(time, y= median), show.legend = FALSE, data = d_preds_ave) +
        ylim(c(0.4,1)) +
        theme_bw(base_size = 12) +
        labs(y = 'Survival probability',
             x = 'Time') +
        theme(strip.background = element_blank(),
              strip.text = element_text(hjust = 0),
              panel.grid = element_blank(),
              axis.text = element_blank(),
              axis.ticks = element_blank()) +
        scale_color_brewer('Plastic', type = 'qual', palette = 6) +
        scale_fill_brewer('Plastic', type = 'qual', palette = 6)

ggsave('figures/survival_curves_graphical_abstract.pdf', height = 2, width = 2.5)

# ok lets check the number that died etc
# calculate summary stats for the location experiment

# for the different locations
summary_stats1 <- group_by(d, site, status) %>%
        summarise(mean = mean(time_of_death),
                  n = n(), .groups = 'drop') %>%
        pivot_wider(names_from = status, values_from = c(mean, n)) %>%
        mutate(n_0 = replace_na(n_0, 0), prop_died = n_1/(n_0 + n_1)) %>%
        select(site, 
               prop_died,
               time_to_death = mean_1) %>%
        mutate(site = case_when(site == 'A' ~ 'dockyard',
                                       site == 'B' ~ 'falmouth wharf',
                                       site == 'C' ~ 'mylor bank'))

# for the different plastics
summary_stats2 <- group_by(d, plastic, status) %>%
        summarise(mean = mean(time_of_death),
                  n = n(), .groups = 'drop') %>%
        pivot_wider(names_from = status, values_from = c(mean, n)) %>%
        mutate(n_0 = replace_na(n_0, 0), prop_died = n_1/(n_0 + n_1)) %>%
        select(plastic, 
               prop_died,
               time_to_death = mean_1) %>%
        mutate(plastic = case_when(plastic == 'T' ~ 'PET',
                                   plastic =='V' ~ 'PVC',
                                   plastic == 'P' ~ 'PP',
                                   plastic == 'L' ~ 'LDPE',
                                   plastic == 'H' ~ 'HDPE'))

arrange(summary_stats1, site)


#------------------------------------------------------------#
# 1. Analyse survival curves from colonisation experiment ####
#------------------------------------------------------------#

# load in data
d2 <- read.csv('data/colonisation_experiment_virulence.csv') %>%
        mutate(week = as.character(week))

# calculate summary stats 
# remove water from the summary stats as nearly none died at all
summary_stats3 <- filter(d2, plastic != 'W') %>%
        group_by(week, status) %>%
        summarise(mean = mean(time_of_death),
                  n = n(), .groups = 'drop') %>%
        pivot_wider(names_from = status, values_from = c(mean, n)) %>%
        mutate(prop_died = n_1/(n_0 + n_1)) %>%
        select(week,
               prop_died,
               time_to_death = mean_1) %>%
        arrange(week)
summary_stats4 <- filter(d2) %>%
        filter(., week == 5) %>%
        group_by(plastic, status) %>%
        summarise(mean = mean(time_of_death),
                  n = n(), .groups = 'drop') %>%
        pivot_wider(names_from = status, values_from = c(mean, n)) %>%
        mutate(prop_died = n_1/(n_0 + n_1)) %>%
        select(plastic, 
               prop_died,
               time_to_death = mean_1)

# model site and plastic differences
mod_week <- stan_surv(Surv(time_of_death, status) ~ week*plastic + (1|sample),
                          data = d2,
                          chains = 3,
                          cores = 3,
                          seed = 42,
                          iter = 3000)
mod_week

# get a list of the variables in the model
tidybayes::get_variables(mod_week)

# grab just the fixed effect parameters of interest
to_plot <- tidybayes::get_variables(mod_week)[1:24]

# check key mcmc trace plots
mcmc_trace(mod_week, pars = to_plot)
# all fuzzy caterpillars

# show names
spread_draws(mod_week, !!!syms(to_plot)) %>%
        janitor::clean_names() %>% names()

# extract draws for fixed effect parameters
params_week <- spread_draws(mod_week, !!!syms(to_plot)) %>%
        janitor::clean_names() %>%
        # calculate hazard for week 5 only
        mutate(H_5 = intercept + week5,
               L_5 = intercept + plastic_l + week5 + week5_plastic_l,
               T_5 = intercept + plastic_t + week5 + week5_plastic_t,
               V_5 = intercept + plastic_v + week5 + week5_plastic_v,
               W_5 = intercept + plastic_w + week5 + week5_plastic_w,
               P_5 = intercept + plastic_p + week5 + week5_plastic_p,
               ave_plastic = (H_5 + L_5 + T_5 + V_5 + P_5)/5) 

# calculate hazard ratios of each plastic from week 5 against the average hazard at week 5
hazard_ratios_week = mutate(params_week,
                       hazard_ratio_H = exp(H_5 - ave_plastic),
                       hazard_ratio_T = exp(`T_5` - ave_plastic),
                       hazard_ratio_L = exp(L_5 - ave_plastic),
                       hazard_ratio_P = exp(P_5 - ave_plastic),
                       hazard_ratio_V = exp(V_5 - ave_plastic))

# calculate 95% credible intervals of each hazard ratio
hazard_ratios_week <- select(hazard_ratios_week, starts_with('hazard')) %>%
        pivot_longer(cols = everything(), names_to = 'variable', values_to = 'estimate') %>%
        group_by(variable) %>%
        median_qi() %>%
        ungroup() %>%
        mutate(var = gsub('hazard_ratio_', '', variable))

# calculate survival curves

# predict population-level estimates (averaging over replicates)
d_preds_time <- select(d2, week, plastic, sample) %>%
        distinct() %>%
        mutate(id = 1:n(),
               plastic2 = plastic,
               week2 = week) %>%
        nest_legacy(-c(plastic2, week2)) %>%
        mutate(., preds = map(data, ~posterior_survfit(mod_week, newdata = .x, times = 0, standardise = TRUE, extrapolate = TRUE, dynamic = TRUE)))

d_preds_time <- unnest(d_preds_time, preds) %>%
        select(-data) %>%
        dplyr::rename(., plastic = plastic2, week = week2)

head(d_preds_time)

# calculate average survival curve for week 5
d_preds_ave <- select(d2, week, plastic, sample) %>%
        distinct() %>%
        filter(week == 5) %>%
        filter(plastic != 'W') %>%
        posterior_survfit(mod_week, newdata = ., times = 0, standardise = TRUE, extrapolate = TRUE) 

# plot all weeks
all_weeks <- ggplot(d_preds_time, aes(time, col = plastic, fill = plastic, group = plastic)) +
        geom_line(aes(y= median)) +
        geom_ribbon(aes(time, ymin = ci_lb, ymax = ci_ub), col = NA, alpha = 0.1, show.legend = FALSE) +
        facet_wrap(~week, labeller = labeller(week = MicrobioUoE::letter_facets)) +
        ylim(c(0,1)) +
        theme_bw(base_size = 14) +
        labs(
             y = 'Survival probability',
             x = 'Time (hours)') +
        theme(strip.background = element_blank(),
              strip.text = element_text(hjust = 0)) +
        scale_color_brewer('Plastic', type = 'qual', palette = 6) +
        scale_fill_brewer('Plastic', type = 'qual', palette = 6)

# plot just week 5
p5 <- ggplot(filter(d_preds_time, week == '5'), aes()) +
        geom_line(aes(y= median, time, col = plastic, group = plastic), alpha = 0.6, show.legend = FALSE) +
        geom_ribbon(aes(time, ymin = ci_lb, ymax = ci_ub, time, fill = plastic, group = plastic), col = NA, show.legend = FALSE, alpha = 0.05) +
        geom_ribbon(aes(time, ymin = ci_lb, ymax = ci_ub), col = NA, alpha = 0.3, show.legend = FALSE, data = d_preds_ave) +
        geom_line(aes(time, y= median), show.legend = FALSE, data = d_preds_ave) +
        ylim(c(0,1)) +
        theme_bw(base_size = 14) +
        labs(title = '(c) virulence across plastics',
             subtitle = 'colonisation experiment',
             y = 'Survival probability',
             x = 'Time (hours)') +
        theme(legend.position = c(0.1, 0.28)) +
        scale_color_brewer('Plastic', type = 'qual', palette = 6) +
        scale_fill_brewer('Plastic', type = 'qual', palette = 6)

# plot hazard ratios for week 5
p6 <- ggplot(hazard_ratios_week, aes(y = var, x = estimate, xmin = .lower, xmax = .upper, col = var)) +
        geom_vline(aes(xintercept = 1), linetype = 2) +
        geom_point(size = 3, show.legend = FALSE) +
        geom_linerange(show.legend = FALSE) +
        theme_bw(base_size = 10) +
        theme(axis.title.y = element_blank()) +
        labs(y = 'Plastic',
             x = 'Hazard ratio') +
        scale_y_discrete(labels = c(H = 'HDPE', L = 'LDPE', P = 'PP', `T` = 'PET', V = 'PVC')) +
        xlim(c(0,6)) +
        scale_color_brewer('Plastic', type = 'qual', palette = 6)

panel_3 <- p5 + inset_element(p6, left = 0.01, bottom = 0.01, right = 0.45, top = 0.45)

# create Figure 6
survival_curves <- panel_1 + panel_2 + panel_3

ggsave('figures/Figure_6.pdf', survival_curves, width = 16, height = 5)
ggsave('figures/Figure_6.png', survival_curves, width = 16, height = 5)

#------------------------------------------------#
# Calculate hazard ratios for each experiment ####
#------------------------------------------------#

# calculate hazard ratios for each pairwise sample
# also calculate probability that hazard ratio is >/< 1

# calculate hazard ratios for each location/site

# HR sites - get each combination of location a b and c
HR_site <- tibble(site = c('a', 'b', 'c')) %>%
        group_by(site) %>%
        do(data.frame(baseline_site = c('a', 'b', 'c'))) %>%
        filter(site != baseline_site) %>%
        ungroup()

# add empty columns
HR_site<- mutate(HR_site, median = NA, .lower = NA, .upper = NA, prob_below_1 = NA, prob_above_1 = NA)

# for loop to calculate necessary parameters
for(i in 1:nrow(HR_site)){
        temp_site <- select(params_location, HR_site$site[i]) %>% pull()
        temp_baseline_site <- select(params_location, HR_site$baseline_site[i]) %>% pull()
        HR_site$median[i] <- median(exp(temp_site - temp_baseline_site))
        HR_site$.lower[i] <- quantile(exp(temp_site - temp_baseline_site), probs = 0.025)
        HR_site$.upper[i] <- quantile(exp(temp_site - temp_baseline_site), probs = 0.975)
        HR_site$prob_above_1[i] <- exp(temp_site - temp_baseline_site) %>% .[. > 1] %>% length(.)/length(temp_site)
        HR_site$prob_below_1[i] <- exp(temp_site - temp_baseline_site) %>% .[. <= 1] %>% length(.)/length(temp_site)
}

# keep just the corresponding contrasts which are above 1
HR_site <- filter(HR_site, median > 1) %>%
        mutate(site = case_when(site == 'a' ~ 'dockyard',
                                site == 'b' ~ 'falmouth wharf',
                                site == 'c' ~ 'mylor bank'),
               baseline_site = case_when(baseline_site == 'a' ~ 'dockyard',
                                         baseline_site == 'b' ~ 'falmouth wharf',
                                         baseline_site == 'c' ~ 'mylor bank'),
               experiment = 'location experiment')

# calculate hazard ratios for each plastic in the location experiment

# HR plastics from location experiment
HR_plastics <- tibble(plastic = c('T', 'V', 'P', 'L', 'H')) %>%
        group_by(plastic) %>%
        do(data.frame(baseline_plastic = c('T', 'V', 'P', 'L', 'H'))) %>%
        filter(plastic != baseline_plastic) %>%
        ungroup() %>%
        mutate(time = '7 weeks')

# add empty columns for for loop to populate
HR_plastics <- mutate(HR_plastics, median = NA, .lower = NA, .upper = NA, prob_below_1 = NA, prob_above_1 = NA)

# run for loop to calculate parameters for each contrast
for(i in 1:nrow(HR_plastics)){
        temp_plastic <- select(params_location, HR_plastics$plastic[i]) %>% pull()
        temp_baseline_plastic <- select(params_location, HR_plastics$baseline_plastic[i]) %>% pull()
        HR_plastics$median[i] <- median(exp(temp_plastic - temp_baseline_plastic))
        HR_plastics$.lower[i] <- quantile(exp(temp_plastic - temp_baseline_plastic), probs = 0.025)
        HR_plastics$.upper[i] <- quantile(exp(temp_plastic - temp_baseline_plastic), probs = 0.975)
        HR_plastics$prob_above_1[i] <- exp(temp_plastic - temp_baseline_plastic) %>% .[. > 1] %>% length(.)/length(temp_site)
        HR_plastics$prob_below_1[i] <- exp(temp_plastic - temp_baseline_plastic) %>% .[. <= 1] %>% length(.)/length(temp_site)
        
        }

# keep only corresponding contrast which is > 1
HR_plastics <- filter(HR_plastics, median > 1) %>%
        mutate(plastic = case_when(plastic == 'T' ~ 'PET',
                                   plastic =='V' ~ 'PVC',
                                   plastic == 'P' ~ 'PP',
                                   plastic == 'L' ~ 'LDPE',
                                   plastic == 'H' ~ 'HDPE'),
               baseline_plastic = case_when(baseline_plastic == 'T' ~ 'PET',
                                            baseline_plastic =='V' ~ 'PVC',
                                            baseline_plastic == 'P' ~ 'PP',
                                            baseline_plastic == 'L' ~ 'LDPE',
                                            baseline_plastic == 'H' ~ 'HDPE'),
               experiment = 'location experiment')

# calculate hazard ratios for plastics in the colonisation experiment

HR_colonisation <- tibble(plastic = c('T_5', 'V_5', 'P_5', 'L_5', 'H_5')) %>%
        group_by(plastic) %>%
        do(data.frame(baseline_plastic = c('T_5', 'V_5', 'P_5', 'L_5', 'H_5'))) %>%
        filter(plastic != baseline_plastic) %>%
        ungroup()

# add empty columns
HR_colonisation <- mutate(HR_colonisation, median = NA, .lower = NA, .upper = NA, prob_below_1 = NA, prob_above_1 = NA)

# run for loop to calculate parameters for each contrast
for(i in 1:nrow(HR_colonisation)){
        temp_plastic <- select(params_week, HR_colonisation$plastic[i]) %>% pull()
        temp_baseline_plastic <- select(params_week, HR_colonisation$baseline_plastic[i]) %>% pull()
        HR_colonisation$median[i] <- median(exp(temp_plastic - temp_baseline_plastic))
        HR_colonisation$.lower[i] <- quantile(exp(temp_plastic - temp_baseline_plastic), probs = 0.025)
        HR_colonisation$.upper[i] <- quantile(exp(temp_plastic - temp_baseline_plastic), probs = 0.975)
        HR_colonisation$prob_above_1[i] <- exp(temp_plastic - temp_baseline_plastic) %>% .[. > 1] %>% length(.)/length(temp_site)
        HR_colonisation$prob_below_1[i] <- exp(temp_plastic - temp_baseline_plastic) %>% .[. <= 1] %>% length(.)/length(temp_site)
}

# remove _5 from plastic column
HR_colonisation <- mutate(HR_colonisation, across(1:2, function(x) gsub('_5', '', x)),
                      time = '5 weeks')

# rename plastics
HR_colonisation <- mutate(HR_colonisation, plastic = case_when(plastic == 'T' ~ 'PET',
                                   plastic =='V' ~ 'PVC',
                                   plastic == 'P' ~ 'PP',
                                   plastic == 'L' ~ 'LDPE',
                                   plastic == 'H' ~ 'HDPE'),
               baseline_plastic = case_when(baseline_plastic == 'T' ~ 'PET',
                                            baseline_plastic =='V' ~ 'PVC',
                                            baseline_plastic == 'P' ~ 'PP',
                                            baseline_plastic == 'L' ~ 'LDPE',
                                            baseline_plastic == 'H' ~ 'HDPE'),
               experiment = 'colonisation experiment') %>%
        filter(., paste(plastic, baseline_plastic, sep = '') %in% paste(HR_plastics$plastic, HR_plastics$baseline_plastic, sep = ''))

# create table for paper
d_table <- bind_rows(select(HR_site, experiment, first = site, baseline = baseline_site, median, .lower, .upper, prob_above_1),
                   bind_rows(select(HR_plastics, experiment, first = plastic, baseline = baseline_plastic, median, .lower, .upper, prob_above_1),
                   select(HR_colonisation, experiment, first = plastic, baseline = baseline_plastic, median, .lower, .upper, prob_above_1)) %>% arrange(first, baseline)) %>%
        mutate(across(where(is.numeric), round, 2),
               contrast = paste(first, 'vs.', baseline, sep = ' '))
        
# create table
table <- flextable(select(d_table, contrast, experiment, median, .lower, .upper, prob_above_1)) %>%
        set_header_labels(median = 'Median',
                          .lower = "2.5%",
                          .upper = "97.5%",
                          prob_above_1 = 'Probability above 1') %>%
        add_header_row(values = c('contrast', 'experiment', 'Hazard ratio', 'Probability above 1'), colwidths = c(1,1,3,1)) %>%
        merge_v(j = c(1,2,6), part = 'header') %>%
        font(fontname = 'Times', part = 'all') %>%
        fontsize(size = 14, part = 'all') %>%
        align(align = 'center', part = 'all') %>%
        merge_v(j = ~ contrast) %>%
        autofit() %>%
        hline(i = c(3,5,7,9,11,13,15,17,19,21), border = fp_border_default()) %>%
        hline_top(part = 'header', border = fp_border(width = 1.5)) %>%
        hline_top(part = 'body', border = fp_border(width = 1.5)) %>%
        hline(i = 23, part = 'body', border = fp_border(width = 1.5)) %>%
        fix_border_issues()
        
# get a png from the html file with webshot
save_as_image(table, 'tables/Table_S5.png', zoom = 3, webshot = 'webshot2')
