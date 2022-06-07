#----------------------------------------------------------------------#
# prevalence filter and do initial screening of 16s sequencing data ####
#----------------------------------------------------------------------#

# what this script does:
# does some quality control on the phyloseq object from the dada2 16S pipeline
# 1. removes any ASV not assigned to the level of phylum
# 2. remove ASVs with <10% prevalence or < 200 reads in total
# 3. remove ASVs with a length > 250bp
# 4. remove ASVs assigned to cyanobacteria or mitochondria
# 5. remove samples with total read numbers <8000
# 6. rarefy and assign reproducible tree root

# load packages
library(patchwork)
library(phyloseq)
library(tidyverse)

# set seed
set.seed(42)

# load data - latest run which we are happy with ####
ps <- readRDS('data/phyloseq_16s.rds')

# look at object
ps

sort(sample_sums(ps))

# first...
rank_names(ps)

#-----------------------------------------------------------#
# 1. removes any ASV not assigned to the level of phylum ####
#-----------------------------------------------------------#

# look number of features for each phyla
table(tax_table(ps)[, "Phylum"], exclude = NULL)
# lots of NAs here. These are probably artifacts and can be removed

# remove NA characterisation at the Phylum level
ps0 <- subset_taxa(ps, !is.na(Phylum))
table(tax_table(ps0)[, "Phylum"], exclude = NULL)

#----------------------------------------------------------------#
# 2. remove ASVs with <10% prevalence or < 200 reads in total ####
#----------------------------------------------------------------#

# explore prevalence of the dataset
# Compute prevalence of each feature, store as data.frame
prevdf = apply(X = otu_table(ps0),
               MARGIN = ifelse(taxa_are_rows(ps0), yes = 1, no = 2),
               FUN = function(x){sum(x > 0)})

# Add taxonomy and total read counts to this data.frame
prevdf = data.frame(prevalence = prevdf,
                    total_abundance = taxa_sums(ps0),
                    tax_table(ps0), stringsAsFactors = FALSE) %>%
  mutate(., feature = row.names(.))

# head(prevdf)

# Are there phyla that are comprised of mostly low-prevalence features? Compute the total and average prevalences of the features in each phylum.
prev_sum <- group_by(prevdf, Phylum) %>%
  summarise(., mean_prevalence = mean(prevalence)/nsamples(ps0)*100,
            total_prevalence = sum(prevalence)) %>%
  ungroup()

prev_sum

# get unique taxa for prevalence and plot
prevdf1 <- filter(prevdf, Phylum %in% get_taxa_unique(ps0, "Phylum"))

# plot prevalence of samples by 
ggplot(prevdf1, aes(total_abundance, prevalence / nsamples(ps0),color=Phylum)) +
  # Include a guess for parameter
  geom_hline(yintercept = 0.05, alpha = 0.1, linetype = 2) + geom_point(size = 2, alpha = 0.7) +
  scale_x_log10() +  xlab("Total Abundance") + ylab("Prevalence [Frac. Samples]") +
  facet_wrap(~ Phylum) + theme(legend.position="none")

# delete things under 10% prevalence
# min total read abundance = 200
ps_prevfilt <- microViz::tax_filter(ps0, min_prevalence = 0.1,
                          min_total_abundance = 200)

# fix ranks - remove Species and only keep Species.1 - as per Ben Callahan's recommendation
tax_table(ps_prevfilt) <- tax_table(ps_prevfilt)[,c(1:6, 8)]
colnames(tax_table(ps_prevfilt)) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

#-----------------------------------------#
# 3. remove ASVs with a length > 250bp ####
#-----------------------------------------#

d_ps <- speedyseq::psmelt(ps_prevfilt) %>%
  janitor::clean_names() %>%
  group_by(sample) %>%
  mutate(rel_abundance = abundance/sum(abundance)) %>%
  ungroup()

# look at distribution of OTU length
ggplot(d_ps, aes(nchar(otu))) +
  geom_histogram()

summary(nchar(d_ps$otu))

# remova OTUs longer than 250 bps
all_taxa <- taxa_names(ps)
to_drop <- all_taxa[nchar(all_taxa) > 250]
to_keep <- all_taxa[!all_taxa %in% to_drop]

ps <- prune_taxa(to_keep, ps_prevfilt)

d_ps <- speedyseq::psmelt(ps_prevfilt) %>%
  janitor::clean_names() %>%
  group_by(sample) %>%
  mutate(rel_abundance = abundance/sum(abundance)) %>%
  ungroup()

# look at distribution of OTU length
ggplot(d_ps, aes(nchar(otu))) +
  geom_histogram()

#-------------------------------------------------------------#
# 4. remove ASVs assigned to cyanobacteria or mitochondria ####
#-------------------------------------------------------------#

# remove cyanobacteria and mitochondria
ps <- subset_taxa(ps, Phylum != 'Cyanobacteria')
ps <- subset_taxa(ps, Family != 'Mitochondria')

#----------------------------------------------------#
# 5. remove samples with total read numbers <8000 ####
#----------------------------------------------------#

# check sample numbers
sample_sums(ps) %>% sort()

to_drop <- tibble(reads = sample_sums(ps), sample = names(sample_sums(ps))) %>%
  filter(reads < 8000)
to_keep <- tibble(reads = sample_sums(ps), sample = names(sample_sums(ps))) %>%
  filter(reads > 8000)

# remove these samples from ps and create a rarefaction
ps2 <- prune_samples(to_keep$sample, ps)
sort(sample_sums(ps2))

#------------------------------------------------#
# 6. rarefy and assign reproducible tree root ####
#------------------------------------------------#

# rarefy these samples down to the lowest number of reads 
ps2_rarefy <- rarefy_even_depth(ps2)
# save out datasets

# assign reproducible tree root to the dataset
# https://john-quensen.com/r/unifrac-and-tree-roots/
# Sebastian Schmidt proposed choosing an out group based on the longest tree branch terminating in a tip
# https://github.com/joey711/phyloseq/issues/597

# function for picking tree root:
pick_new_outgroup <- function(tree.unrooted){
  # tablify parts of tree that we need.
  treeDT <-
    cbind(
      data.table::data.table(tree.unrooted$edge),
      data.table::data.table(length = tree.unrooted$edge.length)
    )[1:ape::Ntip(tree.unrooted)] %>%
    cbind(data.table::data.table(id = tree.unrooted$tip.label))
  # Take the longest terminal branch as outgroup
  new.outgroup <- treeDT[which.max(length)]$id
  return(new.outgroup)}

# set root for rarefied sample then the same one for the none-rarefied sample
out_group <- phy_tree(ps2_rarefy) %>% pick_new_outgroup()

phy_tree(ps2_rarefy) <- ape::root(phy_tree(ps2_rarefy), outgroup=out_group, resolve.root=TRUE)
phy_tree(ps2) <- ape::root(phy_tree(ps2), outgroup=out_group, resolve.root=TRUE)

saveRDS(ps2, 'data/phyloseq_16s_nocyano_nomitoo.rds')
saveRDS(ps2_rarefy, 'data/phyloseq_16s_rarefied_nocyano_nomito.rds')


