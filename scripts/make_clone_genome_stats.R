# make genome summary table

library(tidyverse)
library(flextable)
library(officer)

d_meta <- readxl::read_excel('data/sequencing/clone/clone_metadata.xlsx') %>%
  mutate(sample = paste('sample', CODE, sep = ''),
         site = substr(SAMPLE, 1, 1),
         rep = substr(SAMPLE,2,2),
         plastic = substr(SAMPLE,3,3),
         plastic = case_when(plastic == 'T' ~ 'PET',
                             plastic =='V' ~ 'PVC',
                             plastic == 'P' ~ 'PP',
                             plastic == 'L' ~ 'LDPE',
                             plastic == 'H' ~ 'HDPE')) %>%
  select(sample, site, rep, plastic)

d_mmseq <- read.csv('data/sequencing/clone/unicycler/assignment_mmseq.csv') %>%
  mutate(file = gsub('.fasta', '', file),
         sample = gsub('\\..*', '', file))
d_16s <- read.csv('data/sequencing/clone/unicycler/assignment_16s.csv') %>%
  mutate(file = gsub('.fasta', '', file),
         sample = gsub('\\..*', '', file))
d_checkm <- read.csv('data/sequencing/clone/unicycler/checkm.csv') %>%
  mutate(sample = gsub('\\..*', '', file))
d_amr <- read.csv('data/sequencing/clone/unicycler/amrfinder/sample_amr_genes.csv') %>%
  mutate(sample = gsub('\\..*', '', file))

d_taxa <- left_join(d_meta, d_mmseq) %>%
  left_join(., d_16s) %>%
  left_join(., d_checkm) %>%
  left_join(., d_amr)

# remove doubles
d_taxa <- filter(d_taxa, ! file %in% c('sample11', 'sample31', 'sample20', 'sample6')) %>%
  mutate(assignment_16s = replace_na(assignment_16s, '-')) %>%
  select(-file)

write.csv(d_taxa, 'data/sequencing/clone/unicycler/clone_genome_stats.csv')

d_taxa <- mutate(d_taxa, site = case_when(site == 'A' ~ 'dockyard',
                                          site == 'B' ~ 'falmouth wharf',
                                          site == 'C' ~ 'mylor bank'),
                 across(where(is.character),trimws))

d_table <- dplyr::select(d_taxa, sample, plastic, site, assignment_mmseq, assignment_16s, completeness, contamination, GC, genome_size, n_contigs, N50_contig, number_amr) %>%
  arrange(plastic, sample) %>%
  mutate(across(completeness:GC, round, 2),
         isolate = 1:n(),
         assignment_mmseq = ifelse(isolate == 16, 'Klebsiella/Serratia\nEnterobacter/Enterobacteriaceae', assignment_mmseq)) %>%
  select(isolate, everything())
d_table

table1 <- flextable(select(d_table, -sample)) %>%
  set_header_labels(assignment_mmseq = 'contig taxonomic\nassignment',
                    assignment_16s = '16S taxonomic\nassignment',
                    genome_size = 'genome size (bp)',
                    n_contigs = 'number\nof contigs',
                    N50_contig = 'N50',
                    number_amr = 'number of\nAMR genes'
                    ) %>%
  font(fontname = 'Times', part = 'all') %>%
  fontsize(size = 10, part = 'all') %>%
  align(align = 'center', part = 'all') %>%
  hline(i = c(1,3,4,5,7,9,10,11,12,13,14,16,17,18,19), border = fp_border_default()) %>%
  fix_border_issues() %>%
  italic(j = c('assignment_mmseq', 'assignment_16s')) %>%
  autofit() %>%
  width(j = 'site', 4) %>%
  width(j = 'assignment_mmseq', 1) %>%
  bg(bg = 'white', part = 'all')
table1
save_as_image(table1, 'tables/genome_stats.png', webshot = 'webshot')
save_as_docx(table1, path = 'tables/table_1.docx')

d_table2 <- select(d_taxa, plastic, sample, assignment_16s, amr_genes, number_amr) %>%
  arrange(plastic, sample) %>%
  mutate(isolate = 1:n()) %>%
  select(isolate, amr_genes) %>%
  separate(amr_genes, paste('random', 1:10), sep = ',') %>%
  pivot_longer(starts_with('random'), names_to = 'stuff', values_to = 'gene', values_drop_na = TRUE)

rownames_to_column(d_table2) %>% group_by(isolate) %>% summarise(max(rowname))

d_table3 <- filter(d_table2, isolate <= 13)
d_table4 <- filter(d_table2, isolate > 13) %>% mutate(isolate = as.character(isolate))
#d_table5 <- tibble(isolate = rep(max(d_table4$isolate), times = 6), stuff = NA, gene = '') 
d_table5 <- tibble(isolate = rep('', times = 6), stuff = NA, gene = '') 
d_table4 = bind_rows(d_table4, d_table5)

d_table2 <- cbind(select(d_table3, -stuff), select(d_table4, isolate2 = isolate, gene2 = gene))

table2 <- flextable(d_table2,
                    col_keys = c('isolate', 'gene', 'blank1', 'isolate2', 'gene2')) %>%
  empty_blanks() %>%
  set_header_labels(gene = 'AMR gene',
                    gene2 = 'AMR gene',
                    isolate2 = 'isolate') %>%
  font(fontname = 'Times', part = 'all') %>%
  fontsize(size = 10, part = 'all') %>%
  align(align = 'center', part = 'header') %>%
  align(align = 'left', part = 'body') %>%
  align(align = 'center', part = 'body', j = 1) %>%
  hline(i = c(2,8,9,15,17,18,21,22,25,31,32,34), j = c(1,2), border = fp_border_default()) %>%
  hline(i = c(6,7,17,24,30,32, 34), j = c(4,5), border = fp_border_default()) %>%
  merge_v(., ~isolate) %>%
  merge_v(., ~isolate2) %>%
  valign(valign = 'top', j = 1, part = 'body') %>%
  valign(valign = 'top', j = 4, part = 'body') %>%
  border(j = c(4,5), i = 40, border.bottom = fp_border(color = "white")) %>%
  fix_border_issues() %>%
  autofit()
table2
save_as_image(table2, 'tables/amr_genes.png', zoom = 3, webshot = 'webshot2')
save_as_docx(table2, path = 'tables/table_2.docx')

