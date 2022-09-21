library(tidyverse); packageVersion("tidyverse")
library(patchwork); packageVersion("patchwork")
library(phyloseq); packageVersion("phyloseq")
library(vegan); packageVersion("vegan")
library(broom); packageVersion("broom")
library(purrr); packageVersion("purrr")
library(vegan); packageVersion("vegan")
library(ggmap); packageVersion("ggmap")
library(modelr); packageVersion("modelr")

# custom functions
source("./R/plot_bar2.R")

# LOAD PS OBJECTS ####
ps <- readRDS("./output/partially_cleaned_ps_object.RDS")
ps_silva <- readRDS("output/partially_cleaned_ps_silva_object.RDS")

# look at kingdom-level assignments
ntaxa(ps)
ntaxa(ps_silva)

# export seqs for blastn analysis
# are they in the same order in both ps objects?
identical(phyloseq::taxa_names(ps),phyloseq::taxa_names(ps_silva))

maarjam_seqs <- phyloseq::taxa_names(ps)
maarjam_names <- paste0("ASV_",seq_along(maarjam_seqs))
maarjam_seqs <- Biostrings::DNAStringSet(maarjam_seqs)
maarjam_names <- Biostrings::BStringSet(maarjam_names)
names(maarjam_seqs) <- maarjam_names
ShortRead::writeFasta(maarjam_seqs,"./test.fasta")

# blastn -query ./test.fasta -db maarjam_and_silva -outfmt 6 -max_target_seqs 1 -out ./taxonomy/BLAST_results.txt
# system2("blastn", 
#         args = c("-query","./test.fasta",
#                  "-db","maarjam_and_silva",
#                  "-outfmt",6,
#                  "-max_target_seqs",1,
#                  "-out","./taxonomy/BLAST_results.txt"))
blast <- read_delim("./taxonomy/BLAST_results.txt",col_names = FALSE)
blast_top_hit_table <- 
  blast %>% 
  select(X1,X2,X3) %>% 
  set_names(c("query","top_hit","perc_match")) %>% 
  mutate(top_hit = top_hit %>% str_split("\\|") %>% map_chr(5) %>% 
           str_split(";") %>% map_chr(1) %>% str_remove("k__"))

blast_top_hit_table$top_hit %>% table()

fungal_blast_hits <- 
blast_top_hit_table %>% 
  filter(top_hit == "Fungi") %>% 
  pluck("query") %>% 
  str_remove("ASV_") %>% 
  as.numeric()

tax_table(ps)[fungal_blast_hits,] %>% row.names

nonfungal_blast_hits <- 
  blast_top_hit_table %>% 
  filter(top_hit != "Fungi") %>% 
  pluck("query") %>% 
  str_remove("ASV_") %>% 
  as.numeric()

tax_table(ps)[nonfungal_blast_hits,] %>% row.names

# look at those from maarjam that were identified as fungi
# blast those seqs against ncbi nr
ps_fungal_seqs <- 
  ps@tax_table[grep("Fungi",ps@tax_table[,1]),] %>% 
  row.names()

ps_fungal_names <- paste0("ASV_",seq_along(ps_fungal_seqs))
fungal_seqs <- Biostrings::DNAStringSet(ps_fungal_seqs)
ps_fungal_names <- Biostrings::BStringSet(ps_fungal_names)
names(ps_fungal_seqs) <- ps_fungal_names
ShortRead::writeFasta(ps_fungal_seqs,"./ps_fungal_test.fasta")
# BLAST against NCBI nt database (on HPC) 

ps_silva_fungal_seqs <- 
  ps@tax_table[grep("Fungi",ps_silva@tax_table[,1]),] %>% 
  row.names()

ps_silva_fungal_names <- paste0("ASV_",seq_along(ps_silva_fungal_seqs))
fungal_seqs <- Biostrings::DNAStringSet(ps_silva_fungal_seqs)
ps_silva_fungal_names <- Biostrings::BStringSet(ps_silva_fungal_names)
names(ps_silva_fungal_seqs) <- ps_silva_fungal_names
ShortRead::writeFasta(ps_silva_fungal_seqs,"./ps_silva_fungal_test.fasta")


# read in BLAST results
ps_fungal_blast <- read_delim("./ps_fungal_BLAST_results.txt",col_names = FALSE)
ps_silva_fungal_blast <- read_delim("./ps_silva_fungal_BLAST_results.txt",col_names = FALSE)
maarjam_blast <- read_delim("./weevil_BLAST_results.txt",col_names = FALSE)

head(maarjam_blast)
head(ps_fungal_blast)

# remove duplicate assignments
ps_fungal_blast <- ps_fungal_blast[match(unique(ps_fungal_blast$X1),ps_fungal_blast$X1),]
ps_silva_fungal_blast <- ps_silva_fungal_blast[match(unique(ps_silva_fungal_blast$X1),ps_silva_fungal_blast$X1),]

# paste them together with RDP results

maarjam_assignments <- 
data.frame(
  maarjam_seq = ps_fungal_seqs,
  RDP_Kingdom = ps@tax_table[ps_fungal_seqs,1],
  ncbi_top_hit = ps_fungal_blast$X5[ps_fungal_blast$X1 %in% names(ps_fungal_seqs)]
)

maarjam_and_silva_assignments <- 
  data.frame(
    maarjam_silva_seq = ps_silva_fungal_seqs,
    RDP_Kingdom = ps_silva@tax_table[ps_silva_fungal_seqs,1],
    ncbi_top_hit = ps_silva_fungal_blast$X5[ps_silva_fungal_blast$X1 %in% names(ps_silva_fungal_seqs)]
  )

write_csv(maarjam_assignments,"./taxonomy/maarjam_vs_ncbi_taxonomy.csv")
write_csv(maarjam_and_silva_assignments,"./taxonomy/maarjam-silva_vs_ncbi_taxonomy.csv")


maarjam_and_silva_assignments %>% View

####









# # remove NA kingdoms
# ps <- ps %>%
#   subset_taxa(!is.na(Kingdom))
# ps_silva <- ps_silva %>%
#   subset_taxa(!is.na(Kingdom))
# 
# # add grouping variable
# newvar <- paste0(ps@sam_data$weevil_nw_no_weevil_w_weevil_present,
#                  "_",ps@sam_data$nitrogen_kg_n_ha_yr)
# ps@sam_data$group <- newvar
# ps_silva@sam_data$group <- newvar
# 
# 
# # DEFINITELY NEED TO REVISIT TAXONOMY ASSIGNMENT
# # MAYBE ADD OUTGROUPS TO MAARJAM DATABASE?
# # TRY BLAST INSTEAD?
# 
# # RELABUND
# ps_ra <- ps %>%
#   transform_sample_counts(function(x){x/sum(x)})
# 
# ps_silva_ra <- ps_silva %>%
#   transform_sample_counts(function(x){x/sum(x)})
# 
# # taxa comparison
# ntaxa(ps_ra)
# ntaxa(ps_silva_ra)
# 
# 
# 
# # fix "Kingdom" names
# ps_ra@tax_table[,1] <- ps_ra@tax_table[,1] %>% str_split("\\|") %>% map_chr(5)
# ps_silva_ra@tax_table[,1] <- ps_silva_ra@tax_table[,1] %>% str_split("\\|") %>% map_chr(5)
# 
# 
# # quick barplots
# plot_bar2(ps_ra,fill="Kingdom")
# plot_bar2(ps_silva_ra,fill="Kingdom")
# 
# p1 <- plot_bar2(ps_ra,fill="Family")
# p2 <- plot_bar2(ps_silva_ra,fill="Family")
# 
# # alpha plots
# p1 <- ps %>%
#   plot_richness(x="group",measures = "Shannon") +
#   coord_cartesian(ylim=c(0,2.5)) +
#   ggtitle("maarjam database") +
#   facet_wrap(~weevil_nw_no_weevil_w_weevil_present,scales = "free_x")
# p2 <- ps_silva %>%
#   plot_richness(x="group",measures = "Shannon") +
#   coord_cartesian(ylim=c(0,2.5)) +
#   ggtitle("maarjam + SILVA outgroups") +
#   facet_wrap(~weevil_nw_no_weevil_w_weevil_present,scales = "free_x")
# 
# p1 + p2
# paste(ps_silva@tax_table[,5],
#       ps_silva@tax_table[,6],
#       ps_silva@tax_table[,7])
# 
# plot_bar2(ps_silva_ra,fill="Family",x="group")
# # merge samples by treatment group
# psm <- merge_samples(ps,"group") %>%
#   transform_sample_counts(function(x){x/sum(x)})
# psm_silva <- merge_samples(ps_silva,"group") %>%
#   transform_sample_counts(function(x){x/sum(x)})
# 
# phyloseq::plot_richness(ps,measures = "Shannon")
