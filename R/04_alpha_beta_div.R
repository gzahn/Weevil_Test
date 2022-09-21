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


# remove NA kingdoms
ps <- ps %>%
  subset_taxa(!is.na(Kingdom))
ps_silva <- ps_silva %>%
  subset_taxa(!is.na(Kingdom))

# add grouping variable
newvar <- paste0(ps@sam_data$weevil_nw_no_weevil_w_weevil_present,
                 "_",ps@sam_data$nitrogen_kg_n_ha_yr)
ps@sam_data$group <- newvar
ps_silva@sam_data$group <- newvar


# DEFINITELY NEED TO REVISIT TAXONOMY ASSIGNMENT
# MAYBE ADD OUTGROUPS TO MAARJAM DATABASE?
# TRY BLAST INSTEAD?

# RELABUND
ps_ra <- ps %>%
  transform_sample_counts(function(x){x/sum(x)})

ps_silva_ra <- ps_silva %>%
  transform_sample_counts(function(x){x/sum(x)})

# taxa comparison
ntaxa(ps_ra)
ntaxa(ps_silva_ra)



# fix "Kingdom" names
ps_ra@tax_table[,1] <- ps_ra@tax_table[,1] %>% str_split("\\|") %>% map_chr(5)
ps_silva_ra@tax_table[,1] <- ps_silva_ra@tax_table[,1] %>% str_split("\\|") %>% map_chr(5)


# quick barplots
plot_bar2(ps_ra,fill="Kingdom")
plot_bar2(ps_silva_ra,fill="Kingdom")

p1 <- plot_bar2(ps_ra,fill="Family")
p2 <- plot_bar2(ps_silva_ra,fill="Family")

# alpha plots
p1 <- ps %>%
  plot_richness(x="group",measures = "Shannon") +
  coord_cartesian(ylim=c(0,2.5)) +
  ggtitle("maarjam database") +
  facet_wrap(~weevil_nw_no_weevil_w_weevil_present,scales = "free_x")
p2 <- ps_silva %>%
  plot_richness(x="group",measures = "Shannon") +
  coord_cartesian(ylim=c(0,2.5)) +
  ggtitle("maarjam + SILVA outgroups") +
  facet_wrap(~weevil_nw_no_weevil_w_weevil_present,scales = "free_x")

p1 + p2
paste(ps_silva@tax_table[,5],
      ps_silva@tax_table[,6],
      ps_silva@tax_table[,7])

plot_bar2(ps_silva_ra,fill="Family",x="group")
# merge samples by treatment group
psm <- merge_samples(ps,"group") %>%
  transform_sample_counts(function(x){x/sum(x)})
psm_silva <- merge_samples(ps_silva,"group") %>%
  transform_sample_counts(function(x){x/sum(x)})

phyloseq::plot_richness(ps,measures = "Shannon")
