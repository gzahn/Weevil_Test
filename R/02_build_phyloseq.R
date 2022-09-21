library(tidyverse); packageVersion("tidyverse")
library(dada2); packageVersion("dada2")
library(decontam); packageVersion("decontam")
library(phyloseq); packageVersion("phyloseq")
library(purrr); packageVersion("purrr")
library(Biostrings); packageVersion("Biostrings")
library(janitor); packageVersion("janitor")
library(skimr); packageVersion("skimr")
library(readxl); packageVersion("readxl")


# IMPORT METADATA ####
meta <- readxl::read_xlsx("./data/metadata_amplicon_18s_its2.xlsx") %>% 
  clean_names()

# Clean metadata
meta <- 
meta %>% 
  mutate(sample_id = file_name_amf_18s_forward %>% str_split("_") %>% map_chr(1),
         weevil_treatment = case_when(weevil_nw_no_weevil_w_weevil_present == "NW" ~ "Absent",
                                      weevil_nw_no_weevil_w_weevil_present == "W" ~ "Present"),
         nitrogen = nitrogen_kg_n_ha_yr,
         nitrogen_unit = "kg/ha/yr",
         block = factor(block),
         plot_number = factor(plot_number)) %>% 
  select(site,plot_number,nitrogen,nitrogen_unit,block,sample_id,weevil_treatment)


# IMPORT SEQTABLE ####
seqtab <- readRDS("./output/18S_seq_table.RDS")

# IMPORT TAXONOMY ####
taxa <- readRDS("./output/18S_taxonomy.RDS")
taxa2 <- readRDS("./output/18S_w_silva_taxonomy.RDS")

row.names(meta) <- meta$sample_id


# reorder otu table and metadata to match
row.names(meta) <- meta$sample_id
seqtab <- seqtab[row.names(meta),]

# reorder tax table and otu columns to match
seqtab <- seqtab[,row.names(taxa)]

# check matching orders
identical(colnames(seqtab), row.names(taxa))
identical(row.names(seqtab), row.names(meta))

otu <- otu_table(seqtab,taxa_are_rows = FALSE)
tax <- tax_table(taxa)
tax2 <- tax_table(taxa2)
met <- sample_data(meta)

# fix metadata sample names
sample_names(met) <- meta$sample_id

# combine into phyloseq object ####
ps <- phyloseq(otu,met,tax)
ps_silva <- phyloseq(otu,met,tax2)

# REMOVE CONTAMINANTS ####
# no negative PCR controls to use!
# contams <- decontam::isContaminant(otu_table(ps),neg = ps@sam_data$locationtype == "LAB BLANK",method = "prevalence") # check samdata for appropriate colname!
# ps.noncontam <- prune_taxa(!contams$contaminant, ps)

# REMOVE empty samples/taxa ####
ps <- subset_taxa(ps, taxa_sums(ps) > 0)
ps <- subset_samples(ps, sample_sums(ps) > 0)
ps_silva <- subset_taxa(ps_silva, taxa_sums(ps_silva) > 0)
ps_silva <- subset_samples(ps_silva, sample_sums(ps_silva) > 0)

# SAVE PS OBJECT
saveRDS(ps,"output/partially_cleaned_ps_object.RDS")
saveRDS(ps_silva,"output/partially_cleaned_ps_silva_object.RDS")
