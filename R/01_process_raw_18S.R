
# PACKAGES, SCRIPTS, AND SETUP ####
library(tidyverse); packageVersion("tidyverse")
library(dada2); packageVersion("dada2")
library(decontam); packageVersion("decontam")
library(phyloseq); packageVersion("phyloseq")
library(purrr); packageVersion("purrr")
library(Biostrings); packageVersion("Biostrings")

#################################################################################
#                               Main workflow                                   #
# Filter and trim, denoise, sample inferrence, chimera and contaminant removal, # 
# taxonomic assignment, combine sequence table and metadata                     #
#################################################################################

# PARSE FILE PATHS ####

# File parsing - For this, we will use only the forward illumina reads - make sure to move fwd reads into their own directory for simplest processing
path <- "./data/filtN" # CHANGE to the directory containing your demultiplexed fastq files
filtpath <- file.path(path, "filtered") # Filtered files go into the filtered/ subdirectory
if(!file_test("-d", filtpath)) dir.create(filtpath) # make directory for filtered fqs if not already present
fnFs <- sort(list.files(file.path(path), full.names = TRUE, pattern = "R1_18s.fastq.gz"))
fnRs <- sort(list.files(file.path(path), full.names = TRUE, pattern = "R2_18s.fastq.gz"))


samplenames <- unlist(map(strsplit(basename(fnFs), "_"), 1))



# visualize a couple of fwd read quality profiles to help select reasonable filtration parameters
plotQualityProfile(fnFs[1:2])
plotQualityProfile(fnRs[1:2])

# FILTER AND TRIM ####
filtFs <- file.path(filtpath, paste0(samplenames, "_F_filt.fastq.gz"))
filtRs <- file.path(filtpath, paste0(samplenames, "_R_filt.fastq.gz"))

out <- filterAndTrim(fnFs, filtFs,
                     fnRs, filtRs,
                     maxN=0, maxEE=2, truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE)
head(out)

# sanity check  comparison of before and after filtration
plotQualityProfile(c(fnFs[1:2],filtFs[1:2]))

# LEARN ERROR RATES ####
# Since some samples may have had zero reads pass QC, reassign filts and samplenames
filtsF <- sort(list.files(filtpath, full.names = TRUE, pattern = "_F_"))
filtsR <- sort(list.files(filtpath, full.names = TRUE, pattern = "_R_"))
samplenames <- unlist(map(strsplit(basename(filtsF), "_"),1))

errF <- learnErrors(filtsF, multithread=TRUE, MAX_CONSIST = 20)
errR <- learnErrors(filtsR, multithread=TRUE, MAX_CONSIST = 20)

# sanity check for error model
plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)

# DEREPLICATION ####
derepFs <- derepFastq(filtsF, verbose=TRUE)
derepRs <- derepFastq(filtsR, verbose=TRUE)


# Name the derep-class objects by the sample names
names(derepFs) <- samplenames
names(derepRs) <- samplenames

# SAMPLE INFERRENCE ####
dadaFs <- dada(derepFs, err=errF, multithread=TRUE, selfConsist = TRUE, verbose=TRUE, pool = "pseudo")
dadaRs <- dada(derepRs, err=errR, multithread=TRUE, selfConsist = TRUE, verbose=TRUE, pool = "pseudo")

# MERGE PAIRED READS ####
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)

# MAKE SEQUENCE TABLE ####
seqtab <- makeSequenceTable(mergers)

# REMOVE CHIMERAS ####
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab)

# reassign "out" to remove any missing reads
out = out[as.data.frame(out)$reads.out > 0,]

# Remove all seqs with fewer than 100 nucleotides ####
keeper_esvs <- nchar(names(as.data.frame(seqtab.nochim))) > 99
seqtab.nochim <- seqtab.nochim[,keeper_esvs]

# ASSIGN TAXONOMY ####
taxa <- assignTaxonomy(seqtab.nochim, "./taxonomy/maarjam_database_SSU_reformatted.fasta.gz", multithread=20)

# also assign taxonomy using Maarjam + Silva database
taxa2 <- assignTaxonomy(seqtab.nochim, "./taxonomy/maarjam_and_silva.fasta.gz", multithread=20)
beepr::beep(sound=8)


# Save intermediate files
saveRDS(seqtab.nochim, file = "./output/18S_seq_table.RDS")
saveRDS(taxa, file = "./output/18S_taxonomy.RDS")
saveRDS(taxa2, file = "./output/18S_w_silva_taxonomy.RDS")

