# PACKAGES, SCRIPTS, AND SETUP ####
library(tidyverse); packageVersion("tidyverse")
library(dada2); packageVersion("dada2")
library(decontam); packageVersion("decontam")
library(phyloseq); packageVersion("phyloseq")
library(purrr); packageVersion("purrr")
library(Biostrings); packageVersion("Biostrings")
library(ShortRead); packageVersion("ShortRead")

#################################################################################
#                               Main workflow                                   #
# Filter and trim, denoise, sample inferrence, chimera and contaminant removal, # 
# taxonomic assignment, combine sequence table and metadata                     #
#################################################################################

# PARSE FILE PATHS ####

# File parsing - For this, we will use only the forward illumina reads - make sure to move fwd reads into their own directory for simplest processing
path <- "./data" # CHANGE to the directory containing your demultiplexed fastq files
filtpath <- file.path(path, "filtN") # Filtered files go into the filtered/ subdirectory
if(!file_test("-d", filtpath)) dir.create(filtpath) # make directory for filtered fqs if not already present
list.files(path)

fnFs <- sort(list.files(file.path(path), full.names = TRUE, pattern = "R1_18s.fastq.gz"))
fnRs <- sort(list.files(file.path(path), full.names = TRUE, pattern = "R2_18s.fastq.gz"))
fnFs; fnRs

samplenames <- unlist(map(strsplit(basename(fnFs), "_"), 1))


# CHECK FOR AND REMOVE PRIMER SITES WITH CUTADAPT ####
FWD <- "GCCTCCCTCGCGCCATCAG" # WANDA
REV <- "GAACCCAAACACTTTGGTTTCC"  # AML2

allOrients <- function(primer) {
  # Create all orientations of the input sequence
  require(Biostrings)
  dna <- DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
  orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna), 
               RevComp = reverseComplement(dna))
  return(sapply(orients, toString))  # Convert back to character vector
}
FWD.orients <- allOrients(FWD)
REV.orients <- allOrients(REV)
FWD.orients

# Prefilter to remove reads with ambiguous (N) bases
fnFs.filtN <- file.path(path, "filtN", basename(fnFs)) # Put N-filterd files in filtN/ subdirectory
fnRs.filtN <- file.path(path, "filtN", basename(fnRs)) # Put N-filterd files in filtN/ subdirectory

filterAndTrim(fnFs, fnFs.filtN, 
              fnRs, fnRs.filtN,
              maxN = 0, multithread = TRUE,verbose = TRUE)

# Discover primer matches, regardless of orientation
primerHits <- function(primer, fn) {
  # Counts number of reads in which the primer is found
  nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
  return(sum(nhits > 0))
}
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.filtN[[1]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.filtN[[1]])) # no adaptors in 18S reads

