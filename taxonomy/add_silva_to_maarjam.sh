#!/bin/bash
set -ueo pipefail


zcat SILVA_138.1_SSURef_NR99_tax_silva.fasta.gz | seqtk seq > SILVA_138.1_SSURef_NR99_tax_silva_onleline.fasta
cat SILVA_138.1_SSURef_NR99_tax_silva_onleline.fasta | grep -A1 "Eukaryota" | grep -v "^--$" > SILVA_138.1_SSURef_NR99_tax_silva_onleline_euks.fasta
rm SILVA_138.1_SSURef_NR99_tax_silva_onleline.fasta
cat SILVA_138.1_SSURef_NR99_tax_silva_onleline_euks.fasta | sed 's/ Eukaryota;/|k__/' | sed 's/;.*/;p__NA;c__NA;o__NA;f__NA;g__NA;s__NA/' > SILVA_138.1_SSURef_NR99_tax_silva_euks_formatted.fasta
rm SILVA_138.1_SSURef_NR99_tax_silva_onleline_euks.fasta
cat SILVA_138.1_SSURef_NR99_tax_silva_euks_formatted.fasta | sed 's/.*|/>SILVA_SEQ|NA|NA|NA|/' > SILVA_reformatted_headers
cat <(zcat maarjam_database_SSU_reformatted.fasta.gz) SILVA_reformatted_headers > maarjam_and_silva.fasta
rm SILVA_reformatted_headers
rm SILVA_138.1_SSURef_NR99_tax_silva_euks_formatted.fasta
gzip maarjam_and_silva.fasta
