#!/usr/bin/env Rscript

# step 1
# merge reads

# args 1 is a folder with fastq file pairs from a 16S amplicon dataset
# args 2 is separator used delimiting sample id from the rest of FASTQ filename

require(dada2)
require(rexmap)

args <- commandArgs(trailingOnly = T)

# read fastq files and get sample ids
fq_fwd <- read_files(args[1], 'R1')
fq_rev <- read_files(args[1], 'R2')

sample_ids <- rexmap::sampleids_from_filenames(fq_fwd, separator = args[2]) ###################TODO check if unique
fq_merged <- file.path(paste0(sample_ids, '_merged.fastq'))

# actual merge
mergestats <- rexmap::merge_pairs(fq_fwd, fq_rev, fq_mer, verbose=TRUE, timing = TRUE)

mergestatsdf <- data.frame(
	sample_ids, 
	total = as.numeric(mergestats[1, ]), 
	low_pct_sim = as.numeric(mergestats[2, ]), 
	low_aln_len = as.numeric(mergestats[3, ]))

write.csv(mergestatsdf, file = "mergestats.csv")
