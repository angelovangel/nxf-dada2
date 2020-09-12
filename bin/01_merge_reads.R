#!/usr/bin/env Rscript

# step 1
# merge reads

# args are the fastq files

require(dada2)
require(rexmap)

args <- commandArgs(trailingOnly = T)

# 
fq_fwd <- args[grepl("R1", args)]
fq_rev <- args[grepl("R2", args)]

sample_ids <- rexmap::sampleids_from_filenames(fq_fwd, separator = "_") ###################TODO check if unique
fq_merged <- file.path(paste0(sample_ids, '_merged.fastq'))

# actual merge
mergestats <- rexmap::merge_pairs(fq_fwd, fq_rev, fq_merged, verbose=TRUE, timing = TRUE)

mergestatsdf <- data.frame(
	sample_ids, 
	total_reads = as.numeric(mergestats[1, ]), 
	low_pct_sim = as.numeric(mergestats[2, ]), 
	low_aln_len = as.numeric(mergestats[3, ]))

write.csv(mergestatsdf, file = "mergestats.csv")
