#!/usr/bin/env Rscript

# step 1
# merge reads and trim primers from merged reads

# args are the fastq files
# artefacts are 

#require(dada2)
require(rexmap)
require(argparse)

parser <- ArgumentParser(description='merge and primer-trim fastq files')
parser$add_argument('fastqfiles', metavar = "fastqfiles", nargs = '+', help = 'fastq files')
parser$add_argument('--region', default = 'V3-V4', 
	help = 'Hypervariable region used. Used to automatically retrieve PCR primer sequences for trimming. Can be one of V1-V3, V3-V4 or V4')
# https://docs.python.org/3/library/argparse.html#argparse.ArgumentParser.add_argument

args <- parser$parse_args()

# 
fq_fwd <- args$fastqfiles[grepl("R1", args$fastqfiles)]
fq_rev <- args$fastqfiles[grepl("R2", args$fastqfiles)]
# 
sample_ids <- rexmap::sampleids_from_filenames(fq_fwd, separator = "_") ###################TODO check if unique
fq_merged <- file.path(paste0(sample_ids, '_merged.fastq'))
fq_pcrtrimmed <- file.path(paste0(sample_ids, '_trimmed.fastq'))
# 
# actual merge and pcr-trim
mergestats <- rexmap::merge_pairs(fq_fwd, fq_rev, fq_merged, verbose=TRUE, timing = TRUE)
trimstats  <- rexmap::remove_pcr_primers(fq_merged, fq_pcrtrimmed, region = args$region, verbose = TRUE, timing = TRUE)
# 
mergestatsdf <- data.frame(
	sample_id = sample_ids,
	total_reads = as.numeric(mergestats[1, ]),
	low_pct_sim = as.numeric(mergestats[2, ]),
	low_aln_len = as.numeric(mergestats[3, ])
)
# 
trimstatsdf <- data.frame(
	sample_id = sample_ids,
	total = as.numeric(mergestats[1, ]),
	fwd_trim = as.numeric(trimstats[1,]),
	rev_trim = as.numeric(trimstats[2,])
	)
# 
# 
write.table(mergestatsdf, file = "merge.csv", row.names = FALSE, col.names = FALSE, sep = ",")
write.table(trimstatsdf, file = "primertrim.csv", row.names = FALSE, col.names = FALSE, sep = ",")
