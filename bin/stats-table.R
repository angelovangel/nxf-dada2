#!/usr/bin/env Rscript
# make a stats table with reads going through each step in the pipeline

require(argparse, quietly = TRUE, warn.conflicts = FALSE)
require(dplyr, quietly = TRUE, warn.conflicts = FALSE)
require(data.table, quietly = TRUE, warn.conflicts = FALSE)

parser <- ArgumentParser(description='make a stats table')
parser$add_argument('--stats_merge_file', default = 'stats-merge.csv')
parser$add_argument('--stats_fltrim_file', default = 'stats-fltrim.csv')
args <- parser$parse_args()

mstats <- fread(args$stats_merge_file)
flstats <- fread(args$stats_fltrim_file)
stats <- cbind(mstats, flstats)

stats[,-c(3:5)] %>% 
	mutate(percent_retained = reads.out/total_reads*100) %>%
	fwrite(file = "stats-read-counts.csv")