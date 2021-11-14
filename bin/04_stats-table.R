#!/usr/bin/env Rscript
# make a stats table with reads going through each step in the pipeline

require(argparse, quietly = TRUE, warn.conflicts = FALSE)
require(dplyr, quietly = TRUE, warn.conflicts = FALSE)
require(data.table, quietly = TRUE, warn.conflicts = FALSE)

parser <- ArgumentParser(description='make a stats table')
parser$add_argument('--stats_merge_file', default = 'stats-merge.csv')
parser$add_argument('--stats_fltrim_file', default = 'stats-fltrim.csv')
parser$add_argument('--ab_file', default = 'ab.dt')
parser$add_argument('--osu_ab_file', default = 'osu_abundances.csv')
args <- parser$parse_args()

mstats <- fread(args$stats_merge_file) %>% arrange(sample_id)
flstats <- fread(args$stats_fltrim_file) %>% arrange(file)
abstats <- fread(args$ab_file) %>% group_by(sample_id) %>% summarise_at("raw_count", sum) %>% arrange(sample_id)
osustats <- fread(args$osu_ab_file) %>% group_by(sample_id) %>% summarise_at("osu_count", sum) %>% arrange(sample_id)

if( length(unique(c(nrow(mstats), 
										nrow(flstats), 
										nrow(abstats), 
										nrow(osustats))
									) 
					 ) != 1 ) { 
	stop("something is wrong with the stats files!")
	}
stats <- cbind(mstats, flstats, asv_count = abstats$raw_count, osu_count = osustats$osu_count)

stats %>% 
	mutate(percent_retained = reads.out/total_reads*100) %>%
	dplyr::rename(after_merge = reads.in, after_fixlen_trim = reads.out) %>%
	dplyr::select(sample_id, total_reads, after_merge, after_fixlen_trim, percent_retained, asv_count, osu_count) %>%
	fwrite(file = "stats-read-counts.csv")
