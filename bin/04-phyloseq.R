#!/usr/bin/env Rscript

# go to phyloseq

# args are osu table, tax table, sample table and (optionally) a phylogenetic tree
# return phyloseq object
require(phyloseq)
require(tidyr)
require(dplyr)
require(stringr)
#require(rexmap)
require(argparse)

parser <- ArgumentParser(description='phyloseq')
parser$add_argument('--osu_ab', metavar = "path", help = 'osu abindance table from rexmap')
parser$add_argument('--osu_tax', metavar = "path", help = 'osu taxonomy table from rexmap')
parser$add_argument('--sample_table', metavar = "path", help = 'sample table, the sample_id must match those in the osu table')


args <- parser$parse_args()
osu_ab <- data.table::fread(args$osu_ab)
osu_tax <- data.table::fread(args$osu_tax)

## prepare data for phyloseq ##
# make osu abundance matrix
make_osu_mat <- function(x) {
	osu <- x %>% pivot_wider(names_from = sample_id, values_from = osu_count, values_fill = 0)
	osu_mat <- osu[ ,-c(1:3)] %>% as.matrix()
	rownames(osu_mat) <- osu$osu_id
	osu_mat
}
# make osu tax matrix
make_tax_mat <- function(x) {
	tax_mat <- x[, -c(1:3)] %>% as.matrix()
	rownames(tax_mat) <- x$osu_id
	tax_mat
}

# massage tax table to make it dada2-like 
species <- osu_ab[,c(2,5)] %>% distinct(osu_id, .keep_all = T)
species$species <- str_remove(species$species, ",.+") # remove everything after the comma (i.e. keep first species only)

osu_tax <- osu_tax %>% 
	mutate_if(is.character, stringr::str_replace, "_.+", "") %>% 
	left_join(species)

OSU <- otu_table( make_osu_mat(osu_ab), taxa_are_rows = TRUE )
TAX <- tax_table( make_tax_mat(osu_tax) )

physeq <- phyloseq(
	otu_table(OSU),
	tax_table(TAX)
	)
physeq

