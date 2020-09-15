#!/usr/bin/env Rscript

# step 2
# fixed length trimming

# args are the fltrimmed and fluntrimmed fastq files
#require(dada2)
require(rexmap)
require(argparse)

parser <- ArgumentParser(description='dada2 denoise')
parser$add_argument('--fltrimmed', metavar = "trimmed fastqfiles", nargs = '+', help = 'fixed len trimmed fastq files')
parser$add_argument('--fluntrimmed', metavar = "untrimmed fastqfiles", nargs = '+', help = 'fixed length untrimmed fastq files')

args <- parser$parse_args()

dada_result <- dada_denoise(args$fltrimmed, args$fluntrimmed, verbose = TRUE, timing = TRUE)
ab.dt <- sequence_abundance(dada_result)
ab.dt.nobim <- sequence_abundance(dada_result, remove_bimeras = FALSE)

saveRDS(dada_result, file = "dada2_denoise.rds")