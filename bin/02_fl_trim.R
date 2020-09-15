#!/usr/bin/env Rscript

# step 2
# fixed length trimming

# args are the fastq files


#require(dada2)
require(rexmap)
require(argparse)

parser <- ArgumentParser(description='fixed length trimming')
parser$add_argument('fastqfiles', metavar = "fastqfiles", nargs = '+', help = 'fastq files')
parser$add_argument('--prob', default = 0.01, 
										help = 'The probability to use with the ftquantile() function, which is used to find the minimum seq len above which we have certain fraction of the reads')
# https://docs.python.org/3/library/argparse.html#argparse.ArgumentParser.add_argument

args <- parser$parse_args()

# prepare file names
sample_ids <- rexmap::sampleids_from_filenames(args$fastqfiles, separator = "_")
fq_fltrimmed <- file.path(paste0(sample_ids, '_fltrimmed.fastq'))

# perform fl trim
seqlen.ft <- sequence_length_table(args$fastqfiles)
trimlen <- ftquantile(seqlen.ft, prob = args$prob)
filtstats <- filter_and_trim(args$fastqfiles, fq_fltrimmed, truncLen=trimlen, verbose=T)

