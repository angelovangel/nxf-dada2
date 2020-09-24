#!/usr/bin/env Rscript

# rexmap step

# args are the fltrimmed and fluntrimmed fastq files, and region
#require(dada2)
require(rexmap)
require(argparse)

parser <- ArgumentParser(description='dada2 denoise')
parser$add_argument('--fltrimmed', metavar = "trimmed fastqfiles", nargs = '+', help = 'fixed len trimmed fastq files')
parser$add_argument('--fluntrimmed', metavar = "untrimmed fastqfiles", nargs = '+', help = 'fixed length untrimmed fastq files')
parser$add_argument('--region', default = 'V3-V4', help = 'Hypervariable region used. Used in blast')

args <- parser$parse_args()
### sort the arguments here to make sure the fltrimmed and fluntrimmed files match!
args$fltrimmed <- sort(args$fltrimmed)
args$fluntrimmed <- sort(args$fluntrimmed)
###
#cat(args$fltrimmed)
#cat(args$fluntrimmed)


dada_result <- dada_denoise(args$fltrimmed, args$fluntrimmed, verbose = TRUE, timing = TRUE)
ab.dt <- sequence_abundance(dada_result, fq_prefix_split = "_")
ab.dt.nobim <- sequence_abundance(dada_result, remove_bimeras = FALSE)

bimeras_by_unique <- (max(ab.dt.nobim$qseqid) - max(ab.dt$qseqid)) / max(ab.dt.nobim$qseqid) * 100
bimeras_by_abundance <- (ab.dt.nobim[,sum(raw_count)]-ab.dt[,sum(raw_count)])/ab.dt[,sum(raw_count)] * 100

# blast, osu and taxonomy
blast_output <- blast(ab.dt, region = args$region, verbose = TRUE)

osu_ab.dt <- abundance(abundance_table = ab.dt, blast_object = blast_output, verbose = TRUE)
osu_tax.dt <- taxonomy(osu_ab.dt, verbose = TRUE)
osu_seq.dt <- osu_sequences(osu_ab.dt, blast_output)

# save artefacts
# saveRDS(dada_result, file = "dada2_result.rds")
savetable <- function(obj, file) {
	write.table(obj, file = file, sep = "\t", row.names = FALSE)
}
objs <- list(ab.dt, osu_seq.dt, osu_ab.dt, osu_tax.dt)
files <- c("ab.dt", "osu_sequences.csv", "osu_abundances.csv", "osu_taxonomy.csv")
mapply(savetable, objs, files)


