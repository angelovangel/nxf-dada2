// nextflow-rexmap pipelen
/*
NXF ver 19.08+ needed because of the use of tuple instead of set
*/
if( !nextflow.version.matches('>=19.08') ) {
    println "This workflow requires Nextflow version 19.08 or greater and you are running version $nextflow.version"
    exit 1
}

/* 
 * pipeline input parameters 
 */
params.fastqdir = "fastq"
params.outdir = "$workflow.launchDir/results-dada2"
params.fqpattern = "*_R{1,2}_001.fastq.gz"
params.region = "V3-V4"
params.keep_fastq = false

params.help = ""

/* 
 * pipeline input parameters end
 */

//just in case trailing slash in readsdir not provided...
fastqdir_repaired = "${params.fastqdir}".replaceFirst(/$/, "/") 
fastqfiles = fastqdir_repaired + params.fqpattern

// channel for fastq files --> merge and trim
Channel
    .fromFilePairs( fastqfiles, checkIfExists: true )
    .set { fastq_ch }

// step 1 --> merge reads and trim primers
process merge_and_trim {
    tag "$sample_id"
    publishDir "${params.outdir}/merged_trimmed_reads", enabled: params.keep_fastq, mode: 'copy', pattern: '*_{trimmed,merged}.fastq'
    
    input:
        tuple sample_id, file(x) from fastq_ch
    output:
        file '*_merged.fastq'
        file '*_trimmed.fastq' into mergetrim_ch1
        file '*_trimmed.fastq' into mergetrim_ch2 // 1 goes to fixed_len_trim, 2 goes to dada2
        file 'merge.csv' into mergestats_ch
        file 'primertrim.csv' into trimstats_ch
    script:
    """
    01_merge_and_trim.R --region ${params.region} ${x}
    """
}

// fixed len trimming is done on all files, collected from mergetrim_ch1, therefore a separate process
process fixed_len_trim {
    publishDir "${params.outdir}/merged_trimmed_reads", enabled: params.keep_fastq, mode: 'copy', pattern: '*_fltrimmed.fastq'
    publishDir params.outdir, mode: 'copy', pattern: '*.{pdf,csv}'

    input:
        file x from mergetrim_ch1.collect()
    output:
        file '*_fltrimmed.fastq' into fltrim_ch
        file 'seqlen-fltrim-plot.pdf'
        file 'stats-fltrim.csv'
    script:
    """
    02_fl_trim.R ${x}
    """
}

process merge_and_trim_stats {
    publishDir params.outdir, mode: 'copy'

    input:
        file 'merge*' from mergestats_ch.collect()
        file 'primertrim*' from trimstats_ch.collect()
    output:
        file 'stats*.csv'
    script:
    """
    echo "sample_id, total_reads, low_pct_sim, low_aln_len" > stats-merge.csv
    cat merge* >> stats-merge.csv
    echo "sample_id, total, fwd_trim, rev_trim" > stats-primertrim.csv
    cat primertrim* >> stats-primertrim.csv
    """
}

process dada2_asv {
    publishDir params.outdir, mode: 'copy'

    input:
        file fltrimmed from fltrim_ch
        file fluntrimmed from mergetrim_ch2.collect() // toSortedList?
    output:
        file 'dada2_result.rds'
        file '*.csv' into dada2_ch
    script:
    """
    03_dada.R --fltrimmed ${fltrimmed} --fluntrimmed ${fluntrimmed} --region ${params.region}
    """
}

process phyloseq {
    publishDir params.outdir, mode: 'copy'

    input: 
        file 'osu_*.csv' from dada2_ch.collect()
    output:
        file 'physeq.Rdata'
    script:
    """
    04_phyloseq.R --osu_ab osu_1.csv --osu_tax osu_3.csv
    """
}