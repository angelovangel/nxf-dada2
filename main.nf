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
params.outdir = "$workflow.launchDir/results-rexmap"
params.fqpattern = "*_R{1,2}_001.fastq.gz"

params.help = ""

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
    publishDir "${params.outdir}/merged_trimmed_reads", mode: 'copy', pattern: '*_{trimmed,merged}.fastq'
    
    input:
        tuple sample_id, file(x) from fastq_ch
    output:
        file '*_{trimmed,merged}.fastq' into mergetrim_ch
        file 'merge-stats.csv' into mergestats_ch
        file 'primertrim-stats.csv' into trimstats_ch
    script:
    """
    01_merge_and_trim.R ${x}
    """
}

process merge_and_trim_stats {
    publishDir params.outdir, mode: 'copy'

    input:
        file 'merge-stats*' from mergestats_ch.collect()
        file 'primertrim-stats*' from trimstats_ch.collect()
    output:
        file '*stats.csv'
    script:
    """
    cat merge-stats* > merge-reads-stats.csv
    cat primertrim-stats* > primertrim-stats.csv
    """
}

