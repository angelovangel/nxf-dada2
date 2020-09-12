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

// channel for reads folder --> seqtools
Channel
    .fromFilePairs( fastqfiles, checkIfExists: true )
    .set { fastq_ch }

process merge_reads {
    tag "$sample_id"
    publishDir "${params.outdir}/merged_reads", mode: 'copy'
    
    input:
        tuple sample_id, file(x) from fastq_ch
    
    output:
        file '*_merged.fastq' into merged_reads_ch
    
    script:
    """
    01_merge_reads.R ${x}
    """

}
