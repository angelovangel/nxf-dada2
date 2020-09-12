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

params.readsdir = "fastq"
params.outdir = "${params.readsdir}/results-rexmap"
params.fqpattern = "*_R{1,2}_001.fastq.gz"

params.help = ""

mqc_config = file(params.multiqc_config) // this is needed, otherwise the multiqc config file is not available in the docker image
