process fastqc {
    label 'fastqc'

    if ( params.softlink_results ) { publishDir "${params.output}/${params.fastqc_dir}" }
    else { publishDir "${params.output}/${params.fastqc_dir}", mode: 'copy' }

    input:
    tuple val(name), path(reads), val(mode)

    output:
    path("*_fastqc.zip", emit: zip)
    path("*_fastqc.html", emit: html)

    script:
    """
    fastqc --noextract -t ${task.cpus} ${reads}
    """
}