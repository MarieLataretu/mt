process fastqc {
    label 'fastqc'

    if ( params.softlink_results ) { publishDir "${params.output}/${params.fastqc_dir}", pattern: "*_fastqc.zip" }
    else { publishDir "${params.output}/${params.fastqc_dir}", mode: 'copy', pattern: "*_fastqc.zip" }

    input:
    tuple val(name), path(reads), val(mode)

    output:
    path("*_fastqc.zip", emit: zip)

    script:
    """
    fastqc --noextract -t ${task.cpus} ${reads}
    """
}