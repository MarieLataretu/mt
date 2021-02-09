process multiqc {
    label 'multiqc'
    label 'smallTask'

    if ( params.softlink_results ) { publishDir "${params.output}/${params.multiqc_dir}", pattern: 'multiqc_report.html' }
    else { publishDir "${params.output}/${params.multiqc_dir}", mode: 'copy', pattern: 'multiqc_report.html' }

    input:
    path(fastqcPre)
    path(fastp)
    path(fastqcPost)
    path(kmergenie)
    path(hisat2)
    path(quast)

    output:
    path "multiqc_report.html"
    path "multiqc_data"

    script:
    """
    multiqc . --cl_config "{ sp: { kmergenie: {fn: 'histograms_report_mqc.html'} } }"
    """
}