process multiqc {
    label 'multiqc'
    label 'smallTask'

    if ( params.softlink_results ) { publishDir "${params.output}", pattern: 'multiqc_report.html' }
    else { publishDir "${params.output}", mode: 'copy', pattern: 'multiqc_report.html' }

    input:
    path(config)
    path(fastqcPre)
    path(fastp)
    path(fastqcPost)
    path(kmergenie)
    path(hisat2)
    path(quast_full)
    path(quast_filtered)
    path(mt_result)

    output:
    path "multiqc_report.html"
    path "multiqc_data"

    script:
    """
    multiqc -s . -c ${config}
    """
}

process format_kmergenie_report {
    input:
    path(kmergenie_report)

    output:
    path("*_mqc.html")

    script:
    """
    echo "<!--" > tmp
    echo "id: 'kmergenie'" >> tmp
    echo "id: 'kmergenie'" >> tmp
    echo "section_name: 'KmerGenie'" >> tmp
    echo "-->"  >> tmp
    cat tmp ${kmergenie_report} > ${kmergenie_report.baseName}_mqc.html
    """
}