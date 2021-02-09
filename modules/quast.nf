process quast {
    label 'quast'

    if ( params.softlink_results ) { publishDir "${params.output}/${params.quast_dir}", pattern: 'quast_results/latest/{report,icarus,icarus_viewers/contig_size_viewer}.{pdf,html}' }
    else { publishDir "${params.output}/${params.quast_dir}", mode: 'copy', pattern: 'quast_results/latest/{report,icarus,icarus_viewers/contig_size_viewer}.{pdf,html}' }

    input:
    path(assemblies)

    output:
    path('quast_results/latest/report.tsv'), emit: report_tsv
    tuple path('quast_results/latest/report.pdf'), path('quast_results/latest/report.html'), path('quast_results/latest/icarus.html'), path('quast_results/latest/icarus_viewers/contig_size_viewer.html'), emit: report

    script:
    """
    quast.py -t ${task.cpus} ${assemblies}
    """
}