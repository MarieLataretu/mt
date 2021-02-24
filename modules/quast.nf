process quast {
    label 'quast'

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