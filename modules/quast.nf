process quast {
    label 'quast'

    if ( params.softlink_results ) { publishDir "${params.output}", pattern: "quast_*" }
    else { publishDir "${params.output}", mode: 'copy', pattern: "quast_*" }

    input:
    val(name)
    path(assemblies)
    path(reference)
    path(annotation)

    output:
    path("${name}_report.tsv"), emit: report_tsv
    path("quast_${name}"), emit: output

    script:
    def ref_fn = "${reference.simpleName}" != 'no_ref_genome' ? "-r ${reference}" : ''
    def ref_gtf = "${annotation.simpleName}" != 'no_ref_annotation' ? "-g ${annotation}" : ''
    """
    quast.py -o quast_${name} -t ${task.cpus} ${assemblies} ${ref_fn} ${ref_gtf}
    cp quast_${name}/report.tsv ${name}_report.tsv
    """
}