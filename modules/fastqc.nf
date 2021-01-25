process fastqc {
    label 'fastqc'

    input:
    tuple val(name), path(reads), val(mode)

    output:
    path("*_fastqc.zip", emit: zip)

    script:
    """
    fastqc --noextract -t ${task.cpus} ${reads}
    """
}