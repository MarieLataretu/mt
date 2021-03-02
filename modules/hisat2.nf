/************************************************************************
* HISAT2 INDEX
************************************************************************/
process hisat2index {
    label 'hisat2'

    input:
    path(reference)

    output:
    tuple path(reference), path("${reference.baseName}*.ht2")

    script:
    """
    hisat2-build -p ${task.cpus} ${reference} ${reference.baseName}
    """
}


/************************************************************************
* HISAT2
************************************************************************/
process hisat2 {
    label 'hisat2'
    
    // if ( params.softlink_results ) { publishDir "${params.output}/${params.hisat2_dir}", pattern: "*.sorted.bam" }
    // else { publishDir "${params.output}/${params.hisat2_dir}", mode: 'copy', pattern: "*.sorted.bam" }

    input:
    path(pe1_reads)
    path(pe2_reads)
    path(se_reads)
    tuple path(reference), path(index)
    val(additionalParams)

    output:
    tuple val(reference.baseName), path("${reference.baseName}.sorted.bam"), emit: sample_bam 
    path("${reference.baseName}_summary.log"), emit: log

    script:
    PE1 = "${pe1_reads}".size() > 1 ? '-1 ' + pe1_reads.join(',') : ''
    PE2 = "${pe2_reads}".size() > 1 ? '-2 ' + pe2_reads.join(',') : ''
    SE = "${se_reads}".size() > 1 ? '-U ' + se_reads.join(',') : ''
    """
    hisat2 -x ${reference.baseName} ${PE1} ${PE2} ${SE} -p ${task.cpus} --new-summary --summary-file ${reference.baseName}_summary.log ${additionalParams} | samtools view -bS | samtools sort -o ${reference.baseName}.sorted.bam -T tmp --threads ${task.cpus}
    """
}
process index_bam {
    label 'hisat2'
    label 'smallTask'    
    
    // if ( params.softlink_results ) { publishDir "${params.output}/${params.hisat2_dir}", pattern: "*.bai" }
    // else { publishDir "${params.output}/${params.hisat2_dir}", mode: 'copy', pattern: "*.bai" }

    input:
    tuple val(assembly), path(bam_file)

    output:
    tuple val(assembly), path(bam_file), path("${bam_file}.bai")

    script:
    """
    samtools index ${bam_file}
    """
}