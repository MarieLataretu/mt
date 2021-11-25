process fastp {
    label 'fastp'

    if ( params.softlink_results ) { publishDir "${params.output}/${params.fastp_dir}" }
    else { publishDir "${params.output}/${params.fastp_dir}", mode: 'copy' }

    input:
    tuple val(name), path(reads), val(mode)

    output:
    tuple val(name), path("${name}*_trimmed*.fastq.gz"), val(mode), emit: sample_trimmed
    tuple val(name), path("${name}_fastp.json"), val(mode), emit: json_report

    script:
    if (mode == 'single') {
    """
    fastp -i ${reads[0]} -o ${name}_trimmed.fastq.gz --thread ${task.cpus} --json ${name}_fastp.json ${params.fastp_additional_params}
    """
    }
    else {
    """
    fastp -i ${reads[0]} -I ${reads[1]} -o ${name}_1_trimmed.fastq.gz -O ${name}_2_trimmed.fastq.gz --thread ${task.cpus} --json ${name}_fastp.json ${params.fastp_additional_params}
    """
    }
}

process get_insert_peak_from_fastp {
    label 'python'
    label 'smallTask'
    
    input:
    tuple val(name), path(json_report), val(mode)

    output:
    tuple val(name), env(insert_peak)

    script:

    if ( mode == 'paired' )
        """
        insert_peak=\$(python3 ${projectDir}/bin/get_insert_peak_from_fastp.py ${json_report})
        """
    else
        error "Single-end reads have no insert size."
}

process get_mean_read_length_from_fastp {
    label 'python'
    label 'smallTask'

    input:
    tuple val(name), path(json_report), val(mode)

    output:
    env(mean_read_length)

    script:
    """
    mean_read_length=\$(python3 ${projectDir}/bin/get_mean_read_len_from_fastp.py ${json_report} ${mode})
    """
}
