process fastp {
    label 'fastp'

    // if ( params.softlink_results ) { publishDir "${params.output}/${params.fastp_dir}", pattern: "*.trimmed.fastq.gz" }
    // else { publishDir "${params.output}/${params.fastp_dir}", mode: 'copy', pattern: "*.trimmed.fastq.gz" }

    input:
    tuple val(name), path(reads), val(mode)

    output:
    tuple val(name), path("${name}*.trimmed.fastq.gz"), val(mode), emit: sample_trimmed
    tuple val(name), path("${name}_fastp.json"), val(mode), emit: json_report

    script:
    if (mode == 'single') {
    """
    fastp -i ${reads[0]} -o ${name}.trimmed.fastq.gz --thread ${task.cpus} --json ${name}_fastp.json ${params.fastp_additional_params}
    """
    }
    else {
    """
    fastp -i ${reads[0]} -I ${reads[1]} -o ${name}.R1.trimmed.fastq.gz -O ${name}.R2.trimmed.fastq.gz --thread ${task.cpus} --json ${name}_fastp.json ${params.fastp_additional_params}
    """
    }
}

process get_insert_peak_from_fastp{
    input:
    tuple val(name), path(json_report), val(mode)

    output:
    tuple val(name), env(insert_peak)

    script:
    """
    if [ "\$(grep -c '"insert_size"' ${json_report})" -eq 1 ]; then
        if [ "\$(grep -cP '"peak": [0-9]{1,},' ${json_report})" -eq 1 ]; then
            insert_peak=\$(grep -P '"peak": [0-9]{1,},' ${json_report} | grep -oP '[0-9]{1,}')
        else
            echo "ERROR: Failed to get the insert peak value. No or no unambiguous match for `"peak": [0-9]{1,},`" >&2
            exit 1
        fi
    else
        echo "ERROR: Failed to get the insert peak value. No or no unambiguous match for `"insert_size"`." >&2
        exit 1
    fi
    """
}
