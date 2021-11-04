process spades_input {
    label 'python'
    label 'smallTask'

    input:
    path(pe_reads)
    path(se_reads)

    output:
    path('spades_input.yml')

    script:
    // make lists of strings for the python script, if input is not empty
    pe_reads_py = pe_reads.baseName == 'EMPTY' ? '[]' : pe_reads.collect{ "\"${it}\"" }
    se_reads_py = se_reads.baseName == 'EMPTY' ? '[]' : se_reads.collect{ "\"${it}\"" }
    """
    #!/usr/bin/env python3
    import os
    import yaml
    read_list = []
    for pe_index in range(0, len(${pe_reads_py}), 2):
        read_list.append(({'type': 'paired-end', 'orientation': 'fr', 'right reads': [os.path.abspath(${pe_reads_py}[pe_index])], 'left reads': [os.path.abspath(${pe_reads_py}[pe_index+1])]}))
    for se in ${se_reads_py}:
        read_list.append(({'single reads': [os.path.abspath(se)], 'type': 'single'}))
    with open('spades_input.yml', 'w') as yml:
        yaml.dump(read_list, yml)
    """
}

process spades {
    label 'spades'

    if ( params.softlink_results ) { publishDir "${params.output}/${params.spades_dir}", pattern: "spades/spades.fasta" }
    else { publishDir "${params.output}/${params.spades_dir}", mode: 'copy', pattern: "spades/spades.fasta" }

    input:
    path(input_yaml)
    path(reads)

    output:
    tuple val('spades'), path('spades/spades.fasta')

    script:
    """
    spades.py -o spades -t ${task.cpus} --disable-gzip-output --isolate --dataset ${input_yaml} --memory ${task.memory.toGiga()}
    # removes spades assembly meta data        
    rm -rf K*
    mv spades/scaffolds.fasta spades/spades.fasta
    """
}

process spades_plasmid {
    label 'spades'

    if ( params.softlink_results ) { publishDir "${params.output}/plasmid/${params.spades_dir}", pattern: "spades/spades_plasmid.fasta" }
    else { publishDir "${params.output}/plasmid/${params.spades_dir}", mode: 'copy', pattern: "spades/spades_plasmid.fasta" }

    input:
    path(input_yaml)
    path(reads)

    output:
    tuple val('spades'), path('spades/spades_plasmid.fasta')

    script:
    """
    spades.py --plasmid -o spades -t ${task.cpus} --disable-gzip-output --dataset ${input_yaml} --memory ${task.memory.toGiga()}
    # removes spades assembly meta data        
    rm -rf K*
    mv spades/contigs.fasta spades/spades_plasmid.fasta
    """
}