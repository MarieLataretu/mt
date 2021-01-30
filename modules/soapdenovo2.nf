process kmergenie_input {
    input:
    path(reads)

    output:
    path('kmergenie.config')

    script:
    """
    for READ_FILE in ${reads}; do
        echo \$READ_FILE >> kmergenie.config;
    done
    """
}
process kmergenie {
    label 'kmergenie'

    input:
    path(input)
    path(reads)

    output:
    env(best_kmer)

    script:
    """
    kmergenie ${input} -t ${task.cpus} > kmergenie.log
    best_kmer=\$(grep '^best k: ' kmergenie.log | awk '{print \$3}')
    """
}

process soapdenovo2_input {
    label 'python'

    input:
    val(pe_reads_and_info)
    path(se_reads)

    output:
    path('soapdenovo2.config')

    script:
    // make lists of strings for the python script
    pe_reads_and_info_py =  pe_reads_and_info.collect{ "\"${it}\"" }
    se_reads_py = se_reads.collect { "\"${it}\"" }
    """
    #!/usr/bin/env python3

    with open('soapdenovo2.config', 'a') as txt:
        for pe_index in range(0, len(${pe_reads_and_info_py}), 3):
            txt.write('LIB')
            insert_size = ${pe_reads_and_info_py}[pe_index + 2]
            txt.write(f"avg_ins={insert_size}")
            R1 = ${pe_reads_and_info_py}[pe_index]
            txt.write(f"q1={R1}")
            R2 = ${pe_reads_and_info_py}[pe_index + 1]
            txt.write(f"q2={R2}")
        for se in ${se_reads_py}:
            txt.write('LIB')
            txt.write(f"q={se}")
    """
}



process soapdenovo2 {
    label 'soapdenovo2'

    input:
    path(input_yaml)
    path(pe_reads)
    path(se_reads)

    output:
    path('spades/scaffolds.fasta')

    script:
    """
    spades.py -o spades -t ${task.cpus} --disable-gzip-output --isolate --dataset ${input_yaml}
    # removes spades assembly meta data        
    rm -rf K*
    """
}