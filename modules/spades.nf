process spades_input {
    input:
    path(pe_reads)

    output:
    path('spades_input.yml')

    script:
    """
    echo '[' > spades_input.yml
    echo ']' >> spades_input.yml

    """
}

    //     for i in range(0, len(input.pe1)):
    //         out.write('{\n')
    //         out.write('orientation: "fr",\n')
    //         out.write('type: "paired-end",\n')
    //         out.write('right reads: ["' + os.path.abspath(input.pe1[i]) + '"],\n')
    //         out.write('left reads: ["' + os.path.abspath(input.pe2[i]) + '"]\n')
    //         out.write('}\n')
    //         if(i < len(input.pe1)-1 or len(input.se) > 0):
    //             out.write(',\n')
    //     for i in range(0, len(input.se)):
    //         out.write('{\n')
    //         out.write('type: "single",\n')
    //         out.write('single reads: ["' + os.path.abspath(input.se[i]) + '"],\n')
    //         out.write('}\n')
    //         if(i < len(input.se)-1):
    //             out.write(',')

process spades {
    label 'spades'

    input:
    path(input)

    output:
    path('spades/scaffolds.fasta')

    script:
    """
    spades.py -o spades -t ${task.cpus} --disable-gzip-output --isolate --dataset ${input}
    # removes spades assembly meta data        
    rm -rf K*
    """
}