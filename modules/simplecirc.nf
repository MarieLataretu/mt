process simplecirc {
    label 'simplecirc'
    label 'python'
    // publishDir "${params.output}/${name}/${params.assemblydir}/circularized_contigs", mode: 'copy', pattern: "*_circ.fasta"
    
    input:
    tuple val(name), file(assembly) 
    
    output:
    tuple val(name), file("*_circ.fasta") optional true
    
    script:
    """
    git clone https://github.com/Kzra/Simple-Circularise.git
    BN=\$(basename ${assembly} .fasta)
    python Simple-Circularise/simple_circularise.py ${assembly} "\$BN"_circ.fasta
    id=\$(head -1 ${assembly})
    sed -i "s/>/\$id/g" "\$BN"_circ.fasta

    # the simple-circularise script might produce a fasta file with no sequence and only 'Blank' written
    # if so, use the original fasta 
    if grep -q Blank "\$BN"_circ.fasta; then
        echo 'nothing found'
        cp ${assembly} "\$BN"_circ.fasta
    fi
    """
}

/*
Usage:

python simple_circularise.py [input.fasta] [output.fasta] [-p] [-min] [-max] [-r]
[-min]: set a minimum size for the output sequence.

[-max]: set a maximum size for the output sequence.

[-p]: set the probability threshold for determining repeat size (default 0.005).

[-r]: change the behaviour of the program to maximise repeat size (default is to maximise output sequence size).

Examples:

python simple_circularise.py linear.fasta circular.fasta -p 0.0001 -min 1000 -max 2000
Circularise a genome with a repeat size determined by a probability of co-occurance < 0.0001. Output the largest sequence between 1 - 2kb that can be circularised.

python simple_circularise.py linear.fasta circular.fasta -r 10 -min 1000
Circularise a genome using largest possible repeat size. Start searching at repeat size 10 and increase until largest is found. Ensure output sequence is >1kb.

python simple_circularise.py linear.fasta circular.fasta 
Circularise a genome based on a repeat size determined by p < 0.005. Output the largest sequence that can be circularised (default behaviour).
*/