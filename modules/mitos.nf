process get_mitos_ref{
    label 'smallTask'
    
    output:
    path('refseq89f')

    script:
    """
    curl https://zenodo.org/record/4284483/files/refseq89f.tar.bz2 -o refseq89f.tar.bz2
    tar -xf refseq89f.tar.bz2
    """
}

process mitos{
    label 'mitos'

    input:
    tuple val(assembly_name), path(assembly)
    path(mitos_ref_dir)
    val(genetic_code)

    output:
    path("${assembly_name}")

    script:
    """
    mkdir ${assembly_name}
    runmitos.py -c ${genetic_code} --refdir . --refseqver ${mitos_ref_dir} -i ${assembly} -o ${assembly_name}
    for RES in `ls soapdenovo2k17-cap3`; do
        contig=\$(head -n 1 soapdenovo2k17-cap3/\$RES/sequence.fas-0 | awk '{sub(">", "", \$1); print \$1}')
        mv soapdenovo2k17-cap3/\$RES soapdenovo2k17-cap3/\$contig
    done
    """
}
