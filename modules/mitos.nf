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

    if ( params.softlink_results ) { publishDir "${params.output}/${params.mitos_dir}", pattern: "*" }
    else { publishDir "${params.output}/${params.mitos_dir}", mode: 'copy', pattern: "*" }

    input:
    tuple val(assembly_name), val(sequences)
    path(mitos_ref_dir)
    val(genetic_code)

    output:
    path("${assembly_name}")

    script:
    """
    echo "${sequences}" > ${assembly_name}_single.fa # for long sequences can break
    contig=\$(head -n 1 ${assembly_name}_single.fa | awk '{ sub(">", "", \$1); print \$1 }')
    mkdir -p ${assembly_name}/\$contig
    runmitos.py -c ${genetic_code} --refdir . --refseqver ${mitos_ref_dir} -i ${assembly_name}_single.fa -o ${assembly_name}/\$contig
    # rm ${assembly_name}_single.fa
    """
}