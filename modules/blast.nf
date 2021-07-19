process make_blast_db{
    label 'blast'
    
    input:
    path(assembly)

    output:
    tuple path(assembly), path("*.n??")

    script:
    """
    makeblastdb -in ${assembly} -input_type fasta -dbtype nucl
    """
}

process blast{
    label 'blast'

    if ( params.softlink_results ) { publishDir "${params.output}/${params.blast_dir}", pattern: "*.blast" }
    else { publishDir "${params.output}/${params.blast_dir}", mode: 'copy', pattern: "*.blast" }

    input:
    each path(featureProt)
    tuple path(assembly), path(assembly_blast_db)
    val(genetic_code)
    
    output:
    tuple val(assembly.baseName), path("*.blast")

    script:
    """
    tblastn -num_threads ${task.cpus} -query ${featureProt} -db ${assembly} -db_gencode ${genetic_code} -evalue 1e-10 -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen" > ${featureProt.baseName}_${assembly.baseName}.blast
    """
}