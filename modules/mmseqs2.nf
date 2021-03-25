process mmseqs2_create_target_db_index {
    label 'mmseqs2'

    input:
    path(fasta)

    output:
    tuple val("${fasta.baseName}"), path("*")

    script:
    """
    mmseqs createdb ${fasta} targetDB
    mmseqs createindex targetDB tmp --search-type 2 --threads ${task.cpus}
    rm -rf tmp
    """
}

process mmseqs2_search {
    label 'mmseqs2'

    input:
    each path(query)
    tuple val(target_name), path(target_db)
    val(genetic_code)

    output:
    tuple val(query.baseName), path("*.mmseqs2")

    script:
    """
    mmseqs createdb ${query} queryDB
    mmseqs search queryDB targetDB ${query.baseName}_${target_name}DB tmp --remove-tmp-files --translation-table ${genetic_code} -e 1 --e-profile 1e-10 --threads ${task.cpus}
    mmseqs convertalis queryDB targetDB ${query.baseName}_${target_name}DB ${query.baseName}_${target_name}.mmseqs2 --format-mode 0 --format-output 'query,target,pident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,qlen,tlen' --threads ${task.cpus}
    rm -rf tmp
    """
}
