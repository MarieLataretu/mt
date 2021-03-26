process mmseqs2_create_target_db_index {
    label 'mmseqs2'

    input:
    path(assembly)

    output:
    tuple path(assembly), path("*")

    script:
    """
    mmseqs createdb ${assembly} ${assembly.baseName}DB
    mmseqs createindex ${assembly.baseName}DB tmp --search-type 2 --threads ${task.cpus}
    rm -rf tmp
    """
}

process mmseqs2_search {
    label 'mmseqs2'

    input:
    each path(query)
    tuple path(assembly), path(assembly_db)
    val(genetic_code)

    output:
    tuple val(assembly.baseName), path("*.mmseqs2")

    script:
    """
    mmseqs createdb ${query} queryDB
    mmseqs search queryDB ${assembly.baseName}DB ${query.baseName}_${assembly.baseName}DB tmp --remove-tmp-files --translation-table ${genetic_code} -e 1 --e-profile 1e-10 -s 7.5 --search-type 4 --threads ${task.cpus}
    mmseqs convertalis queryDB ${assembly.baseName}DB ${query.baseName}_${assembly.baseName}DB ${query.baseName}_${assembly.baseName}.mmseqs2 --format-mode 0 --format-output 'query,target,pident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,qlen,tlen' --translation-table 4 --search-type 4 --threads ${task.cpus}
    rm -rf tmp
    """
}
