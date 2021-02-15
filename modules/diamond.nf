process make_diamond_db{
    label 'diamond'

    input:
    path(featureProt)

    output:
    path("${featureProt.baseName}.dmnd")

    script:
    """
    diamond makedb --in ${featureProt} --db ${featureProt.baseName} --threads ${task.cpus} --log
    """
}

process diamond{
    label 'diamond'

    input:
    path(featureProt_db)
    each path(assembly)

    output:
    tuple val(assembly.baseName), path("*.diamond")

    script:
    """
    diamond blastx --query ${assembly} --db ${featureProt_db} --out ${assembly.baseName}_${featureProt_db.baseName}.diamond --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen --threads ${task.cpus} --log
    # --evalue default 0.001
    """
}