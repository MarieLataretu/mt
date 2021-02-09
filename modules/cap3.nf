process cap3{
    label 'cap3'

    input:
    path(assembly)

    output:
    path('*_cap3.fasta')

    script:
    """
    cap3 ${assembly} > cap3.log 2> cap3.log
    cat ${assembly}.cap.contigs ${assembly}.cap.singlets > ${assembly.baseName}_cap3.fasta
    """
}