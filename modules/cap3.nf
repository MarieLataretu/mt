process cap3{
    label 'cap3'

    input:
    tuple val(assembly_name), path(assembly)

    output:
    tuple val("${assembly_name}-cap3"), path("${assembly_name}-cap3.fasta") optional true

    script:
    """
    if [ -s ${assembly} ]
    then
        cap3 ${assembly} > cap3.log 2> cap3.log
        cat ${assembly}.cap.contigs ${assembly}.cap.singlets > ${assembly_name}-cap3.fasta
    fi
    """
}