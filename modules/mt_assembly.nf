process extract_contigs {
    label 'python'
    label 'smallTask'

    if ( params.softlink_results ) { publishDir "${params.output}/${params.mtContigs_dir}", pattern: '*' }
    else { publishDir "${params.output}/${params.mtContigs_dir}", mode: 'copy', pattern: '*' }


    input:
    tuple val(assembly_name), path(full_assembly), path(contig_list)

    output:
    tuple val(assembly_name), path("${assembly_name}.filtered.fasta")

    script:
    """
    #!/usr/bin/env python3
    from Bio import SeqIO
    
    contig_ids = []
    with open("${contig_list}", 'r') as tsv:
        for line in tsv:
            contig_ids.append(line.strip())
    contigs = []
    for seq_record in SeqIO.parse("${full_assembly}", 'fasta'):
        id = seq_record.description.split(' ')[0].strip('>')
        if id in contig_ids:
            contigs.append(seq_record)
    SeqIO.write(contigs, "${assembly_name}.filtered.fasta", 'fasta')
    """
}