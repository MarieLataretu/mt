process filter_featureProt{
    label 'python'

    input:
    path(featureProt)
    val(genus)
    
    output:
    path("*.filtered")

    script:
    """
    #!/usr/bin/env python3
    from Bio import SeqIO

    filtered_list = []
    for seq_record in SeqIO.parse("${featureProt}", 'fasta'):
        current_genus = seq_record.description.split('[')[1].split(']')[0].split(' ')[0]
        if not "${genus}".lower() == current_genus.lower():
            filtered_list.append(seq_record)

    with open("${featureProt}.filtered", 'w') as out:
        SeqIO.write(filtered_list, out, 'fasta')
    """
}