process get_bed{
    label 'hisat2'

    input:
    tuple path(bam), path(bam_index)

    output:
    tuple path(bam), path(bam_index), emit: bam
    path("${bam.baseName}.bed"), emit: bed

    script:
    """
    samtools view -H ${bam} | grep '^@SQ' | awk -F ':|\t' -v OFS='\t' '{{print \$3, 0, \$5}}' > ${bam.baseName}.bed
    """
    
}

process get_coverage {
    label 'megadepth'

    input:
    tuple path(bam), path(bam_index)
    path(bed)

    output:
    path("*_coverage.tsv")

    script:
    """
    megadepth ${bam} --annotation ${bed} | awk -F '\\t' -v OFS='\\t' '{print \$0, \$4/\$3}' > ${bam.baseName.strip('.sorted')}_coverage.tsv
    """
}