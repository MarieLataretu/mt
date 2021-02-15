process get_bed{
    label 'hisat2'

    input:
    tuple val(assembly_name), path(bam), path(bam_index)

    output:
    tuple path(bam), path(bam_index), emit: bam
    tuple val(assembly_name), path("${bam.baseName}.bed"), emit: bed

    script:
    """
    samtools view -H ${bam} | grep '^@SQ' | awk -F ':|\t' -v OFS='\t' '{{print \$3, 0, \$5}}' > ${bam.baseName}.bed
    """
    
}

process get_coverage {
    label 'megadepth'

    input:
    tuple path(bam), path(bam_index)
    tuple val(assembly_name), path(bed)

    output:
    tuple val(assembly_name), path("*_coverage.tsv")

    script:
    """
    megadepth ${bam} --annotation ${bed} | awk -F '\\t' -v OFS='\\t' '{print \$1, \$3, \$4/\$3}' > ${bam.baseName.split('.sorted')[0]}_coverage.tsv
    """
}

process get_95th_percentile {
    label 'python'

    input:
    tuple val(assembly_name), path(read_coverage_tsv)

    output:
    tuple val(assembly_name), path("*.tsv")

    script:
    """
    #!/usr/bin/env python3
    import pandas as pd

    df = pd.read_csv("${read_coverage_tsv}", sep='\\t', names=['contig', 'length', 'norm_read_coverage'], header=0)
    print(df.head())
    quantile_95 = df['norm_read_coverage'].quantile(q=.95)
    df.loc[df.norm_read_coverage >= quantile_95, 'in_95_quantile_read_coverage'] = True
    df.loc[df.norm_read_coverage < quantile_95, 'in_95_quantile_read_coverage'] = False   
    df.to_csv("${read_coverage_tsv.baseName}_95th.tsv", sep='\\t', index=False)
    """
}

process pident_filter {
    label 'python'
    
    input:
    tuple val(assembly_name), path(homo_result)
    val(threshold)

    output:
    tuple val(assembly_name), path("${homo_result}.pident")

    script:
    """
    #!/usr/bin/env python3
    import pandas as pd

    df = pd.read_csv("${homo_result}", sep='\\t', header=None, names=['s/qseqid', 'q/sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore', 'qlen', 'slen'])
    df = df[ df['pident'] >= int(${threshold}) ]
    df.to_csv("${homo_result}.pident", sep='\\t', index=False, header=False)
    """
}

process get_features {
    label 'python'
    
    input:
    val(tool)
    tuple val(assembly_name), path(results)

    output:
    tuple val(assembly_name), path("${assembly_name}_${tool}-cov.tsv")

    script:
    result_list =  results.collect{ "\"${it}\"" }
    if ( tool == 'blast')
        """
        #!/usr/bin/env python3
        import pandas as pd

        results = []
        for feature_prot in ${result_list}:
            df = pd.read_csv(feature_prot, sep='\\t', header=None)
            results.append(df)
        df = pd.concat(results, axis=0, ignore_index=True)
        df.columns = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore', 'qlen', 'slen']

        contig_grouped = df.groupby(by='sseqid')

        dict_contig_cov_raw = {}
        dict_contig_cov = {}
        
        for name, group in contig_grouped:
            if name not in dict_contig_cov_raw:
                dict_contig_cov_raw[name] = set()
            for i, start in group['sstart'].iteritems():
                dict_contig_cov_raw[name] = dict_contig_cov_raw[name].union(set(range(min(int(start), int(group.at[i, 'send'])), max(int(start), int(group.at[i, 'send']))+1)))

            dict_contig_cov[str(name)] = [len(dict_contig_cov_raw[name])/group['slen'].iloc[0], len(group)]

        df_cov = pd.DataFrame.from_dict(dict_contig_cov, orient='index', columns=['${tool}_cov', '#${tool}_hits'])
        df_cov.to_csv("${assembly_name}_${tool}-cov.tsv", sep='\\t', index_label='contig')
        """
    else if( tool == 'diamond' )
        """
        #!/usr/bin/env python3
        import pandas as pd

        results = []
        for feature_prot in ${result_list}:
            df = pd.read_csv(feature_prot, sep='\\t', header=None)
            results.append(df)
        df = pd.concat(results, axis=0, ignore_index=True)
        df.columns = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore', 'qlen', 'slen']

        contig_grouped = df.groupby(by='qseqid')

        dict_contig_cov_raw = {}
        dict_contig_cov = {}
        
        for name, group in contig_grouped:
            if name not in dict_contig_cov_raw:
                dict_contig_cov_raw[name] = set()
            for i, start in group['qstart'].iteritems():
                dict_contig_cov_raw[name] = dict_contig_cov_raw[name].union(set(range(min(int(start), int(group.at[i, 'qend'])), max(int(start), int(group.at[i, 'qend']))+1)))

            dict_contig_cov[str(name)] = [len(dict_contig_cov_raw[name])/group['qlen'].iloc[0], len(group)]

        df_cov = pd.DataFrame.from_dict(dict_contig_cov, orient='index', columns=['${tool}_cov', '#${tool}_hits'])
        df_cov.to_csv("${assembly_name}_${tool}-cov.tsv", sep='\\t', index_label='contig')
        """
    else
        error "Unknown tool: ${tool}"
}

process collect_features {
    label 'python'

    input:
    tuple val(assembly_name), path(read_coverage), path(blast_features), path(diamond_features)

    output:
    tuple val(assembly_name), path("${assembly_name}.tsv")

    script:
    """
    #!/usr/bin/env python3
    import pandas as pd

    df_cov = pd.read_csv("${read_coverage}", sep='\\t', index_col='contig')
    df_blast = pd.read_csv("${blast_features}", sep='\\t', index_col='contig')
    df_diamond = pd.read_csv("${diamond_features}", sep='\\t', index_col='contig')

    df_cov.index = df_cov.index.map(str)
    df_blast.index = df_blast.index.map(str)
    df_diamond.index = df_diamond.index.map(str)

    df = df_cov.join(df_blast, on='contig').join(df_diamond, on='contig')
    df.sort_values(by=['norm_read_coverage', 'in_95_quantile_read_coverage', 'length'], ascending=[False, False, False], inplace=True)

    df = pd.concat([df], keys=["${assembly_name}"], names=['assembly'])

    df.to_csv("${assembly_name}.tsv", sep='\t', na_rep='NA')
    """
}

process result_table {
    label 'python'
    echo true

    input:
    path(assembly_result)

    output:
    path("result.tsv")

    script:
    result_list =  assembly_result.collect{ "\"${it}\"" }
    """
    #!/usr/bin/env python3
    import pandas as pd

    results = []
    for assembly_result in ${result_list}:
        df = pd.read_csv(assembly_result, sep='\\t')
        results.append(df)
    df = pd.concat(results, axis=0)
    df.to_csv('result.tsv', sep='\t', na_rep='NA', index=False)
    """
}