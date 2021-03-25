#!/usr/bin/env nextflow

XX = "20"
YY = "07"
ZZ = "1"

if ( nextflow.version.toString().tokenize('.')[0].toInteger() < XX.toInteger() ) {
println "\033[0;33mRNAflow requires at least Nextflow version " + XX + "." + YY + "." + ZZ + " -- You are using version $nextflow.version\u001B[0m"
exit 1
}
else if ( nextflow.version.toString().tokenize('.')[1].toInteger() < YY.toInteger() ) {
println "\033[0;33mRNAflow requires at least Nextflow version " + XX + "." + YY + "." + ZZ + " -- You are using version $nextflow.version\u001B[0m"
exit 1
}

nextflow.enable.dsl=2

// Parameters sanity checking

Set valid_params = ['max_cores', 'cores', 'memory', 'profile', 'help', 'genus', 'se_reads', 'pe_reads', 'reference_genome', 'reference_annotation', 'fastp_additional_params', 'hisat2_additional_params', 'genetic_code', 'contig_len_filter', 'contig_high_read_cov_filter', 'output', 'condaCacheDir', 'softlink_results', 'conda-cache-dir'] // don't ask me why there is 'conda-cache-dir'
def parameter_diff = params.keySet() - valid_params
if (parameter_diff.size() != 0){
    exit 1, "ERROR: Parameter(s) $parameter_diff is/are not valid in the pipeline!\n"
}

// terminal prints
if (params.help) { exit 0, helpMSG() }
if ( params.profile ) { exit 1, "--profile is WRONG use -profile" }
if ( ! params.pe_reads && ! params.se_reads ) { exit 1, "Read data is required." }
if ( params.reference_annotation && ! params.reference_genome ) { exit 1, "If an annotation file is provided, a reference file needs to be there, too." }

include { fastqc as fastqcPre; fastqc as fastqcPost } from './modules/fastqc'
include { fastp; get_insert_peak_from_fastp; get_mean_read_length_from_fastp } from './modules/fastp'
include { spades_input; spades } from './modules/spades'
include { kmergenie_input; kmergenie; soapdenovo2_input; soapdenovo2 } from './modules/soapdenovo2'
include { cap3 } from './modules/cap3'
include { hisat2index; hisat2; index_bam } from './modules/hisat2'
include { filter_featureProt } from './modules/filter_featureProt'
include { make_blast_db; blast } from './modules/blast'
include { mmseqs2_create_target_db_index ; mmseqs2_search } from './modules/mmseqs2'
include { get_bed; get_coverage; get_95th_percentile; pident_filter as blast_pident_filter; pident_filter as mmseqs2_pident_filter ; get_features as get_blast_features; get_features as get_mmseqs2_features; collect_features; result_table } from './modules/features'
include { extract_contigs } from './modules/mt_assembly'
include { get_mitos_ref; mitos as mitos; mitos as mitos_ref } from './modules/mitos'
include { quast as quast_complete_assembly; quast as quast_mt_assemblys } from './modules/quast'
include { format_kmergenie_report; multiqc } from './modules/multiqc'

if ( params.pe_reads ) {
    paired_reads_ch = Channel.fromFilePairs( params.pe_reads, size: 2, checkIfExists: true ).map {it -> it + ['paired']}
} else {
    paired_reads_ch = Channel.empty()
}
if ( params.se_reads ) {
    single_reads_ch = Channel.fromFilePairs( params.se_reads, size: 1, checkIfExists: true ).map {it -> it + ['single']}
} else {
    single_reads_ch = Channel.empty()
}

reference_genome = params.reference_genome ? Channel.fromPath( params.reference_genome, checkIfExists: true ) : file( "${params.output}/no_ref_genome" )
reference_annotation = params.reference_annotation ? Channel.fromPath( params.reference_annotation, checkIfExists: true ) : file( "${params.output}/no_ref_annotation" )

featureProt_ch = Channel.fromPath( workflow.projectDir + '/assets/featureProt/*.faa', checkIfExists: true )
multiqc_config = Channel.fromPath( workflow.projectDir + '/assets/multiqc_config.yml', checkIfExists: true )

def get_mean_two_third_read_length (mean_read_lengths) {
    mean_read_len = ( mean_read_lengths.sum() / mean_read_lengths.count() )
    mean_read_len.map{ it -> Math.round((2/3) * it) }.map{ it % 2 == 0 ? it+1 : it }
}

workflow {
    // preprocessing
    fastqcPre(paired_reads_ch.concat(single_reads_ch))
    fastp(paired_reads_ch.concat(single_reads_ch))
    fastqcPost(fastp.out.sample_trimmed)

    // get insert peak for paired-end data
    get_insert_peak_from_fastp(fastp.out.json_report.filter { it[2] == 'paired' })
    // join insert peak value to read Channel
    trimmed_paired_reads = fastp.out.sample_trimmed.filter { it[2] == 'paired' }.join(get_insert_peak_from_fastp.out, by: [0,0])
    trimmed_single_reads = fastp.out.sample_trimmed.filter { it[2] == 'single' }

    all_trimmed_paired_read_paths = trimmed_paired_reads.map{it -> it[1]}.collect().ifEmpty { file( "${params.output}/EMPTY") }
    all_trimmed_single_read_paths = trimmed_single_reads.map{it -> it[1]}.collect().ifEmpty { file( "${params.output}/EMPTY") }
    all_trimmed_read_paths = trimmed_paired_reads.map{it -> it[1]}.collect().concat(trimmed_single_reads.map{it -> it[1]}.collect()).collect()

    // Assemblies
    // SPAdes
    spades_input(all_trimmed_paired_read_paths, all_trimmed_single_read_paths)
    spades(spades_input.out, all_trimmed_read_paths)
    
    // get kmers for SOAPdenovo2
    kmergenie_input(all_trimmed_read_paths)
    kmergenie(kmergenie_input.out, all_trimmed_read_paths)

    get_mean_read_length_from_fastp(fastp.out.json_report)
    mean_read_len = get_mean_two_third_read_length(get_mean_read_length_from_fastp.out.toFloat())
    kmers = mean_read_len.concat(kmergenie.out.best_kmer).collect().map { it.unique() }
    // SOAPdenovo2
    soapdenovo2_input(trimmed_paired_reads.map{it -> it[1]+it[3]}.collect().ifEmpty { [file( "${params.output}/EMPTY")] }, all_trimmed_single_read_paths)
    soapdenovo2(kmers, soapdenovo2_input.out, all_trimmed_read_paths)

    assemblies = spades.out.concat(soapdenovo2.out)

    // Scaffolding
    cap3(assemblies)

    assemblies_scaffolds = assemblies.concat(cap3.out)

    // map reads back to assembly
    hisat2index(assemblies_scaffolds.map{it -> it[1]})
    hisat2(trimmed_paired_reads.map{it -> it[1][0]}.collect().ifEmpty { file( "${params.output}/EMPTY1")}, trimmed_paired_reads.map{it -> it[1][1]}.collect().ifEmpty { file( "${params.output}/EMPTY2")}, all_trimmed_single_read_paths, hisat2index.out, params.hisat2_additional_params)
    index_bam(hisat2.out.sample_bam)

    // Read coverage
    get_bed(index_bam.out)
    get_coverage(get_bed.out.bam, get_bed.out.bed)
    get_95th_percentile(get_coverage.out)

    // filter feature proteins (exclude genus)
    if ( params.genus ){
        filter_featureProt(featureProt_ch, params.genus)
        featureProt_filtered = filter_featureProt.out
    } else {
        featureProt_filtered = featureProt_ch
    }

    // blast
    make_blast_db(assemblies_scaffolds.map{it -> it[1]})
    blast(featureProt_filtered.collect(), make_blast_db.out, params.genetic_code)
    // mmseqs2
    mmseqs2_create_target_db_index(assemblies_scaffolds.map{it -> it[1]})
    mmseqs2_search(featureProt_filtered.collect(), mmseqs2_create_target_db_index.out, params.genetic_code)

    // blast features
    blast_pident_filter(blast.out, 70)
    get_blast_features('blast', blast_pident_filter.out.groupTuple())

    // mmseqs2 features
    mmseqs2_pident_filter(mmseqs2_search.out, 70)
    get_mmseqs2_features('mmseqs2', blast_pident_filter.out.groupTuple())

    // collect features
    collect_features(get_95th_percentile.out.join(get_blast_features.out.join(get_mmseqs2_features.out)), params.contig_len_filter, params.contig_high_read_cov_filter == 'true' ? 'True' : 'False')
    result_table(collect_features.out.table.map{ it -> it[1] }.collect())

    extract_contigs(assemblies_scaffolds.join(collect_features.out.contigs))

    // Annotate
    get_mitos_ref()
    mitos(extract_contigs.out, get_mitos_ref.out, params.genetic_code)
    if ( params.reference_genome ) {
        mitos_ref(reference_genome.map{ it -> [it.baseName, it] }, get_mitos_ref.out, params.genetic_code)
    }

    // Summary
    // QUAST full assembly
    quast_complete_assembly('full_assembly', assemblies_scaffolds.map{it -> it[1]}.collect(), file( "${params.output}/no_ref_genome" ), file( "${params.output}/no_ref_annotation"))
    // QUAST filtered assembly
    quast_mt_assemblys('filtered_assembly', extract_contigs.out.map{it -> it[1]}.collect(), reference_genome, reference_annotation)

    // format stuff for MultiQC
    format_kmergenie_report(kmergenie.out.report)
    // run MultiQC
    multiqc(multiqc_config, fastqcPre.out.collect(), fastp.out.json_report.map{ it -> it[1] }.collect(), fastqcPost.out.collect(), format_kmergenie_report.out, hisat2.out.log.collect(), quast_complete_assembly.out.report_tsv, quast_mt_assemblys.out.report_tsv,  result_table.out)
}

def helpMSG() {
    c_green = "\033[0;32m";
    c_reset = "\033[0m";
    c_yellow = "\033[0;33m";
    c_blue = "\033[0;34m";
    c_dim = "\033[2m";
    log.info """
    """.stripIndent()
}