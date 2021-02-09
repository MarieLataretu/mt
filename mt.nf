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

// Set valid_params = ['max_cores', 'cores', 'memory', 'profile', 'help', 'reads', 'genome', 'annotation', 'deg', 'autodownload', 'pathway', 'species', 'include_species', 'strand', 'mode', 'tpm', 'fastp_additional_params', 'hisat2_additional_params', 'featurecounts_additional_params', 'feature_id_type', 'busco_db', 'dammit_uniref90', 'skip_sortmerna', 'assembly', 'output', 'fastp_dir', 'sortmerna_dir', 'hisat2_dir', 'featurecounts_dir', 'tpm_filter_dir', 'annotation_dir', 'deseq2_dir', 'assembly_dir', 'rnaseq_annotation_dir', 'uniref90_dir', 'multiqc_dir', 'nf_runinfo_dir', 'permanentCacheDir', 'condaCacheDir', 'singularityCacheDir', 'softlink_results', 'cloudProcess', 'permanent-cache-dir', 'conda-cache-dir', 'singularity-cache-dir', 'cloud-process'] // don't ask me why there is 'permanent-cache-dir', 'conda-cache-dir', 'singularity-cache-dir', 'cloud-process'
// def parameter_diff = params.keySet() - valid_params
// if (parameter_diff.size() != 0){
//     exit 1, "ERROR: Parameter(s) $parameter_diff is/are not valid in the pipeline!\n"
// }

// terminal prints
if (params.help) { exit 0, helpMSG() }

include { fastqc as fastqcPre; fastqc as fastqcPost } from './modules/fastqc'
include { fastp; get_insert_peak_from_fastp; get_mean_read_length_from_fastp } from './modules/fastp'
include { spades_input; spades } from './modules/spades'
include { kmergenie_input; kmergenie; soapdenovo2_input; soapdenovo2 } from './modules/soapdenovo2'
include { cap3 } from './modules/cap3'
include { hisat2index; hisat2; index_bam } from './modules/hisat2'
include { quast } from './modules/quast'
include { multiqc } from './modules/multiqc'

if ( params.genus ) {
    genus_ch = Channel.value( params.genus )
}
if ( params.pe_reads ) {
    paired_reads_ch = Channel.fromFilePairs( params.pe_reads, size: 2, checkIfExists: true ).map {it -> it + ['paired']}
}
if ( params.se_reads ) {
    single_reads_ch = Channel.fromFilePairs( params.se_reads, size: 1, checkIfExists: true ).map {it -> it + ['single']}
}

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

    all_trimmed_paired_read_paths = trimmed_paired_reads.map{it -> it[1]}.collect()
    all_trimmed_single_read_paths = trimmed_single_reads.map{it -> it[1]}.collect()
    all_trimmed_read_paths = all_trimmed_paired_read_paths.concat(all_trimmed_single_read_paths).collect()
    
    // assemblies
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
    soapdenovo2_input(trimmed_paired_reads.map{it -> it[1]+it[3]}.collect(), all_trimmed_single_read_paths)
    soapdenovo2(kmers, soapdenovo2_input.out, all_trimmed_read_paths)

    assemblies = spades.out.concat(soapdenovo2.out)

    cap3(assemblies)

    // map reads back to assembly
    hisat2index(assemblies)
    hisat2(trimmed_paired_reads.map{it -> it[1][0]}.collect(), trimmed_paired_reads.map{it -> it[1][1]}.collect(), all_trimmed_single_read_paths, hisat2index.out, params.hisat2_additional_params)
    index_bam(hisat2.out.sample_bam)

    // summary
    quast(assemblies.mix(cap3.out).collect())

    multiqc(fastqcPre.out.collect(), fastp.out.json_report.map{ it -> it[1] }.collect(), fastqcPost.out.collect(), kmergenie.out.report, hisat2.out.log.collect(), quast.out.report_tsv)

    // blasten
    // ranken & filtern
    
    // mit ref: vgl
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