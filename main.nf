#!/usr/bin/env nextflow

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ALMA - AMF LSU METABARCODING ANALYSIS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Includes: Quality control, taxonomic classification, alignment, phylogenetic placement
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

nextflow.enable.dsl = 2

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PARAMETER VALIDATION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Input parameters
params.input              = null  // CSV samplesheet OR directory with FASTQ files
params.input_dir          = null  // Alternative: directory with *_R{1,2}.fastq.gz files
params.input_pattern      = '*_R{1,2}_001.fastq.gz'  // Pattern for paired-end files
params.outdir             = './results'

// Reference databases
params.reference_db       = null  // Mothur-formatted reference database
params.taxonomy_db        = null  // Taxonomy file for classification
params.alignment_ref      = null  // Reference alignment file
params.reference_mafrax   = null  // Reference alignmente for mafft and Raxml (without problematic sequences)
params.phylo_tree         = null  // Reference phylogenetic tree
params.taxonomy_ALMA      = null  // Taxonomy dictionary for GAPPA

// Mothur parameters
params.min_length         = 200
params.max_length         = 600
params.max_ambig          = 18
params.max_homop          = 22
params.classify_cutoff    = 90
params.align_start        = 4
params.align_end          = 530
params.cluster_cutoff     = 0.02
params.cluster_method     = 'dgc'
params.precluster_diffs = 3
params.cutoff           = 1  // For split.abund

// Post-Mothur parameters
params.gappa_lwr_threshold = 0.9

// Workflow control
params.skip_chimera_check = false
params.skip_phylogenetic  = false

// Help message
params.help = false

def helpMessage() {
    log.info"""
    =========================================
      ALMA - AMF LSU METABARCODING ANALYSIS
    =========================================
    
    Usage (Option 1 - CSV):
      nextflow run main.nf --input samplesheet.csv --reference_db refs.fasta --taxonomy_db taxonomy.txt
    
    Usage (Option 2 - Directory):
      nextflow run main.nf --input_dir data/ --reference_db refs.fasta --taxonomy_db taxonomy.txt
    
    Required Arguments (choose ONE input method):
      --input              Path to samplesheet CSV with columns: sample,fastq_1,fastq_2
                           OR
      --input_dir          Directory containing paired FASTQ files (e.g., *_R1.fastq.gz, *_R2.fastq.gz)
      --input_pattern      Pattern for FASTQ files [default: '*_R{1,2}.fastq.gz']
      
      --reference_db       Mothur-formatted reference database (FASTA)
      --taxonomy_db        Taxonomy file for classification
      --alignment_ref      Reference alignment for Mothur
      --reference_mafrax   Reference alignmente for Mafft and RaxML      
      --phylo_tree         Reference phylogenetic tree for RAxML-EPA
      
    
    Optional Arguments:
      --outdir             Output directory [default: ./results]
      --min_length         Minimum sequence length [default: 300]
      --max_length         Maximum sequence length [default: 580]
      --max_ambig          Maximum ambiguous bases [default: 0]
      --max_homop          Maximum homopolymer length [default: 8]
      --classify_cutoff    Classification cutoff percentage [default: 90]
      --align_start        Alignment start position [default: 4]
      --align_end          Alignment end position [default: 561]
      --cluster_cutoff     Clustering cutoff [default: 0.02]
      --cluster_method     Clustering method [default: 'dgc']
      --gappa_lwr_threshold GAPPA likelihood weight ratio threshold [default: 0.9]
      --skip_chimera_check Skip chimera detection [default: false]
      --skip_phylogenetic  Skip phylogenetic placement [default: false]
    
    """.stripIndent()
}

if (params.help) {
    helpMessage()
    exit 0
}

// Validate required parameters
if (!params.input && !params.input_dir) {
    log.error "ERROR: Must provide either --input (CSV) or --input_dir (directory)!"
    helpMessage()
    exit 1
}

if (params.input && params.input_dir) {
    log.error "ERROR: Cannot use both --input and --input_dir. Choose one!"
    helpMessage()
    exit 1
}

if (!params.reference_db || !params.taxonomy_db) {
    log.error "ERROR: --reference_db and --taxonomy_db are required!"
    helpMessage()
    exit 1
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { MOTHUR_MAKE_CONTIGS } from './modules/mothur_make_contigs'
include { MOTHUR_SCREEN_SEQS as MOTHUR_SCREEN_SEQS_1 } from './modules/mothur_screen_seqs'
include { MOTHUR_UNIQUE_SEQS as MOTHUR_UNIQUE_SEQS_1 } from './modules/mothur_unique_seqs'
include { MOTHUR_COUNT_SEQS } from './modules/mothur_count_seqs'
//include { MOTHUR_COUNT_SEQS as MOTHUR_COUNT_SEQS_2 } from './modules/mothur_count_seqs'
include { MOTHUR_CLASSIFY_SEQS } from './modules/mothur_classify_seqs'
include { MOTHUR_GET_LINEAGE } from './modules/mothur_get_lineage'
include { MOTHUR_ALIGN_SEQS } from './modules/mothur_align_seqs'
include { MOTHUR_SCREEN_SEQS_2 } from './modules/mothur_screen_seqs_2'
include { MOTHUR_FILTER_SEQS } from './modules/mothur_filter_seqs'
include { MOTHUR_UNIQUE_SEQS as MOTHUR_UNIQUE_SEQS_2 } from './modules/mothur_unique_seqs'
include { MOTHUR_PRE_CLUSTER } from './modules/mothur_pre_cluster'
include { MOTHUR_SPLIT_ABUND } from './modules/mothur_split_abund'
include { MOTHUR_CHIMERA_VSEARCH } from './modules/mothur_chimera_vsearch'
include { MOTHUR_REMOVE_SEQS } from './modules/mothur_remove_seqs'
include { MOTHUR_CLUSTER } from './modules/mothur_cluster'
include { MOTHUR_MAKE_SHARED } from './modules/mothur_make_shared'
include { MOTHUR_GET_OTUREP } from './modules/mothur_get_oturep'
include { MAFFT_ALIGN } from './modules/mafft_align'
include { RAXML_EPA } from './modules/raxml_epa'
include { GAPPA_ASSIGN } from './modules/gappa_assign'
include { GAPPA_MERGE } from './modules/gappa_merge'
include { GAPPA_TABLE } from './modules/gappa_table.nf'
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {
    
    //
    // STAGE 1: Create input channel (supports CSV or directory)
    //
    def input_ch
    
    if (params.input_dir) {
        // Option 1: Read from directory with pattern matching
        input_ch = channel
            .fromFilePairs("${params.input_dir}/${params.input_pattern}", checkIfExists: true)
            .map { sample_id, reads ->
                def meta = [id: sample_id]
                return tuple(meta, reads[0], reads[1])
            }
    } else if (params.input) {
        // Option 2: Read from CSV samplesheet
        input_ch = channel
            .fromPath(params.input, checkIfExists: true)
            .splitCsv(header: true)
            .map { row ->
                def meta = [id: row.sample]
                def fastq_1 = file(row.fastq_1, checkIfExists: true)
                def fastq_2 = file(row.fastq_2, checkIfExists: true)
                return tuple(meta, fastq_1, fastq_2)
            }
    } else {
        error "ERROR: Must provide either --input (CSV) or --input_dir (directory)"
    }
    
    //
    // STAGE 2: Mothur processing pipeline
    //
    
    // Make contigs from paired-end reads
    MOTHUR_MAKE_CONTIGS(input_ch)
    
    // Screen sequences (length, ambiguous bases, homopolymers)
    MOTHUR_SCREEN_SEQS_1(
        MOTHUR_MAKE_CONTIGS.out.contigs,
        params.min_length,
        params.max_length,
        params.max_ambig,
        params.max_homop
    )
    
    // Dereplicate sequences
    MOTHUR_UNIQUE_SEQS_1(MOTHUR_SCREEN_SEQS_1.out.seqs)
    
    // Count sequence occurrences
    MOTHUR_COUNT_SEQS(
        MOTHUR_UNIQUE_SEQS_1.out.fasta,
        MOTHUR_UNIQUE_SEQS_1.out.names
    )
    
    // Classify sequences
    classify_input = MOTHUR_UNIQUE_SEQS_1.out.fasta.join(MOTHUR_COUNT_SEQS.out.count, by: 0)
    // Resultado: [meta, fasta, count] - P1 com P1, P3 com P3

    MOTHUR_CLASSIFY_SEQS(
        classify_input,
        file(params.reference_db),
        file(params.taxonomy_db),
        params.classify_cutoff
    )
   
    // Join channels by meta.id to ensure all files are from the same sample
    lineage_input = MOTHUR_UNIQUE_SEQS_1.out.fasta
        .join(MOTHUR_COUNT_SEQS.out.count, by: 0)
        .join(MOTHUR_CLASSIFY_SEQS.out.taxonomy, by: 0)
 
    // Get Glomeromycota lineage
    MOTHUR_GET_LINEAGE(
    lineage_input.map { meta, fasta, count, taxonomy -> [meta, fasta] },
    lineage_input.map { meta, fasta, count, taxonomy -> [meta, count] },
    lineage_input.map { meta, fasta, count, taxonomy -> [meta, taxonomy] }
    )
    
    // Align to reference
    MOTHUR_ALIGN_SEQS(
        MOTHUR_GET_LINEAGE.out.fasta,
        file(params.alignment_ref)
    )
       
    //join the count table and aligned fasta
    ch_aligned_sync = MOTHUR_ALIGN_SEQS.out.aligned.join(MOTHUR_GET_LINEAGE.out.count, by: 0)
    // Screen aligned sequences
    
    MOTHUR_SCREEN_SEQS_2(
        ch_aligned_sync,
        params.align_start,
        params.align_end
    )
        
    // Filter alignment
    MOTHUR_FILTER_SEQS(MOTHUR_SCREEN_SEQS_2.out.seqs)
    
    // Second dereplication
    MOTHUR_UNIQUE_SEQS_2(MOTHUR_FILTER_SEQS.out.filtered)

    // The Pre-Cluster Detox
    precluster_sync = MOTHUR_UNIQUE_SEQS_2.out.fasta.join(MOTHUR_UNIQUE_SEQS_2.out.names, by: 0)
    MOTHUR_PRE_CLUSTER(precluster_sync)

    //The Split Abund Cleanse (Removing Singletons)
    split_abund_sync = MOTHUR_PRE_CLUSTER.out.fasta.join(MOTHUR_PRE_CLUSTER.out.count, by: 0)
    MOTHUR_SPLIT_ABUND(split_abund_sync)

    // Chimera detection (optional)
        chimera_sync = MOTHUR_SPLIT_ABUND.out.fasta.join(MOTHUR_SPLIT_ABUND.out.count, by: 0)

    if (!params.skip_chimera_check) {
        
        MOTHUR_CHIMERA_VSEARCH(chimera_sync)

        remove_sync = chimera_sync.join(MOTHUR_CHIMERA_VSEARCH.out.accnos, by: 0)

        MOTHUR_REMOVE_SEQS(
            remove_sync.map { meta, fasta, count, accnos -> tuple(meta, fasta, accnos, count) }
        )

        cluster_input_fasta = MOTHUR_REMOVE_SEQS.out.fasta
        cluster_input_count = MOTHUR_REMOVE_SEQS.out.count
        
    } else {
        // If skipping chimera, grab straight from SPLIT_ABUND
        cluster_input_fasta = MOTHUR_SPLIT_ABUND.out.fasta
        cluster_input_count = MOTHUR_SPLIT_ABUND.out.count
    }    
    // Cluster sequences
    MOTHUR_CLUSTER(
        cluster_input_fasta.join(cluster_input_count, by: 0),
        params.cluster_cutoff,
        params.cluster_method
    )

    // Make OTU table
    shared_sync = MOTHUR_CLUSTER.out.list.join(cluster_input_count, by: 0)

    MOTHUR_MAKE_SHARED(shared_sync)

    // Get representative sequences
    oturep_sync = MOTHUR_CLUSTER.out.list
        .join(cluster_input_fasta, by: 0)
        .join(cluster_input_count, by: 0)

    MOTHUR_GET_OTUREP(oturep_sync)

    //
    // STAGE 3: Post-Mothur phylogenetic analysis (optional)
    //
    if (!params.skip_phylogenetic) {
        
        // Align to reference with MAFFT
        MAFFT_ALIGN(
            MOTHUR_GET_OTUREP.out.rep_seqs,
            file(params.reference_mafrax)
        )
        
        // Phylogenetic placement with RAxML-EPA
        RAXML_EPA(
            MAFFT_ALIGN.out.alignment,
            file(params.phylo_tree),
            file(params.reference_mafrax)
        )
        
        // Taxonomic assignment with GAPPA
        GAPPA_ASSIGN(
            RAXML_EPA.out.jplace,
            file(params.taxonomy_ALMA),
            params.gappa_lwr_threshold
        )
	// Coleta todas as árvores que saíram do RAxML
        all_jplaces = RAXML_EPA.out.jplace.map { meta, jplace -> jplace }.collect()

        // Merge all individual jplace of each sample into one jplace
        GAPPA_MERGE(all_jplaces)

        // Table abundance count for ecological analysis
        
        all_per_query = GAPPA_ASSIGN.out.per_query.map { meta, file -> file }.collect()
        all_shared = MOTHUR_MAKE_SHARED.out.shared.map { meta, file -> file }.collect()
        GAPPA_TABLE(all_per_query, all_shared)


    }
    
    //
    // STAGE 4: Workflow completion
    //
    
    // Completion handlers will be added in workflow block below
}

workflow.onComplete {
    def duration_mins = workflow.duration.toMinutes()
    def status = workflow.success ? "SUCCESS ✅" : "FAILED ❌"
    
    log.info """
    =========================================
     Pipeline Execution Summary
    =========================================
    Status:    ${status}
    Duration:  ${duration_mins} minutes
    Output:    ${params.outdir}
    =========================================
    """.stripIndent()
}

workflow.onError {
    log.error """
    =========================================
     Pipeline Execution Error
    =========================================
    Error:     ${workflow.errorMessage}
    =========================================
    """.stripIndent()
}
