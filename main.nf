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
params.reference_db       = null  // ALMAdb (Curated Glomeromycota Database)
params.sintax_db          = null  // Eukaryome SINTAX Database (.udb)
params.reference_mafrax   = null  // Reference alignmente for mafft and Raxml (without problematic sequences)
params.phylo_tree         = null  // Reference phylogenetic tree
params.taxonomy_ALMA      = null  // Taxonomy dictionary for GAPPA

params.min_length         = 200
params.max_length         = 600
params.max_ambig          = 0     
params.cluster_cutoff     = 0.03  // 97% similarity for OTU clustering

// Post-Clustering parameters
params.gappa_lwr_threshold = 0.9

// Workflow control
params.skip_phylogenetic  = false

// Help message
params.help = false

def helpMessage() {
    log.info"""
    =========================================
      ALMA - AMF LSU METABARCODING ANALYSIS
    =========================================
    
    Usage (Option 1 - CSV):
      nextflow run main.nf --input samplesheet.csv --reference_db ALMAdb.fasta --sintax_db Eukaryome.udb
    
    Usage (Option 2 - Directory):
      nextflow run main.nf --input_dir data/ --reference_db ALMAdb.fasta --sintax_db Eukaryome.udb
    
    Required Arguments (choose ONE input method):
      --input              Path to samplesheet CSV with columns: sample,fastq_1,fastq_2
                           OR
      --input_dir          Directory containing paired FASTQ files (e.g., *_R1.fastq.gz, *_R2.fastq.gz)
      --input_pattern      Pattern for FASTQ files [default: '*_R{1,2}.fastq.gz']
      
      --reference_db       ALMAdb formatted reference database (FASTA)
      --sintax_db          Eukaryome SINTAX database (.udb format)
      --reference_mafrax   Reference alignment for MAFFT and RAxML      
      --phylo_tree         Reference phylogenetic tree for RAxML-EPA
      --taxonomy_ALMA      Taxonomy dictionary for GAPPA
      
    Optional Arguments:
      --outdir             Output directory [default: ./results]
      --min_length         Minimum sequence length [default: 200]
      --max_length         Maximum sequence length [default: 600]
      --max_ambig          Maximum ambiguous bases allowed [default: 0]
      --cluster_cutoff     Clustering divergence cutoff (0.03 = 97%) [default: 0.03]
      --gappa_lwr_threshold GAPPA likelihood weight ratio threshold [default: 0.9]
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

if (!params.reference_db || !params.sintax_db) {
    log.error "ERROR: --reference_db (ALMAdb) and --sintax_db (Eukaryome) are required!"
    helpMessage()
    exit 1
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// VSEARCH Modules (Phase 1 & 2)
include { VSEARCH_MERGE }          from './modules/vsearch_merge.nf'
include { VSEARCH_FILTER }         from './modules/vsearch_filter.nf'
include { VSEARCH_DEREP }          from './modules/vsearch_derep.nf'
include { VSEARCH_CHIMERA }        from './modules/vsearch_chimera.nf'
include { VSEARCH_SINTAX_FILTER }  from './modules/vsearch_sintax_filter.nf'
include { VSEARCH_GLOBAL_DEREP }   from './modules/vsearch_global_derep.nf'
include { VSEARCH_CLOSED_REF_CLUSTER }     from './modules/vsearch_closed_ref_cluster.nf'
include { VSEARCH_DENOVO_CLUSTER } from './modules/vsearch_denovo_cluster.nf'
include { VSEARCH_MAP_OTUS }       from './modules/vsearch_map_otu.nf'

// Phylogenetic Modules (Phase 3)
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
    // STAGE 2: VSEARCH Processing Pipeline (Quality & Taxonomic Filtering)
    //
    
    VSEARCH_MERGE(input_ch)
    
    VSEARCH_FILTER(
        VSEARCH_MERGE.out.fastq,
        params.min_length,
        params.max_length,
        params.max_ambig
    )
    
    VSEARCH_DEREP(VSEARCH_FILTER.out.fasta)
    
    VSEARCH_CHIMERA(VSEARCH_DEREP.out.fasta)
    
    // The Bouncer: Filter out plants, bacteria, and non-AMF fungi using Eukaryome
    VSEARCH_SINTAX_FILTER(
        VSEARCH_CHIMERA.out.fasta,
        file(params.sintax_db)
    )

    // Collect all valid Glomeromycota FASTAs into a single channel for global clustering.
    // The .map is crucial here to extract just the files and leave the meta IDs behind.
    ch_clean_glomeros = VSEARCH_SINTAX_FILTER.out.fasta
        .map { meta, fasta -> fasta }
        .collect()

    //
    // STAGE 3: Open Reference Clustering (The Jeweler)
    //
    
    VSEARCH_GLOBAL_DEREP(ch_clean_glomeros)

    VSEARCH_CLOSED_REF(
        VSEARCH_GLOBAL_DEREP.out.global_derep,
        file(params.reference_db),
        params.cluster_cutoff
    )

    VSEARCH_DENOVO_CLUSTER(
        VSEARCH_CLOSED_REF.out.unmapped,
        params.cluster_cutoff
    )

    VSEARCH_MAP_OTUS(
        VSEARCH_CLOSED_REF.out.centroids,
        VSEARCH_DENOVO_CLUSTER.out.centroids,
        VSEARCH_GLOBAL_DEREP.out.all_samples,
        params.cluster_cutoff
    )

    //
    // STAGE 4: Post-Clustering Phylogenetic Analysis
    //
    if (!params.skip_phylogenetic) {
        
        // Align to reference with MAFFT
        MAFFT_ALIGN(
            VSEARCH_MAP_OTUS.out.rep_seqs,
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

        GAPPA_TABLE(all_per_query, VSEARCH_MAP_OTUS.out.shared)


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
