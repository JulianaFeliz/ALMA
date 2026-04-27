#!/usr/bin/env Rscript

# =============================================================================
# TraitAM Report Generation Script
# =============================================================================
# Description: Compile RMarkdown report from analysis results
# Author: Seqera AI
# Date: 2026-04-17
# Version: 1.0.0
# =============================================================================

suppressPackageStartupMessages({
    library(rmarkdown)
    library(knitr)
})

# =============================================================================
# CONFIGURATION
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
    cat("
USAGE: Rscript generate_report.R <sample_id> <results_dir> [traitam_dir]

ARGUMENTS:
  sample_id     Sample identifier
  results_dir   Path to results directory
  traitam_dir   Path to TraitAM data directory (default: data/)

EXAMPLE:
  Rscript generate_report.R sample1 results/ data/

")
    quit(status = 1)
}

sample_id <- args[1]
results_dir <- args[2]
traitam_dir <- ifelse(length(args) >= 3, args[3], "data")

cat("\n=============================================================================\n")
cat("TraitAM Report Generation\n")
cat("=============================================================================\n\n")
cat("Sample:", sample_id, "\n")
cat("Results:", results_dir, "\n")
cat("TraitAM:", traitam_dir, "\n\n")

# =============================================================================
# VALIDATE FILES
# =============================================================================

cat("📋 Validating required files...\n")

# Check for Rmd template
rmd_template <- "scripts/report_template.Rmd"

if (!file.exists(rmd_template)) {
    stop("ERROR: Report template not found: ", rmd_template)
}

# Check for key result files
required_files <- c(
    file.path(results_dir, "01_matched_taxa", paste0(sample_id, "_matched_species.tsv")),
    file.path(results_dir, "02_trait_data", paste0(sample_id, "_traits_abundance.tsv")),
    file.path(results_dir, "03_functional_diversity", paste0(sample_id, "_FD_metrics.tsv"))
)

missing_files <- c()
for (f in required_files) {
    if (!file.exists(f)) {
        missing_files <- c(missing_files, f)
    }
}

if (length(missing_files) > 0) {
    stop("ERROR: Missing required result files:\n  ",
         paste(missing_files, collapse = "\n  "),
         "\n\nRun analysis first: Rscript run_traitam_analysis.R")
}

cat("  ✓ All required files present\n\n")

# =============================================================================
# RENDER REPORT
# =============================================================================

cat("📝 Rendering HTML report...\n")

output_file <- file.path(
    results_dir, 
    "07_reports", 
    paste0(sample_id, "_analysis_report.html")
)

# Create output directory
dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)

# Render the report
tryCatch({
    render(
        input = rmd_template,
        output_file = basename(output_file),
        output_dir = dirname(output_file),
        params = list(
            sample_id = sample_id,
            results_dir = results_dir,
            traitam_dir = traitam_dir
        ),
        quiet = FALSE
    )
    
    cat("\n✅ Report successfully generated!\n\n")
    cat("📄 Report location:", output_file, "\n")
    cat("📊 File size:", round(file.info(output_file)$size / 1024 / 1024, 2), "MB\n\n")
    
}, error = function(e) {
    cat("\n❌ ERROR generating report:\n")
    cat(conditionMessage(e), "\n\n")
    quit(status = 1)
})

# =============================================================================
# CREATE TEXT SUMMARY
# =============================================================================

cat("📋 Creating text summary...\n")

# Load key results for summary
library(tidyverse, quietly = TRUE, warn.conflicts = FALSE)

matched_taxa <- read_tsv(
    file.path(results_dir, "01_matched_taxa", paste0(sample_id, "_matched_species.tsv")),
    show_col_types = FALSE
)

fd_metrics <- read_tsv(
    file.path(results_dir, "03_functional_diversity", paste0(sample_id, "_FD_metrics.tsv")),
    show_col_types = FALSE
)

# Create summary text
summary_text <- paste0(
    "TraitAM Functional Analysis Summary\n",
    "===================================\n\n",
    "Sample: ", sample_id, "\n",
    "Date: ", Sys.Date(), "\n",
    "Pipeline Version: 1.0.0\n\n",
    "TAXONOMIC SUMMARY\n",
    "-----------------\n",
    "Total OTUs analyzed: ", nrow(matched_taxa), "\n",
    "Species matched to TraitAM: ", sum(!is.na(matched_taxa$match_type)), "\n",
    "  - Exact species matches: ", sum(matched_taxa$match_type == "exact_species", na.rm = TRUE), "\n",
    "  - Fuzzy matches: ", sum(grepl("fuzzy", matched_taxa$match_type), na.rm = TRUE), "\n",
    "  - Genus-level matches: ", sum(matched_taxa$match_type == "genus_level", na.rm = TRUE), "\n",
    "Unmatched: ", sum(is.na(matched_taxa$match_type)), "\n\n",
    "FUNCTIONAL DIVERSITY METRICS\n",
    "----------------------------\n",
    "FRic (Functional Richness): ", round(fd_metrics$FRic, 4), "\n",
    "FEve (Functional Evenness): ", round(fd_metrics$FEve, 4), "\n",
    "FDiv (Functional Divergence): ", round(fd_metrics$FDiv, 4), "\n",
    "FDis (Functional Dispersion): ", round(fd_metrics$FDis, 4), "\n",
    "RaoQ (Rao's Quadratic Entropy): ", round(fd_metrics$RaoQ, 4), "\n\n",
    "PHYLOGENETIC ANALYSIS\n",
    "---------------------\n"
)

# Add phylogenetic metrics if available
phylo_div_file <- file.path(results_dir, "04_phylogenetic_analysis", 
                             paste0(sample_id, "_phylo_diversity.tsv"))
if (file.exists(phylo_div_file)) {
    phylo_div <- read_tsv(phylo_div_file, show_col_types = FALSE)
    summary_text <- paste0(
        summary_text,
        "Faith's PD: ", round(phylo_div$PD, 2), "\n",
        "Species Richness: ", phylo_div$SR, "\n",
        "MPD: ", round(phylo_div$MPD, 3), "\n",
        "MNTD: ", round(phylo_div$MNTD, 3), "\n\n"
    )
} else {
    summary_text <- paste0(summary_text, "Not calculated (insufficient species in phylogeny)\n\n")
}

summary_text <- paste0(
    summary_text,
    "OUTPUT FILES\n",
    "------------\n",
    "HTML Report: ", output_file, "\n",
    "Visualizations: ", file.path(results_dir, "06_visualizations"), "\n",
    "Trait Data: ", file.path(results_dir, "02_trait_data"), "\n",
    "FD Metrics: ", file.path(results_dir, "03_functional_diversity"), "\n\n",
    "CITATION\n",
    "--------\n",
    "TraitAM Database: Chaudhary et al. (2024)\n",
    "FD Package: Laliberté & Legendre (2010) Ecology\n",
    "picante: Kembel et al. (2010) Bioinformatics\n\n",
    "For full methods and references, see HTML report.\n"
)

# Save summary
summary_file <- file.path(results_dir, "07_reports", paste0(sample_id, "_summary.txt"))
writeLines(summary_text, summary_file)

cat("  ✓ Text summary saved:", summary_file, "\n\n")

# =============================================================================
# COMPLETION
# =============================================================================

cat("=============================================================================\n")
cat("✅ REPORT GENERATION COMPLETE\n")
cat("=============================================================================\n\n")

cat("📊 Generated files:\n")
cat("  1. HTML report:", output_file, "\n")
cat("  2. Text summary:", summary_file, "\n\n")

cat("🌐 To view the HTML report:\n")
cat("  Open in browser: file://", normalizePath(output_file), "\n\n")

cat("Next steps:\n")
cat("  • Review HTML report for complete analysis\n")
cat("  • Examine visualizations in 06_visualizations/\n")
cat("  • Launch interactive dashboard: Rscript launch_dashboard.R\n\n")

cat("🎉 Analysis complete!\n\n")
