#!/usr/bin/env Rscript

# =============================================================================
# TraitAM Phylogenetic Analysis Module
# =============================================================================
# Description: Phylogenetic signal, diversity, and tree comparison
# Author: Seqera AI
# Date: 2026-04-17
# Version: 1.0.0
# =============================================================================

suppressPackageStartupMessages({
    library(tidyverse)
    library(ape)
    library(phytools)
    library(picante)
    library(ggtree)
    library(ggpubr)
})

# =============================================================================
# CONFIGURATION
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)

# Default paths (can be overridden with command line args)
sample_id <- ifelse(length(args) >= 1, args[1], "sample1")
results_dir <- ifelse(length(args) >= 2, args[2], "results")
traitam_dir <- ifelse(length(args) >= 3, args[3], "data")

cat("\n=============================================================================\n")
cat("TraitAM Phylogenetic Analysis\n")
cat("=============================================================================\n\n")
cat("Sample:", sample_id, "\n")
cat("Results dir:", results_dir, "\n")
cat("TraitAM dir:", traitam_dir, "\n\n")

# =============================================================================
# LOAD DATA
# =============================================================================

cat("📂 Loading data...\n")

# Load trait data with abundance
traits_file <- file.path(results_dir, "02_trait_data", 
                         paste0(sample_id, "_traits_abundance.tsv"))

if (!file.exists(traits_file)) {
    stop("Trait file not found. Run run_traitam_analysis.R first.")
}

traits_data <- read_tsv(traits_file, show_col_types = FALSE)

# Load TraitAM phylogenetic tree
tree_file <- file.path(traitam_dir, "DataRecord_7_FinalTree_Oct2023.tre")
traitam_tree <- read.tree(tree_file)

cat("  - Traits loaded:", nrow(traits_data), "records\n")
cat("  - TraitAM tree loaded:", Ntip(traitam_tree), "tips\n\n")

# =============================================================================
# PHYLOGENETIC SIGNAL ANALYSIS
# =============================================================================

cat("🌳 Analyzing phylogenetic signal of traits...\n\n")

# Prepare data: unique species with trait values
trait_values <- traits_data %>%
    group_by(species_name) %>%
    summarise(
        Spore_vol_um3 = first(Spore_vol_um3),
        Aspect_ratio = first(Aspect_ratio),
        Color_score = first(Color_score),
        Wall_investment = first(Wall_investment),
        Ornamentation_height_um = first(Ornamentation_height_um),
        .groups = "drop"
    )

# Match species to tree
# Tree tips might have different format (e.g., Genus_species)
# Need to match carefully

tree_tips <- traitam_tree$tip.label

# Try to match species names to tree tips
# Common formats: "Genus_species", "Genus species", accession numbers
match_to_tree <- function(species_name, tree_tips) {
    # Try exact match
    if (species_name %in% tree_tips) {
        return(species_name)
    }
    
    # Try with underscore
    species_underscore <- str_replace(species_name, " ", "_")
    if (species_underscore %in% tree_tips) {
        return(species_underscore)
    }
    
    # Try partial match (genus + species)
    genus_species <- str_extract(species_name, "^\\S+ \\S+")
    matches <- grep(genus_species, tree_tips, value = TRUE, fixed = TRUE)
    if (length(matches) > 0) {
        return(matches[1])
    }
    
    return(NA_character_)
}

trait_values$tree_tip <- sapply(trait_values$species_name, match_to_tree, tree_tips)

# Filter to species in tree
trait_values_in_tree <- trait_values %>%
    filter(!is.na(tree_tip))

n_matched <- nrow(trait_values_in_tree)
n_total <- nrow(trait_values)

cat("Species matched to tree:", n_matched, "/", n_total, "\n\n")

if (n_matched < 3) {
    cat("⚠️  WARNING: Too few species matched to tree (n=", n_matched, ")\n")
    cat("   Skipping phylogenetic signal analysis\n\n")
    phylo_signal_results <- NULL
} else {
    # Prune tree to matched species
    tree_pruned <- keep.tip(traitam_tree, trait_values_in_tree$tree_tip)
    
    cat("Pruned tree:", Ntip(tree_pruned), "tips\n\n")
    
    # Calculate phylogenetic signal for each trait
    phylo_signal_results <- data.frame(
        Trait = character(),
        Blomberg_K = numeric(),
        K_pvalue = numeric(),
        Pagel_lambda = numeric(),
        Lambda_pvalue = numeric(),
        Interpretation = character(),
        stringsAsFactors = FALSE
    )
    
    traits_to_test <- c("Spore_vol_um3", "Aspect_ratio", "Color_score", 
                        "Wall_investment", "Ornamentation_height_um")
    
    for (trait_name in traits_to_test) {
        cat("Testing:", trait_name, "\n")
        
        # Prepare trait vector
        trait_vec <- trait_values_in_tree[[trait_name]]
        names(trait_vec) <- trait_values_in_tree$tree_tip
        
        # Remove NAs
        trait_vec <- trait_vec[!is.na(trait_vec)]
        
        if (length(trait_vec) < 3) {
            cat("  ⚠️  Skipping (too few values)\n\n")
            next
        }
        
        # Prune tree to available tips
        tree_for_trait <- keep.tip(tree_pruned, names(trait_vec))
        
        # Blomberg's K
        tryCatch({
            K_result <- phylosig(tree_for_trait, trait_vec, method = "K", test = TRUE)
            K_value <- K_result$K
            K_pval <- K_result$P
        }, error = function(e) {
            K_value <- NA
            K_pval <- NA
        })
        
        # Pagel's lambda
        tryCatch({
            lambda_result <- phylosig(tree_for_trait, trait_vec, method = "lambda", test = TRUE)
            lambda_value <- lambda_result$lambda
            lambda_pval <- lambda_result$P
        }, error = function(e) {
            lambda_value <- NA
            lambda_pval <- NA
        })
        
        # Interpretation
        interpretation <- if (!is.na(K_value)) {
            if (K_value < 0.5) {
                "Weak phylogenetic signal"
            } else if (K_value < 1) {
                "Moderate phylogenetic signal"
            } else {
                "Strong phylogenetic signal (conserved)"
            }
        } else {
            "Could not calculate"
        }
        
        cat("  K =", round(K_value, 3), ", p =", round(K_pval, 4), "\n")
        cat("  λ =", round(lambda_value, 3), ", p =", round(lambda_pval, 4), "\n")
        cat("  →", interpretation, "\n\n")
        
        # Add to results
        phylo_signal_results <- rbind(phylo_signal_results, data.frame(
            Trait = trait_name,
            Blomberg_K = K_value,
            K_pvalue = K_pval,
            Pagel_lambda = lambda_value,
            Lambda_pvalue = lambda_pval,
            Interpretation = interpretation
        ))
    }
    
    # Save results
    write_tsv(
        phylo_signal_results,
        file.path(results_dir, "04_phylogenetic_analysis", 
                  paste0(sample_id, "_phylo_signal.tsv"))
    )
    
    cat("✓ Phylogenetic signal results saved\n\n")
}

# =============================================================================
# PHYLOGENETIC DIVERSITY METRICS
# =============================================================================

cat("📊 Calculating phylogenetic diversity metrics...\n")

if (n_matched >= 3) {
    # Prepare community matrix (samples × species)
    comm_matrix <- traits_data %>%
        select(sample, species_name, abundance) %>%
        filter(species_name %in% trait_values_in_tree$species_name) %>%
        left_join(
            trait_values_in_tree %>% select(species_name, tree_tip),
            by = "species_name"
        ) %>%
        select(sample, tree_tip, abundance) %>%
        pivot_wider(names_from = tree_tip, values_from = abundance, values_fill = 0) %>%
        column_to_rownames("sample") %>%
        as.matrix()
    
    # Ensure tree and community matrix match
    tree_for_pd <- keep.tip(traitam_tree, colnames(comm_matrix))
    
    # Faith's Phylogenetic Diversity
    pd_result <- pd(comm_matrix, tree_for_pd)
    
    # Mean Pairwise Distance (MPD)
    mpd_result <- mpd(comm_matrix, cophenetic(tree_for_pd))
    
    # Mean Nearest Taxon Distance (MNTD)
    mntd_result <- mntd(comm_matrix, cophenetic(tree_for_pd))
    
    # Combine results
    phylo_diversity <- data.frame(
        Sample = sample_id,
        PD = pd_result$PD,
        SR = pd_result$SR,  # Species richness
        MPD = mpd_result,
        MNTD = mntd_result
    )
    
    cat("  ✓ Phylogenetic diversity calculated\n")
    cat("    - Faith's PD:", round(phylo_diversity$PD, 2), "\n")
    cat("    - Species richness:", phylo_diversity$SR, "\n")
    cat("    - MPD:", round(phylo_diversity$MPD, 3), "\n")
    cat("    - MNTD:", round(phylo_diversity$MNTD, 3), "\n\n")
    
    # Save results
    write_tsv(
        phylo_diversity,
        file.path(results_dir, "04_phylogenetic_analysis", 
                  paste0(sample_id, "_phylo_diversity.tsv"))
    )
} else {
    cat("  ⚠️  Skipping phylogenetic diversity (too few species)\n\n")
}

# =============================================================================
# TREE VISUALIZATION
# =============================================================================

cat("🎨 Creating phylogenetic tree visualization...\n")

if (n_matched >= 3) {
    # Highlight detected species on TraitAM tree
    tree_plot <- ggtree(tree_pruned, layout = "circular") +
        geom_tiplab(size = 2, offset = 0.01) +
        theme_tree2()
    
    # Add trait heatmap
    # Prepare trait data for ggtree
    trait_heatmap_data <- trait_values_in_tree %>%
        select(tree_tip, Spore_vol_um3, Aspect_ratio, Color_score, 
               Wall_investment, Ornamentation_height_um) %>%
        column_to_rownames("tree_tip")
    
    # Normalize traits for visualization (0-1 scale)
    trait_heatmap_norm <- apply(trait_heatmap_data, 2, function(x) {
        (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
    })
    
    # Create tree + heatmap plot
    p_tree <- gheatmap(tree_plot, trait_heatmap_norm, 
                       offset = 0.05, width = 0.3, 
                       colnames_angle = 90, hjust = 1) +
        scale_fill_gradient2(low = "blue", mid = "white", high = "red", 
                            midpoint = 0.5, name = "Trait value\n(normalized)") +
        labs(title = paste("Phylogenetic Tree with Trait Data -", sample_id),
             subtitle = paste(n_matched, "species detected in sample"))
    
    # Save plot
    ggsave(
        file.path(results_dir, "04_phylogenetic_analysis", 
                  paste0(sample_id, "_tree_comparison.pdf")),
        p_tree,
        width = 12, height = 12
    )
    
    cat("  ✓ Tree visualization saved\n\n")
} else {
    cat("  ⚠️  Skipping tree visualization (too few species)\n\n")
}

# =============================================================================
# SUMMARY
# =============================================================================

cat("=============================================================================\n")
cat("✅ Phylogenetic Analysis Complete\n")
cat("=============================================================================\n\n")

cat("Results saved to:", file.path(results_dir, "04_phylogenetic_analysis"), "\n\n")

if (!is.null(phylo_signal_results)) {
    cat("Key findings:\n")
    cat("  • Phylogenetic signal analysis:", nrow(phylo_signal_results), "traits tested\n")
    
    strong_signal <- sum(phylo_signal_results$Blomberg_K > 1, na.rm = TRUE)
    if (strong_signal > 0) {
        cat("  • Traits with strong phylogenetic signal:", strong_signal, "\n")
    }
    
    cat("  • Phylogenetic diversity (PD):", round(phylo_diversity$PD, 2), "\n")
    cat("  • Species in phylogeny:", n_matched, "\n\n")
}

cat("Next steps:\n")
cat("  1. Review phylogenetic signal results\n")
cat("  2. Examine tree comparison plot\n")
cat("  3. Run visualizations: Rscript visualizations.R\n")
cat("  4. Generate report: Rscript generate_report.R\n\n")
