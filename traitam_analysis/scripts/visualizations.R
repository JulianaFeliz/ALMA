#!/usr/bin/env Rscript

# =============================================================================
# TraitAM Visualization Module
# =============================================================================
# Description: Comprehensive visualization suite for functional trait analysis
# Author: Seqera AI
# Date: 2026-04-17
# Version: 1.0.0
# =============================================================================

suppressPackageStartupMessages({
    library(tidyverse)
    library(ggpubr)
    library(corrplot)
    library(factoextra)
    library(viridis)
    library(patchwork)
})

# =============================================================================
# CONFIGURATION
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)

sample_id <- ifelse(length(args) >= 1, args[1], "sample1")
results_dir <- ifelse(length(args) >= 2, args[2], "results")

cat("\n=============================================================================\n")
cat("TraitAM Visualization Suite\n")
cat("=============================================================================\n\n")
cat("Sample:", sample_id, "\n")
cat("Results dir:", results_dir, "\n\n")

# Output directory
viz_dir <- file.path(results_dir, "06_visualizations")
dir.create(viz_dir, recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# LOAD DATA
# =============================================================================

cat("📂 Loading data...\n")

# Load traits with abundance
traits_file <- file.path(results_dir, "02_trait_data", 
                         paste0(sample_id, "_traits_abundance.tsv"))
traits_data <- read_tsv(traits_file, show_col_types = FALSE)

# Load FD metrics
fd_file <- file.path(results_dir, "03_functional_diversity", 
                     paste0(sample_id, "_FD_metrics.tsv"))
fd_metrics <- read_tsv(fd_file, show_col_types = FALSE)

# Load CWM
cwm_file <- file.path(results_dir, "03_functional_diversity", 
                      paste0(sample_id, "_CWM.tsv"))
cwm_data <- read_tsv(cwm_file, show_col_types = FALSE)

cat("  ✓ All data loaded\n\n")

# =============================================================================
# 1. TRAIT DISTRIBUTION PLOTS
# =============================================================================

cat("📊 Creating trait distribution plots...\n")

# Prepare data for plotting
trait_long <- traits_data %>%
    select(species_name, Family, Spore_vol_um3, Aspect_ratio, Color_score, 
           Wall_investment, Ornamentation_height_um, abundance) %>%
    distinct(species_name, .keep_all = TRUE) %>%
    pivot_longer(
        cols = c(Spore_vol_um3, Aspect_ratio, Color_score, 
                Wall_investment, Ornamentation_height_um),
        names_to = "Trait",
        values_to = "Value"
    )

# Create faceted density plots
p_distributions <- ggplot(trait_long, aes(x = Value, fill = Family)) +
    geom_density(alpha = 0.6) +
    facet_wrap(~Trait, scales = "free", ncol = 2) +
    scale_fill_viridis_d() +
    theme_bw() +
    labs(
        title = paste("Trait Distributions -", sample_id),
        subtitle = "Colored by taxonomic family",
        x = "Trait Value",
        y = "Density"
    ) +
    theme(
        legend.position = "bottom",
        strip.background = element_rect(fill = "lightblue")
    )

ggsave(
    file.path(viz_dir, paste0(sample_id, "_trait_distribution.pdf")),
    p_distributions,
    width = 10, height = 8
)

cat("  ✓ Trait distribution plot saved\n")

# =============================================================================
# 2. PCA BIPLOT
# =============================================================================

cat("📊 Creating PCA biplot...\n")

# Prepare trait matrix for PCA
trait_matrix <- traits_data %>%
    select(species_name, Spore_vol_um3, Aspect_ratio, Color_score, 
           Wall_investment, Ornamentation_height_um) %>%
    distinct(species_name, .keep_all = TRUE) %>%
    column_to_rownames("species_name") %>%
    na.omit()

# Perform PCA
pca_result <- prcomp(trait_matrix, scale. = TRUE)

# Create biplot
p_pca <- fviz_pca_biplot(
    pca_result,
    geom.ind = "point",
    pointshape = 21,
    pointsize = 3,
    fill.ind = "steelblue",
    col.ind = "black",
    col.var = "red",
    repel = TRUE,
    title = paste("PCA Biplot of Functional Traits -", sample_id),
    subtitle = paste("Explained variance: PC1 =", 
                    round(summary(pca_result)$importance[2,1]*100, 1), 
                    "%, PC2 =",
                    round(summary(pca_result)$importance[2,2]*100, 1), "%")
) +
    theme_bw()

ggsave(
    file.path(viz_dir, paste0(sample_id, "_PCA_biplot.pdf")),
    p_pca,
    width = 10, height = 8
)

cat("  ✓ PCA biplot saved\n")

# =============================================================================
# 3. TRAIT CORRELATION MATRIX
# =============================================================================

cat("📊 Creating trait correlation plot...\n")

# Calculate correlations
cor_matrix <- cor(trait_matrix, use = "complete.obs")

# Save as plot
pdf(file.path(viz_dir, paste0(sample_id, "_trait_correlations.pdf")), 
    width = 8, height = 8)

corrplot(
    cor_matrix, 
    method = "circle",
    type = "upper",
    order = "hclust",
    tl.col = "black",
    tl.srt = 45,
    addCoef.col = "black",
    number.cex = 0.7,
    title = paste("Trait Correlations -", sample_id),
    mar = c(0, 0, 2, 0)
)

dev.off()

cat("  ✓ Trait correlation plot saved\n")

# =============================================================================
# 4. ABUNDANCE-TRAIT RELATIONSHIPS
# =============================================================================

cat("📊 Creating abundance-trait relationship plots...\n")

# Create scatterplots for each trait vs abundance
trait_names <- c("Spore_vol_um3", "Aspect_ratio", "Color_score", 
                 "Wall_investment", "Ornamentation_height_um")

abundance_plots <- list()

for (trait_name in trait_names) {
    # Log-transform abundance for better visualization
    plot_data <- traits_data %>%
        distinct(species_name, .keep_all = TRUE) %>%
        mutate(log_abundance = log10(abundance + 1))
    
    p <- ggplot(plot_data, aes_string(x = trait_name, y = "log_abundance")) +
        geom_point(aes(color = Family), size = 3, alpha = 0.7) +
        geom_smooth(method = "lm", se = TRUE, color = "black", linetype = "dashed") +
        scale_color_viridis_d() +
        theme_bw() +
        labs(
            x = str_replace_all(trait_name, "_", " "),
            y = "log10(Abundance + 1)",
            title = str_replace_all(trait_name, "_", " ")
        ) +
        theme(legend.position = "none")
    
    abundance_plots[[trait_name]] <- p
}

# Combine plots
p_abundance_combined <- wrap_plots(abundance_plots, ncol = 2) +
    plot_annotation(
        title = paste("Abundance-Trait Relationships -", sample_id),
        subtitle = "Do abundant species have specific trait values?"
    )

ggsave(
    file.path(viz_dir, paste0(sample_id, "_abundance_traits.pdf")),
    p_abundance_combined,
    width = 12, height = 10
)

cat("  ✓ Abundance-trait plots saved\n")

# =============================================================================
# 5. FUNCTIONAL DIVERSITY BARPLOT
# =============================================================================

cat("📊 Creating functional diversity barplot...\n")

# Reshape FD metrics for plotting
fd_long <- fd_metrics %>%
    select(Sample, FRic, FEve, FDiv, FDis, RaoQ) %>%
    pivot_longer(-Sample, names_to = "Metric", values_to = "Value")

p_fd_barplot <- ggplot(fd_long, aes(x = Metric, y = Value, fill = Metric)) +
    geom_col(width = 0.7) +
    geom_text(aes(label = round(Value, 3)), vjust = -0.5) +
    scale_fill_viridis_d() +
    theme_bw() +
    labs(
        title = paste("Functional Diversity Metrics -", sample_id),
        x = "Functional Diversity Index",
        y = "Value"
    ) +
    theme(
        legend.position = "none",
        axis.text.x = element_text(angle = 0, hjust = 0.5)
    ) +
    ylim(0, max(fd_long$Value) * 1.2)

ggsave(
    file.path(viz_dir, paste0(sample_id, "_functional_barplot.pdf")),
    p_fd_barplot,
    width = 8, height = 6
)

cat("  ✓ Functional diversity barplot saved\n")

# =============================================================================
# 6. COMMUNITY WEIGHTED MEANS (CWM) PLOT
# =============================================================================

cat("📊 Creating community weighted means plot...\n")

# Reshape CWM data
cwm_long <- cwm_data %>%
    pivot_longer(-Sample, names_to = "Trait", values_to = "CWM")

p_cwm <- ggplot(cwm_long, aes(x = Trait, y = CWM, fill = Trait)) +
    geom_col(width = 0.7) +
    geom_text(aes(label = round(CWM, 2)), vjust = -0.5, size = 3) +
    scale_fill_viridis_d() +
    theme_bw() +
    labs(
        title = paste("Community Weighted Mean Trait Values -", sample_id),
        subtitle = "Abundance-weighted average trait values",
        x = "Trait",
        y = "Community Weighted Mean"
    ) +
    theme(
        legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1)
    ) +
    ylim(0, max(cwm_long$CWM) * 1.2)

ggsave(
    file.path(viz_dir, paste0(sample_id, "_CWM_barplot.pdf")),
    p_cwm,
    width = 10, height = 6
)

cat("  ✓ CWM plot saved\n")

# =============================================================================
# 7. TRAIT-TRAIT SCATTERPLOT MATRIX
# =============================================================================

cat("📊 Creating trait scatterplot matrix...\n")

# Create pairs plot with taxonomic family colors
plot_data <- traits_data %>%
    distinct(species_name, .keep_all = TRUE) %>%
    select(species_name, Family, Spore_vol_um3, Aspect_ratio, Color_score, 
           Wall_investment, Ornamentation_height_um)

# Custom pairs plot function
pdf(file.path(viz_dir, paste0(sample_id, "_trait_pairs.pdf")), 
    width = 12, height = 12)

pairs(
    plot_data[, c("Spore_vol_um3", "Aspect_ratio", "Color_score", 
                  "Wall_investment", "Ornamentation_height_um")],
    col = as.factor(plot_data$Family),
    pch = 19,
    main = paste("Trait Scatterplot Matrix -", sample_id),
    cex = 1.5
)

dev.off()

cat("  ✓ Trait pairs plot saved\n")

# =============================================================================
# 8. SUMMARY FIGURE
# =============================================================================

cat("📊 Creating summary figure...\n")

# Combine key plots into one summary figure
p_summary <- (p_distributions + p_pca) / (p_fd_barplot + p_cwm) +
    plot_annotation(
        title = paste("TraitAM Functional Analysis Summary -", sample_id),
        subtitle = paste("Generated:", Sys.Date())
    )

ggsave(
    file.path(viz_dir, paste0(sample_id, "_summary_figure.pdf")),
    p_summary,
    width = 16, height = 12
)

cat("  ✓ Summary figure saved\n\n")

# =============================================================================
# SUMMARY
# =============================================================================

cat("=============================================================================\n")
cat("✅ Visualization Complete\n")
cat("=============================================================================\n\n")

cat("Plots created:\n")
cat("  1. Trait distributions (by family)\n")
cat("  2. PCA biplot\n")
cat("  3. Trait correlation matrix\n")
cat("  4. Abundance-trait relationships\n")
cat("  5. Functional diversity barplot\n")
cat("  6. Community weighted means\n")
cat("  7. Trait scatterplot matrix\n")
cat("  8. Summary figure (combined)\n\n")

cat("All visualizations saved to:", viz_dir, "\n\n")

cat("Next step:\n")
cat("  Generate HTML report: Rscript generate_report.R", sample_id, results_dir, "\n\n")
