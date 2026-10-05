# =============================================================================
# ALMA - Advanced Ecological Analysis (Article Strategies)
# Following Boshuizen et al. 2021 & Leite & Kuramae 2025
# =============================================================================

# Package verification (Bioconductor)
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
# BiocManager::install(c("microbiome", "zCompositions")) # Uncomment if you need to install

library(vegan)
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggrepel)
library(RColorBrewer)
library(patchwork)
library(microbiome)

# -----------------------------------------------------------------------------
# 1. UNIVERSAL DATA LOADING AND METADATA INTEGRATION
# -----------------------------------------------------------------------------
# 1.1 Load OTU abundance table
otu <- read.table("ALMA_OTU_abundance.tsv", header = TRUE, sep = "\t", 
                  row.names = 2, check.names = FALSE)
comm <- otu[, !colnames(otu) %in% c("label", "numOtus")]

# 1.2 Load Metadata table (User must provide this file)
# Expected format: Column 1 = SampleID, Column 2 = Group
metadata_raw <- read.table("metadata.tsv", header = TRUE, sep = "\t", row.names = 1)

# 1.3 Sync and Order Data
# Find samples that exist in BOTH the abundance table and the metadata
common_samples <- intersect(rownames(metadata_raw), rownames(comm))

if(length(common_samples) == 0) {
  stop("ERROR: Sample names in ALMA_OTU_abundance.tsv do not match metadata.tsv!")
}

# Subset and force the abundance table to follow the EXACT order of the metadata file
comm <- comm[common_samples, ]
metadata <- metadata_raw[common_samples, , drop = FALSE]

# Ensure the 'Group' column is explicitly named for the rest of the script
colnames(metadata)[1] <- "Group" 

# 1.4 Dynamic Factors and Colors
# Lock the Sample order for plotting based on the metadata row order
metadata$Sample <- factor(rownames(metadata), levels = rownames(metadata))
metadata$Group  <- as.factor(metadata$Group)

# Generate a dynamic color palette based on the number of groups
unique_groups <- levels(metadata$Group)
default_palette <- c("#E07B54", "#4A9E8A", "#3B82F6", "#F59E0B", "#8B5CF6", "#EC4899")

# Assign colors (if more than 6 groups, it repeats, but usually it's 2 to 4)
col_groups <- setNames(default_palette[1:length(unique_groups)], unique_groups)

cat("Check: Universal integration successful.\n")
cat("Samples loaded:", nrow(comm), "| Taxa:", ncol(comm), "\n")
cat("Groups detected:", paste(unique_groups, collapse = ", "), "\n")

# -----------------------------------------------------------------------------
# 2. EXPLORATION (DE) & SPARSITY 
# -----------------------------------------------------------------------------

#Rarefaction Curve 
rare_data <- rarecurve(comm, step = 50, tidy = TRUE)

# We merge the curve data with your metadata to fetch the Groups and Colors
rare_data_merged <- rare_data %>%
  left_join(metadata, by = c("Site" = "Sample"))

# Draw the plot using ggplot2
p_rar <- ggplot(rare_data_merged, aes(x = Sample, y = Species, group = Site, color = Group)) +
  geom_line(linewidth = 0.8, alpha = 0.8) +
  geom_vline(xintercept = min(rowSums(comm)), linetype = "dashed", color = "red") +
  scale_color_manual(values = col_groups) +
  theme_bw(base_size = 12) +
  labs(title = "Rarefaction Curves",
       subtitle = "Saturation Check",
       x = "Sequencing Depth (Reads)",
       y = "OTU Richness") +
  theme(legend.position = "none")


# Cleveland Dotplot
df_dot <- data.frame(Sample = metadata$Sample, Group = metadata$Group, Depth = rowSums(comm))
df_dot$Sample <- factor(df_dot$Sample, levels = rev(levels(df_dot$Sample)))

p_cleveland <- ggplot(df_dot, aes(x = Depth, y = Sample, color = Group)) +
  geom_segment(aes(yend = Sample, xend = 0), color = "grey80", linetype = "dotted") +
  geom_point(size = 3) +
  facet_grid(Group ~ ., scales = "free_y", space = "free_y") +
  scale_color_manual(values = col_groups) +
  theme_bw(base_size = 12) +
  labs(title = "Cleveland Dotplot",
       subtitle = "Sequencing depth outlier check",
       x = "Total Reads (Depth)", y = NULL) +
  theme(legend.position = "bottom", # Moved the legend here for the final panel!
        strip.background = element_rect(fill = "grey90"),
        strip.text = element_text(face = "bold", size = 11))

p_cleveland + p_rar
ggsave("plot_DE_Combined.pdf", p_cleveland + p_rar, width = 14, height = 6)
cat("Data Exploration (DE) panel saved successfully!\n")


# -----------------------------------------------------------------------------
# 3. ALPHA DIVERSITY (Indices & Asymptotic Estimators) [Ref 210]
# -----------------------------------------------------------------------------
alpha_est <- estimateR(comm) # Chao1 and ACE
alpha <- data.frame(
  Sample   = rownames(comm),
  Group    = metadata$Group,
  Richness = specnumber(comm),
  Chao1    = alpha_est[2, ], 
  Shannon  = diversity(comm, index = "shannon"),
  Simpson  = diversity(comm, index = "simpson"),
  Evenness = diversity(comm, index = "shannon") / log(specnumber(comm))
)

# Alpha Diversity Plot
alpha_long <- alpha %>% 
  pivot_longer(cols = c(Richness, Shannon, Simpson, Evenness), names_to = "Metric", values_to = "Value")

p_alpha <- ggplot(alpha_long, aes(x = Group, y = Value, fill = Group)) +
  geom_boxplot(alpha = 0.8, outlier.shape = NA, width = 0.5) +
  geom_jitter(aes(color = Group), width = 0.15, size = 2, alpha = 0.6) +
  facet_wrap(~Metric, scales = "free_y", nrow = 1) +
  scale_fill_manual(values = col_groups) +
  scale_color_manual(values = col_groups) +
  theme_bw(base_size = 12) + theme(legend.position = "none", plot.title = element_text(face="bold")) +
  labs(title = "Alpha Diversity Comparison", subtitle = "Wilcoxon Rank Sum Test", x = NULL)

p_alpha
ggsave("plot_alpha_diversity.pdf", p_alpha, width = 12, height = 5)
ggsave("plot_alpha_diversity.png", p_alpha, width = 12, height = 5, dpi = 300)
cat("Saved: plot_alpha_diversity\n")

# =============================================================================
# 4. TRANSFORMATIONS & BETA DIVERSITY (THE COMPARATIVE TEST)
# =============================================================================

# -----------------------------------------------------------------------------
# STRATEGY A: Hellinger + Bray-Curtis (Classic)
# -----------------------------------------------------------------------------
comm_hel <- decostand(comm, method = "hellinger")
dist_hel_bc <- vegdist(comm_hel, method = "bray")

# PERMDISP & PERMANOVA (Hellinger)
disper_hel <- betadisper(dist_hel_bc, metadata$Group)
perm_hel <- adonis2(dist_hel_bc ~ Group, data = metadata, permutations = 999)

cat("\n=== [STRATEGY A: HELLINGER] ===\n")
cat("PERMDISP (Dispersion) p-value:", anova(disper_hel)$`Pr(>F)`[1], "\n")
cat("PERMANOVA R²:", perm_hel$R2[1], " | p-value:", perm_hel$`Pr(>F)`[1], "\n")

# -----------------------------------------------------------------------------
# STRATEGY B: CLR + Euclidean (Aitchison Distance) - [Ref 201, 212]
# -----------------------------------------------------------------------------
# CLR does not accept exact zeros. We use the 'microbiome' package function 
# which internally adds a minimal pseudocount before applying the log.
comm_clr <- microbiome::transform(comm, transform = "clr")

# Euclidean distance on CLR data = Aitchison Distance!
dist_aitchison <- vegdist(comm_clr, method = "euclidean")

# PERMDISP & PERMANOVA (CLR / Aitchison)
disper_clr <- betadisper(dist_aitchison, metadata$Group)
perm_clr <- adonis2(dist_aitchison ~ Group, data = metadata, permutations = 999)

cat("\n=== [STRATEGY B: CLR / AITCHISON] ===\n")
cat("PERMDISP (Dispersion) p-value:", anova(disper_clr)$`Pr(>F)`[1], "\n")
cat("PERMANOVA R²:", perm_clr$R2[1], " | p-value:", perm_clr$`Pr(>F)`[1], "\n\n")

# =============================================================================
# 5. NMDS VISUALIZATION (Comparing both worlds side-by-side)
# =============================================================================

# NMDS - Hellinger
set.seed(42)
nmds_hel <- metaMDS(dist_hel_bc, trace = FALSE)
df_hel <- data.frame(NMDS1 = nmds_hel$points[, 1], NMDS2 = nmds_hel$points[, 2], 
                     Group = metadata$Group, Sample = metadata$Sample)

p_nmds_hel <- ggplot(df_hel, aes(x = NMDS1, y = NMDS2, color = Group, fill = Group)) +
  stat_ellipse(aes(group = Group), type = "t", level = 0.95, alpha = 0.15, geom = "polygon") +
  geom_point(size = 3, alpha = 0.8) +
  scale_color_manual(values = col_groups) + scale_fill_manual(values = col_groups) +
  annotate("text", x = Inf, y = Inf, hjust = 1.1, vjust = 1.5, size = 3, fontface = "italic",
           label = sprintf("PERMANOVA R² = %.3f (p = %.3f)\nPERMDISP p = %.3f", 
                           perm_hel$R2[1], perm_hel$`Pr(>F)`[1], anova(disper_hel)$`Pr(>F)`[1])) +
  theme_bw() + labs(title = "NMDS - Hellinger (Bray-Curtis)", x = "NMDS1", y = "NMDS2") +
  theme(legend.position = "none")

# NMDS - CLR (Aitchison)
set.seed(42)
nmds_clr <- metaMDS(dist_aitchison, trace = FALSE)
df_clr <- data.frame(NMDS1 = nmds_clr$points[, 1], NMDS2 = nmds_clr$points[, 2], 
                     Group = metadata$Group, Sample = metadata$Sample)

p_nmds_clr <- ggplot(df_clr, aes(x = NMDS1, y = NMDS2, color = Group, fill = Group)) +
  stat_ellipse(aes(group = Group), type = "t", level = 0.95, alpha = 0.15, geom = "polygon") +
  geom_point(size = 3, alpha = 0.8) +
  scale_color_manual(values = col_groups) + scale_fill_manual(values = col_groups) +
  annotate("text", x = Inf, y = Inf, hjust = 1.1, vjust = 1.5, size = 3, fontface = "italic",
           label = sprintf("PERMANOVA R² = %.3f (p = %.3f)\nPERMDISP p = %.3f", 
                           perm_clr$R2[1], perm_clr$`Pr(>F)`[1], anova(disper_clr)$`Pr(>F)`[1])) +
  theme_bw() + labs(title = "NMDS - CLR (Aitchison/Euclidean)", x = "NMDS1", y = "NMDS2")

p_nmds_hel + p_nmds_clr

# -----------------------------------------------------------------------------
# 6. TAXONOMIC COMPOSITION (TOP 15) & DIFFERENTIAL ABUNDANCE
# -----------------------------------------------------------------------------

# Proportional (%)
comm_rel <- sweep(comm, 1, rowSums(comm), "/") * 100
taxa_mean <- colMeans(comm_rel)
top15 <- names(sort(taxa_mean, decreasing = TRUE))[1:min(15, ncol(comm_rel))]

# Barplot preparation
comm_plot <- comm_rel %>%
  as.data.frame() %>%
  mutate(Sample = rownames(.), Group = metadata$Group) %>%
  pivot_longer(cols = -c(Sample, Group), names_to = "Taxon", values_to = "Abundance") %>%
  mutate(Taxon = ifelse(Taxon %in% top15, Taxon, "Other")) %>%
  group_by(Sample, Group, Taxon) %>%
  summarise(Abundance = sum(Abundance), .groups = "drop") %>%
  mutate(
    Taxon = factor(Taxon, levels = c(top15, "Other")),
    # PULLING THE CORRECT NUMERICAL ORDER HERE:
    Sample = factor(Sample, levels = levels(metadata$Sample)) 
  )

p_bar <- ggplot(comm_plot, aes(x = Sample, y = Abundance, fill = Taxon)) +
  geom_bar(stat = "identity") +
  facet_grid(~ Group, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = c(brewer.pal(12, "Paired"), brewer.pal(4, "Set3"), "grey80")) +
  theme_bw() + theme(axis.text.x = element_blank()) +
  labs(y = "Relative Abundance (%)", title = "Taxonomic Profile")

p_bar

# Differential Abundance (D vs P)
diff_df <- data.frame(Taxon = colnames(comm_rel), 
                      Diff = colMeans(comm_rel[metadata$Group == "Protected", ]) - 
                        colMeans(comm_rel[metadata$Group == "Degraded", ])) %>%
  filter(Taxon %in% top15) %>% arrange(Diff) %>%
  mutate(Taxon = factor(Taxon, levels = Taxon))

p_diff <- ggplot(diff_df, aes(x = Diff, y = Taxon, fill = Diff > 0)) +
  geom_col() +
  scale_fill_manual(values = c("TRUE" = "#4A9E8A", "FALSE" = "#E07B54"), 
                    labels = c("Higher in Degraded", "Higher in Protected")) +
  theme_bw() + labs(title = "Differential Abundance", x = "Difference in Mean % (Protected - Degraded)")

p_diff

# =============================================================================
# 7. COMPLETE EXPORT (Merging Comparative Plots)
# =============================================================================

# Now the panel shows the transformation comparison in the middle row!
p_painel <- (p_alpha) / 
  (p_nmds_hel | p_nmds_clr) / 
  (p_bar | p_diff) + 
  plot_layout(heights = c(1, 1.2, 1.5))

p_painel

ggsave("ALMA_Ecological_Comparison_Plot.pdf", p_painel, width = 18, height = 16)
cat("Analysis complete! The comparative panel was saved in 'ALMA_Ecological_Comparison_Plot.pdf'.\n")
