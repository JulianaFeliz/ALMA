# TraitAM Analysis - Complete User Guide

**Version:** 1.0.0  
**Last Updated:** April 17, 2026

---

## Table of Contents

1. [Quick Start](#quick-start)
2. [Complete Workflow](#complete-workflow)
3. [Input Requirements](#input-requirements)
4. [Output Files](#output-files)
5. [Interpreting Results](#interpreting-results)
6. [Interactive Dashboard](#interactive-dashboard)
7. [Advanced Options](#advanced-options)
8. [Troubleshooting](#troubleshooting)
9. [Examples](#examples)

---

## Quick Start

### One-Command Analysis

Run the complete pipeline from GAPPA output to final HTML report:

```bash
./run_complete_analysis.sh \
    --input-dir /path/to/gappa_output \
    --traitam-dir data \
    --output-dir results \
    --sample-name my_sample \
    --launch-dashboard
```

This executes:
1. ✅ Taxonomic matching (OTUs → TraitAM species)
2. ✅ Trait extraction from TraitAM database
3. ✅ Functional diversity calculation (FRic, FEve, FDiv, FDis, RaoQ)
4. ✅ Phylogenetic analysis (if tree available)
5. ✅ Comprehensive visualizations
6. ✅ HTML report generation
7. ✅ Interactive Shiny dashboard launch

**Typical runtime:** 5-10 minutes for 1000 OTUs

---

## Complete Workflow

### Step 1: Prepare Data

**Required files:**

1. **GAPPA output** (from phylogenetic placement):
   ```
   gappa_output/
   ├── sample1.tsv    # Taxonomy table with OTU placements
   └── ...
   ```

2. **TraitAM database** (download from Dryad):
   ```bash
   # Download TraitAM data
   cd data/
   wget https://datadryad.org/stash/downloads/file_stream/2891899 \
        -O DataRecord_1.txt
   wget https://datadryad.org/stash/downloads/file_stream/2891900 \
        -O Glomeromycota_tree.tre
   ```

3. **Abundance data** (optional, if not in GAPPA output):
   ```
   abundance.tsv with columns:
   - name: OTU identifier
   - abundance: Read count
   ```

### Step 2: Run Analysis

**Option A: Complete automated pipeline**
```bash
./run_complete_analysis.sh \
    -i gappa_output/ \
    -t data/ \
    -o results/ \
    -s sample1
```

**Option B: Step-by-step execution**
```bash
# 1. Main analysis
Rscript scripts/run_traitam_analysis.R \
    --gappa-dir gappa_output/ \
    --traitam-dir data/ \
    --output-dir results/ \
    --sample-id sample1

# 2. Generate visualizations
Rscript scripts/generate_visualizations.R sample1 results/

# 3. Create HTML report
Rscript scripts/generate_report.R sample1 results/ data/

# 4. Launch interactive dashboard
Rscript scripts/launch_dashboard.R results/
```

### Step 3: Explore Results

**Primary outputs:**
- `results/07_reports/sample1_analysis_report.html` - Open in browser
- `results/06_visualizations/` - All publication-ready plots
- `results/03_functional_diversity/sample1_FD_metrics.tsv` - FD statistics

**Interactive exploration:**
```bash
Rscript scripts/launch_dashboard.R results/
```

---

## Input Requirements

### GAPPA Output Format

**Expected columns:**
```
name          taxopath                                          LWR
OTU_0001      Fungi; Glomeromycota; Glomerales; Glomeraceae    0.95
OTU_0002      Fungi; Glomeromycota; Diversisporales; ...       0.87
```

**Required:**
- `name`: Unique OTU/ASV identifier
- `taxopath`: Full taxonomy string with semicolon separators
- `LWR`: Likelihood Weight Ratio (placement confidence, 0-1)

**Optional:**
- `abundance`: Read counts (if not provided separately)
- Additional metadata columns (preserved in output)

### TraitAM Data Files

**DataRecord_1.txt** format:
```
taxon                 Family           Genus           Spore_vol_um3  ...
Acaulospora bireticulata   Acaulosporaceae  Acaulospora     15234.5       ...
Rhizophagus irregularis    Glomeraceae      Rhizophagus     8456.2        ...
```

**Key trait columns:**
- `Spore_vol_um3`: Spore volume (μm³)
- `Aspect_ratio`: Length/width ratio
- `Color_score`: Color intensity (0-4 scale)
- `Wall_investment`: Wall thickness index
- `Ornamentation_height_um`: Surface ornamentation height (μm)

**Glomeromycota_tree.tre**: Newick format phylogenetic tree with species names matching `taxon` column

---

## Output Files

### Directory Structure

```
results/
├── 01_matched_taxa/              # Taxonomic matching results
│   ├── sample1_matched_species.tsv      # OTU → species mapping
│   └── sample1_match_summary.txt        # Match statistics
│
├── 02_trait_data/                # Trait + abundance integration
│   └── sample1_traits_abundance.tsv     # Complete trait dataset
│
├── 03_functional_diversity/      # FD metrics
│   ├── sample1_FD_metrics.tsv           # FRic, FEve, FDiv, FDis, RaoQ
│   ├── sample1_CWM.tsv                  # Community weighted means
│   └── sample1_functional_space.pdf     # Functional space plot
│
├── 04_phylogenetic_analysis/     # Phylogenetic metrics
│   ├── sample1_phylo_signal.tsv         # Blomberg's K, Pagel's λ
│   ├── sample1_phylo_diversity.tsv      # PD, MPD, MNTD
│   ├── sample1_tree_comparison.pdf      # Phylogenetic trees
│   └── sample1_pruned_tree.tre          # Sample-specific tree
│
├── 05_rarefaction/               # Rarefaction analysis
│   ├── sample1_rarefaction_data.tsv
│   └── sample1_rarefaction_plot.pdf
│
├── 06_visualizations/            # Publication-ready figures
│   ├── sample1_trait_distribution.pdf   # Trait histograms
│   ├── sample1_trait_correlations.pdf   # Correlation heatmap
│   ├── sample1_PCA_biplot.pdf           # PCA of traits
│   ├── sample1_abundance_traits.pdf     # Trait-abundance plots
│   └── sample1_family_composition.pdf   # Taxonomic composition
│
└── 07_reports/                   # Final reports
    ├── sample1_analysis_report.html     # Interactive HTML report
    └── sample1_summary.txt              # Text summary
```

### Key Output Files Explained

**`matched_species.tsv`:**
```
name       taxopath              LWR   species_name              match_type
OTU_0001   Fungi; Glomeromycota  0.95  Rhizophagus irregularis  exact_species
OTU_0002   Fungi; Glomeromycota  0.87  Funneliformis mosseae    fuzzy_species_0.90
```

**`FD_metrics.tsv`:**
```
Sample    FRic      FEve   FDiv   FDis    RaoQ    nb_sp
sample1   0.8234    0.712  0.845  2.456   3.234   42
```

**`CWM.tsv`:**
```
Sample    Spore_vol_um3  Aspect_ratio  Color_score  Wall_investment
sample1   12456.7        1.89          2.34         0.156
```

---

## Interpreting Results

### Functional Diversity Metrics

#### FRic (Functional Richness)
- **Range:** 0 to ∞
- **Meaning:** Total volume of functional trait space occupied by the community
- **Interpretation:**
  - High FRic → Diverse trait combinations present
  - Low FRic → Species occupy similar functional niches

#### FEve (Functional Evenness)
- **Range:** 0 to 1
- **Meaning:** How evenly species are distributed in trait space
- **Interpretation:**
  - FEve = 1 → Perfectly even distribution
  - FEve < 0.5 → Some regions of trait space are over-represented

#### FDiv (Functional Divergence)
- **Range:** 0 to 1
- **Meaning:** Degree to which abundance is distributed at extremes of trait space
- **Interpretation:**
  - High FDiv → Abundant species have divergent traits
  - Low FDiv → Abundant species cluster in trait space

#### FDis (Functional Dispersion)
- **Range:** 0 to ∞
- **Meaning:** Mean distance of species to the centroid in multidimensional trait space
- **Interpretation:**
  - High FDis → Large spread of trait values
  - Low FDis → Species cluster around mean trait values

#### RaoQ (Rao's Quadratic Entropy)
- **Range:** 0 to ∞
- **Meaning:** Expected trait dissimilarity between two randomly selected individuals
- **Interpretation:**
  - Accounts for both trait differences AND species abundances
  - Higher values → More functionally diverse community

### Community Weighted Means (CWM)

**Formula:**
```
CWM_trait = Σ (trait_value_i × abundance_i) / Σ abundance_i
```

**What it tells you:**
- The **dominant trait value** in your community
- Weighted by species abundance (common species have more influence)

**Examples:**
- CWM spore volume = 15,000 μm³ → Large-spored species dominate
- CWM wall investment = 0.20 → Thick-walled species are abundant
- CWM color score = 3.5 → Dark-colored spores are common

**Ecological significance:**
- CWM reflects the **functional identity** of the community
- Can be linked to environmental conditions
- Useful for comparing communities across gradients

### Phylogenetic Signal

#### Blomberg's K
- **K < 0.5:** Weak phylogenetic signal (traits evolve independently of phylogeny)
- **K ≈ 1.0:** Brownian motion model (neutral evolution, traits follow phylogeny)
- **K > 1.0:** Strong conservatism (closely related species are very similar)

**p-value interpretation:**
- p < 0.05 → Significant phylogenetic signal (traits are NOT randomly distributed across phylogeny)
- p > 0.05 → No significant signal (trait evolution independent of phylogeny)

#### Pagel's λ
- **λ = 0:** No phylogenetic signal (star phylogeny model)
- **λ = 1:** Traits perfectly follow Brownian motion along phylogeny
- **0 < λ < 1:** Intermediate signal

### Phylogenetic Diversity

#### Faith's PD (Phylogenetic Diversity)
- **What it measures:** Total branch length in the phylogenetic tree for your community
- **Higher values → More evolutionary history captured**
- **Compare PD to species richness:** High PD relative to species richness = phylogenetically diverse

#### MPD (Mean Pairwise Distance)
- **What it measures:** Average phylogenetic distance between all species pairs
- **High MPD → Phylogenetically overdispersed** (distantly related species coexist)
- **Low MPD → Phylogenetically clustered** (closely related species coexist)

#### MNTD (Mean Nearest Taxon Distance)
- **What it measures:** Average distance to the closest relative
- **High MNTD → Species lack close relatives** (competitive exclusion?)
- **Low MNTD → Closely related species co-occur** (niche conservatism?)

---

## Interactive Dashboard

### Launch Dashboard

```bash
# From results directory
Rscript scripts/launch_dashboard.R results/

# Custom port
Rscript scripts/launch_dashboard.R results/ --port 8080
```

Dashboard opens in your default web browser at `http://127.0.0.1:3838`

### Dashboard Features

#### 1. Overview Tab
- 📊 Summary statistics (species count, families, total reads)
- 📈 FD metrics overview
- 🥧 Family composition pie chart
- 📊 Community weighted means bar plot

#### 2. Species Explorer Tab
- 🔍 Searchable species table
- 📊 Species abundance ranking
- 📋 Detailed trait information per species
- 🎨 Interactive plots (click to zoom, hover for details)

#### 3. Trait Analysis Tab
- 📊 Trait distribution histograms
- 🔗 Trait correlation matrix (interactive heatmap)
- 📈 Trait vs. abundance scatterplots
- 📊 PCA biplot of functional space

#### 4. Functional Diversity Tab
- 📊 FD metrics comparison (if multiple samples)
- 🎨 Functional space visualization (3D PCA)
- 📊 Convex hull plots

#### 5. Community Comparison Tab (multi-sample)
- 📊 Side-by-side FD metrics comparison
- 📊 CWM heatmap across samples
- 📈 Trait distribution overlays

#### 6. Data Tables Tab
- 📋 Complete trait dataset (downloadable)
- 🔍 Filtered and sortable tables
- 📥 Export to CSV

### Dashboard Tips

✅ **Use filters** to subset data by family, abundance threshold, etc.

✅ **Hover over points** in plots to see species names and values

✅ **Click and drag** to zoom in on plot regions

✅ **Double-click** to reset zoom

✅ **Download plots** using the camera icon (plotly figures)

---

## Advanced Options

### Custom Phylogenetic Tree

Use your own phylogenetic tree instead of TraitAM's:

```bash
./run_complete_analysis.sh \
    -i gappa_output/ \
    -t data/ \
    -o results/ \
    --phylogeny my_custom_tree.tre
```

**Requirements:**
- Tree must be in Newick format
- Tip labels must match species names in TraitAM database
- Tree should include all or most species in your dataset

### Skip Components

**Skip phylogenetic analysis** (faster):
```bash
./run_complete_analysis.sh ... --skip-phylo
```

**Skip visualizations** (minimal output):
```bash
./run_complete_analysis.sh ... --skip-visualizations
```

**Skip report generation**:
```bash
./run_complete_analysis.sh ... --skip-report
```

### Batch Processing Multiple Samples

```bash
#!/bin/bash
# batch_analysis.sh

SAMPLES=(sample1 sample2 sample3 sample4)

for sample in "${SAMPLES[@]}"; do
    echo "Processing $sample..."
    ./run_complete_analysis.sh \
        -i gappa_output/${sample}/ \
        -t data/ \
        -o results/${sample}/ \
        -s ${sample} \
        --skip-dashboard
done

echo "All samples processed. Launching combined dashboard..."
Rscript scripts/launch_dashboard.R results/
```

### Parallel Processing

```bash
# Use GNU parallel for faster processing
parallel ./run_complete_analysis.sh \
    -i gappa_output/{}/  \
    -t data/ \
    -o results/{}/ \
    -s {} \
    ::: sample1 sample2 sample3
```

---

## Troubleshooting

### Problem: No species matched

**Symptoms:**
```
WARNING: Only 2 species matched to TraitAM database
```

**Diagnosis:**
```bash
# Check GAPPA taxonomy format
head -n 5 gappa_output/taxonomy.tsv

# Verify TraitAM loaded correctly
Rscript -e "
source('scripts/01_load_traitam.R')
traitam <- load_traitam_data('data/DataRecord_1.txt')
cat('TraitAM species:', nrow(traitam), '\n')
head(traitam$taxon)
"
```

**Solutions:**
1. **Check taxonomy path format** - should be: `Fungi; Glomeromycota; Order; Family; Genus; species`
2. **Enable fuzzy matching** - already enabled by default, threshold = 0.85
3. **Lower genus-level matching** - edit `scripts/02_match_taxonomy.R` to allow more genus matches

### Problem: Phylogenetic analysis fails

**Symptoms:**
```
Error: Tree does not contain species from dataset
```

**Diagnosis:**
```bash
# Check tree tip labels
Rscript -e "
library(ape)
tree <- read.tree('data/Glomeromycota_tree.tre')
cat('Tree tips:', length(tree\$tip.label), '\n')
cat('Sample species:', tree\$tip.label[1:10], '\n')
"

# Check species overlap
Rscript -e "
tree <- ape::read.tree('data/Glomeromycota_tree.tre')
traits <- readr::read_tsv('results/02_trait_data/sample1_traits_abundance.tsv')
overlap <- sum(traits\$species_name %in% tree\$tip.label)
cat('Species overlap:', overlap, '/', length(unique(traits\$species_name)), '\n')
"
```

**Solutions:**
1. **Skip phylogenetic analysis** if tree-species overlap is low (<5 species)
2. **Use custom tree** that includes your detected species
3. **Check species name format** - must match exactly between data and tree

### Problem: R package errors

**Symptoms:**
```
Error: package 'FD' is not installed
```

**Solution:**
```r
# Install missing packages
install.packages(c(
    "tidyverse", "data.table", "stringdist",
    "FD", "vegan", "ape", "phytools", "picante",
    "ggplot2", "ggpubr", "ggtree", "patchwork",
    "rmarkdown", "knitr", "DT", 
    "shiny", "shinydashboard", "plotly"
))

# ggtree requires Bioconductor
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("ggtree")
```

### Problem: Memory errors

**Symptoms:**
```
Error: vector memory exhausted (limit reached?)
```

**Solutions:**

1. **Increase R memory limit:**
```bash
Rscript --max-mem-size=16G scripts/run_traitam_analysis.R ...
```

2. **Reduce sample size for testing:**
```bash
# Take first 500 OTUs
head -n 501 gappa_output/taxonomy.tsv > gappa_output/test_taxonomy.tsv
```

3. **Filter low-abundance OTUs:**
```r
# Add to analysis script
filtered_data <- gappa_data %>% filter(abundance > 10)
```

### Problem: Dashboard won't launch

**Symptoms:**
```
Error: Shiny app failed to start
```

**Diagnosis:**
```bash
# Check if results exist
ls -lh results/02_trait_data/

# Test loading data
Rscript -e "
library(tidyverse)
traits <- read_tsv('results/02_trait_data/sample1_traits_abundance.tsv')
cat('Loaded', nrow(traits), 'records\n')
"
```

**Solutions:**
1. Ensure analysis completed successfully
2. Check that all required files exist in results/
3. Try running dashboard with explicit path:
```bash
Rscript scripts/launch_dashboard.R $(pwd)/results/
```

---

## Examples

### Example 1: Basic Analysis

```bash
# Download example data
wget https://example.com/gappa_example.tar.gz
tar -xzf gappa_example.tar.gz

# Run complete pipeline
./run_complete_analysis.sh \
    -i gappa_example/ \
    -t data/ \
    -o results_example/ \
    -s example_sample

# View report
firefox results_example/07_reports/example_sample_analysis_report.html
```

### Example 2: Multiple Samples Comparison

```bash
# Analyze 3 forest soil samples
for sample in forest1 forest2 forest3; do
    ./run_complete_analysis.sh \
        -i data/${sample}_gappa/ \
        -t data/traitam/ \
        -o results/${sample}/ \
        -s ${sample}
done

# Launch dashboard for comparison
Rscript scripts/launch_dashboard.R results/
```

### Example 3: Custom Workflow

```r
# Custom R script combining steps

library(tidyverse)
source("scripts/01_load_traitam.R")
source("scripts/02_match_taxonomy.R")

# Load data
traitam <- load_traitam_data("data/DataRecord_1.txt")
gappa <- read_tsv("gappa_output/sample1.tsv")

# Custom matching (e.g., only exact matches)
matched <- match_taxonomy_exact_only(gappa, traitam)

# Continue with analysis
source("scripts/03_extract_traits.R")
traits <- extract_traits(matched, traitam)

# Calculate FD
source("scripts/04_functional_diversity.R")
fd_results <- calculate_fd_metrics(traits)

print(fd_results)
```

### Example 4: Publication-Ready Figures

```r
# Generate custom high-resolution figures

library(tidyverse)
library(ggpubr)

# Load results
traits <- read_tsv("results/02_trait_data/sample1_traits_abundance.tsv")
fd <- read_tsv("results/03_functional_diversity/sample1_FD_metrics.tsv")

# Custom publication plot
p1 <- ggplot(traits, aes(x = Spore_vol_um3, y = log10(abundance + 1))) +
    geom_point(aes(color = Family), size = 3, alpha = 0.7) +
    geom_smooth(method = "lm", se = TRUE) +
    theme_pubr() +
    labs(
        x = "Spore Volume (μm³)",
        y = "log₁₀(Abundance + 1)",
        title = "Trait-Abundance Relationship"
    )

ggsave("publication_figure.pdf", p1, width = 8, height = 6, dpi = 300)
```

---

## FAQ

**Q: Can I use this pipeline with non-AMF fungi?**  
A: No, this pipeline is specifically designed for AMF (Glomeromycota) using the TraitAM database. For other fungi, you would need a different trait database.

**Q: What if I don't have GAPPA output?**  
A: You can provide any taxonomy table with `name`, `taxopath`, and `abundance` columns. The GAPPA-specific columns (`LWR`) are optional.

**Q: Can I add custom traits to the analysis?**  
A: Yes! Modify `scripts/03_extract_traits.R` to include additional trait columns from your own data or from TraitAM's supplementary tables.

**Q: How do I cite this pipeline?**  
A: Cite the TraitAM database (Chaudhary et al. 2024) and the FD package (Laliberté & Legendre 2010). Mention this pipeline was created by Seqera AI.

**Q: Can I run this on a cluster?**  
A: Yes! Submit as an R batch job. Use `--threads` option to parallelize within R:
```bash
#!/bin/bash
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G

Rscript scripts/run_traitam_analysis.R \
    --threads 8 \
    ...other options...
```

---

## Additional Resources

**TraitAM Database:**
- Paper: Chaudhary et al. (2024)
- Data: https://datadryad.org/dataset/doi:10.5061/dryad.6hdr7sr8z

**Functional Diversity:**
- FD package vignette: https://cran.r-project.org/package=FD
- Laliberté & Legendre (2010) Ecology

**Phylogenetic Diversity:**
- picante package: https://cran.r-project.org/package=picante
- Kembel et al. (2010) Bioinformatics

**R Shiny:**
- Dashboard tutorial: https://shiny.rstudio.com/tutorial/

---

**Questions? Issues? Suggestions?**  
Open an issue on GitHub or contact seqera-ai@example.com

**Happy analyzing! 🍄🔬📊**
