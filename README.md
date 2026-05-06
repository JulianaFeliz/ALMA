# 🍄 ALMA - AMF LSU METABARCODING ANALYSIS

A reproducible and scalable **Nextflow DSL2** pipeline for analyzing arbuscular mycorrhizal fungi (AMF) amplicon sequences (LSU rDNA). Based on a curated database of AMF sequences 

This pipeline integrates **Mothur** for core sequence processing and **RAxML/GAPPA** for phylogenetic placement.

---

## ⚡ Quick Start

Get your analysis running in 3 simple steps:

### 1. Prepare Your Data
The pipeline automatically detects paired-end samples based on file names.
* Option A: Directory
Organize your files in a single folder like this:
    ```text
    data/
    ├── sample1_R1.fastq.gz
    ├── sample1_R2.fastq.gz
    ├── sample2_R1.fastq.gz
    └── sample2_R2.fastq.gz
    ```
    *Note: The pipeline uses the prefix before `_R1` as the sample ID.*

* Option B: Samplesheet
If you have files in different locations, use a CSV:
sample,fastq_1,fastq_2

### 2. Run the Pipeline
```bash
nextflow run main.nf \
    -input_dir "data/" \
    -outdir "./results" \
    -profile docker \
    -resume

     or -profile singularity
```
#### 🔧 Execution Profiles
You do not need to install Mothur or RAxML manually. The pipeline uses profiles to automatically pull the required containers or environments. Append one of these to your run command:

* `-profile docker` **(Recommended):** Runs processes in isolated Docker containers.
* `-profile singularity`: For HPC clusters where Docker is not allowed.


**Pro-Tip:** Add `-resume` to your command if a run is interrupted. Nextflow will pick up exactly where it left off!

---
## 🏷️ Main Parameters
The default parameters 
| Parameter | Description | Default |
| :--- | :--- | :--- |
| `--min_length` | Minimum sequence length allowed | `300` |
| `--max_length` | Maximum sequence length allowed | `580` |
| `--cluster_cutoff` | Distance threshold for OTU clustering | `0.03` |
| `--skip_chimera_check` | Set to `true` to skip VSEARCH chimera removal | `false` |
| `--skip_phylogenetic` | Set to `true` to skip RAxML/GAPPA placement | `false` |

### ⚙️ Fine-tuning your Analysis
You can customize the biological thresholds by adding these flags to your run command:

* **Filter by Length:** Use `--min_length 350` and `--max_length 600` if your amplicons are outside the standard LSU range.
   * This may be checked in by summary with this command: 
   ```bash

   ```
* **Clustering Identity:** Change `--cluster_cutoff 0.02` to group sequences at 98% similarity (default is 0.03 for 97%).
   * If you want to try the ASV approach, implement 0.00 for 100% clustering.
* **Save Time:** Use `--skip_phylogenetic` if you only need the OTU table and taxonomic classification, skipping the heavy RAxML placement.
* **HPC Clusters:** If running on a cluster with Slurm, just use `-profile slurm,singularity`.

---

## 📊 Output Structure

Once completed, the `--outdir` (default: `results/`) will contain several numbered folders following the pipeline's logic:

```text
results/
├── mothur/                  # Core sequence processing outputs
│   ├── shared/              # The OTU Abundance Table (.shared) ⭐
│   ├── taxonomy/            # Taxonomic classification per sequence
│   └── cluster/             # OTU clustering lists and binning data
├── raxml_results/           # Files used for phylogenetic inference
│   ├── query_only.fasta     # Your processed sequences (queries)
│   └── ref_only.fasta       # Reference sequences used for placement
├── gappa_results/           # High-accuracy phylogenetic classification
│   ├── ALMA_with_samples.jplace.gz  # Tree data with placements ⭐
│   ├── ALMA_OTU_abundance.tsv       # Final abundance matrix (TSV format)
│   └── [Sample].taxonomy.tsv        # Detailed taxonomy per sample
└── pipeline_info/           # Nextflow execution reports (runtime, memory)
```

---

# 🔬 Interpreting Your Results


## 1. The Final Matrix (GAPPA)
📁 Location: results/gappa_results/ALMA_OTU_abundance.tsv

This is the most important file for your downstream analysis in R or Python. It combines the OTU counts with the refined phylogenetic taxonomy.

| label | Group	| numOtus	| Acaulospora_mellea | Gigaspora_sp	| Sclerocarpum_amazonicum |
|-------|---------|-----------|---------|----------|---------|
| 0.02	| UFRN_D1	| 150	| 450	|0	| 120 |
|0.02	| UFRN_D2	| 150	| 320 |	15	| 90 |
| 0.02	| UFRN_P1	| 150 |	10 | 3600 | 0 |

Use this for: Diversity indices, NMDS, and Differential Abundance.

## 2. Phylogenetic Tree Data
📁 Location: results/gappa_results/ALMA_with_samples.jplace.gz

Contains the reference tree and the precise location where each of your OTUs was placed.

Use this for: Creating the circular trees with ggtree in R or vizualition on iTOL interface (https://itol.embl.de/)

---
## 🧪 Downstream Analysis: Advanced Ecological Statistics

This project includes a comprehensive, publication-ready R script (ALMA_ecology.R`) to perform advanced statistical ecology on your amplicon data. It bridges raw bioinformatic outputs with rigorous ecological theory, properly handling the compositional nature of microbiome datasets.

### 🌟 Features
* **Data Exploration (DE):** Automated Rarefaction curves and Cleveland dotplots for sequencing depth saturation and outlier detection.
* **Alpha Diversity:** Calculates Richness, Shannon, Simpson, Evenness, and asymptotic estimators (Chao1), validated by Wilcoxon rank-sum tests.
* **Compositional Beta Diversity:** Compares classical approaches (Hellinger + Bray-Curtis) against compositionally-aware methods (Centered Log-Ratio [CLR] + Aitchison distance).
* **Robust Hypothesis Testing:** Controls for Type I errors by performing **PERMDISP** (homogeneity of multivariate dispersions) prior to **PERMANOVA**.
* **Publication-ready Viz:** Automatically generates multipanel figures (using `ggplot2` and `patchwork`) including NMDS ordinations, top 15 taxa barplots, and differential abundance plots.

### 🛠️ How to Run

**1. Prepare your Metadata:**
The script is universally applicable and does not rely on hardcoded sample names. You must provide a simple tab-separated file named `metadata.tsv` in your working directory. 

Example of `metadata.tsv`:
```text
SampleID	Group
Sample_1	Degraded
Sample_2	Degraded
Sample_3	Protected
Sample_4	Protected
```
(Note: The exact order of samples in this file will dictate the plotting order in the final figures).

**2. Required Inputs:**
Place these files in the same directory as the R script:
`ALMA_OTU_abundance.tsv` (Output from the pipeline's GAPPA step)
 
 `metadata.tsv` (Created by you)

### 📊 Expected Outputs
The script outputs high-quality PDFs and PNGs: 

`plot_DE_Combined.pdf`: Quality control.

`plot_alpha_diversity.pdf`: Boxplots of alpha diversity metrics.

`ALMA_Ecological_Comparison_Plot.pdf`: A massive, fully arranged panel comparing Hellinger vs. CLR ordinations, alongside taxonomic profiles and differential abundance.

---

## 📚 Citations

If you use this pipeline in your research, please cite:

Felix, JRB; Magurno, F; Queiroz, MB; Goto, BT; Lima, JPSM. ALMA - AMF LSU METABARCODING ANALYSIS. Refining Community Analysis: Strategies and Challenges for Environmental Metabarcoding of general fungi and specific recommendations on Glomeromycota phylum, 2026.  (UPDATE TO ARTICLE TITLE and DOI)
 
 or also:
* **Nextflow:** Di Tommaso, P., et al. (2017). *Nature Biotechnology*.
* **Mothur:** Schloss, P.D., et al. (2009). *Applied and Environmental Microbiology*.
* **RAxML-NG:** Kozlov, A. M., et al. (2019). *Bioinformatics*.
* **GAPPA:** Czech, L., et al. (2020). *Bioinformatics*.

# 🐛 Troubleshooting

### Common Issues

#### ❌ "Cannot find input files"

**Solution:**
- Check that file paths in samplesheet are correct
- Ensure files follow `*_R1.fastq.gz` and `*_R2.fastq.gz` naming
- Use absolute paths or paths relative to work directory

#### ❌ "Reference files not found"

**Solution:**
```bash
# Check if databases exist
ls -lh refs/

# Should show:
# 1806DB.ng.fasta
# 1806DB.tax
# 629ref_aln_cut.fasta
```


#### ❌ "Container not found" or "Docker daemon not running"

**Solution:**
```bash
# Start Docker
sudo systemctl start docker

# Or use Singularity
nextflow run main.nf -profile singularity
```

#### ❌ "Out of memory" error

**Solution:**
```bash
# Reduce parallel processes
nextflow run main.nf --max_cpus 4 --max_memory '8.GB'
```

### Getting Help

1. **Check logs:**
   - Main log: `.nextflow.log`
   - Process logs: `work/` directory

2. **Detailed error reporting:**
   ```bash
   nextflow run main.nf -with-report report.html -with-timeline timeline.html
   ```

3. **Open an issue:**
   - Provide Nextflow version: `nextflow -version`
   - Include error message and command used
