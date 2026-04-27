
# 🍄 Glomeromycota LSU rDNA Analysis Pipeline

A reproducible and scalable **Nextflow DSL2** pipeline for analyzing arbuscular mycorrhizal fungi (AMF) amplicon sequences (LSU rDNA). 

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

### 2. Download Reference Databases
Download the database files and place them in the `refs/` directory:
* `1806DB.ng.fasta` (Reference sequences)
* `1806DB.tax` (Taxonomy file)
* `629ref_aln_cut.fas` (Alignment reference)


### 3. Run the Pipeline
```bash
nextflow run main.nf \
    --input_dir "data/" \
    --reference_db "refs/1806DB.ng.fasta" \
    --taxonomy_db "refs/1806DB.tax" \
    --alignment_ref "refs/629ref_aln_cut.fas" \
    --outdir "./results" \
    -profile docker 

     or -profile singularity
```

---

## 🔧 Execution Profiles

You do not need to install Mothur or RAxML manually. The pipeline uses profiles to automatically pull the required containers or environments. Append one of these to your run command:

* `-profile docker` **(Recommended):** Runs processes in isolated Docker containers.
* `-profile singularity`: For HPC clusters where Docker is not allowed.


**Pro-Tip:** Add `-resume` to your command if a run is interrupted. Nextflow will pick up exactly where it left off!

---

## ⚙️ Fine-tuning your Analysis

You can adjust the pipeline's behavior by adding these flags to your command:

* **Filter by Length:** Use `--min_length 350` and `--max_length 600` if your amplicons are outside the standard LSU range.
* **Clustering Identity:** Change `--cluster_cutoff 0.03` to group sequences at 97% similarity (default is 0.02 for 98%).
* **Save Time:** Use `--skip_phylogenetic` if you only need the OTU table and taxonomic classification, skipping the heavy RAxML placement.
* **HPC Clusters:** If running on a cluster with Slurm, just use `-profile slurm,singularity`.

---


## 🧪 Downstream Analysis: TRAITAM

This project includes a dedicated module for trait-based analysis of Glomeromycota communities using **TRAITAM**. This analysis is performed as a post-processing step after the main Nextflow pipeline completes.

### What it does:
* Cross-references your OTU table with known ecological traits of AMF from Chaudhary, V.B., Nokes, L.F., González, J.B. et al. TraitAM, a global spore trait database for arbuscular mycorrhizal fungi. Sci Data 12, 588 (2025). https://doi.org/10.1038/s41597-025-04940-x.
* Generates functional reports and advanced visualizations.
* Provides a dashboard for interactive data exploration.

### How to run:
The TRAITAM scripts are located in the `traitam_analysis/` directory. Since it requires the outputs from the main pipeline (OTU table and Taxonomy), you should run it after the Nextflow execution is finished.

📖 **For detailed instructions and requirements, please read:** 👉 `traitam_analysis/USER_GUIDE.md`

## 🏷️ Main Parameters

You can customize the biological thresholds by adding these flags to your run command:

| Parameter | Description | Default |
| :--- | :--- | :--- |
| `--min_length` | Minimum sequence length allowed | `300` |
| `--max_length` | Maximum sequence length allowed | `580` |
| `--cluster_cutoff` | Distance threshold for OTU clustering | `0.02` |
| `--skip_chimera_check` | Set to `true` to skip VSEARCH chimera removal | `false` |
| `--skip_phylogenetic` | Set to `true` to skip RAxML/GAPPA placement | `false` |

---

## 📊 Output Structure

Once completed, the `--outdir` (default: `results/`) will contain several numbered folders following the pipeline's logic:

```text
results/
├── 01_contigs/         # Assembled paired-end reads
├── 05_classify/        # Taxonomic summary per sequence
├── 11_cluster/         # OTU grouping lists
├── 12_shared/          # The OTU Table (Shared format) ⭐
├── 13_oturep/          # Representative sequences for each OTU
└── 16_gappa/           # (If enabled) Tree-based taxonomic assignments
```

---

## 🔬 Interpreting Your Results

After a successful run, focus on these key files for your downstream ecological analysis (like Phyloseq in R):

### 1. The OTU Table (Shared File)
**Location:** `results/12_shared/*.shared`
This is your main abundance matrix. 
* Rows represent samples.
* Columns represent OTUs.
* Values are the read counts.

### 2. Taxonomic Classification
**Location:** `results/05_classify/*.taxonomy`
Shows the lineage (Phylum, Class, Order, Family, Genus, Species) assigned to your sequences based on the MaarjAM database, along with bootstrap confidence scores.

### 3. Representative Sequences
**Location:** `results/13_oturep/*.rep.fasta`
Contains the single most abundant sequence for each generated OTU, useful for BLASTing novel sequences or aligning with external datasets.

### 4. Phylogenetic Placement (GAPPA)
**Location:** `results/16_gappa/*.taxonomy.tsv`
If you didn't skip the phylogenetic step, this file provides highly accurate taxonomic assignments by placing your representative sequences onto a reference tree using maximum likelihood.

---

## 📚 Citations

If you use this pipeline in your research, please cite:

* **Nextflow:** Di Tommaso, P., et al. (2017). *Nature Biotechnology*.
* **Mothur:** Schloss, P.D., et al. (2009). *Applied and Environmental Microbiology*.
* **RAxML-NG:** Kozlov, A. M., et al. (2019). *Bioinformatics*.
* **GAPPA:** Czech, L., et al. (2020). *Bioinformatics*.
```
## 🐛 Troubleshooting

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

---