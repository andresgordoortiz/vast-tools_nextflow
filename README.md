# VAST-TOOLS Splicing Analysis Pipeline

A Nextflow pipeline for quantifying alternative splicing in RNA-Seq data using [VAST-tools](https://github.com/vastgroup/vast-tools).

---

## 📋 Quick Start Guide (CRG Cluster)

This guide will help you run the pipeline step by step. **No programming experience required!**

### Step 1: Prepare Your Sample Sheet

Create a CSV file (like an Excel spreadsheet saved as CSV) that describes your samples.

**📁 Example file:** See [data/samples_example.csv](data/samples_example.csv)

```csv
sample,fastq_1,fastq_2,type,group
control_rep1,control_rep1_R1.fastq.gz,control_rep1_R2.fastq.gz,paired,control
control_rep2,control_rep2_R1.fastq.gz,control_rep2_R2.fastq.gz,paired,control
treatment_rep1,treatment_rep1_R1.fastq.gz,treatment_rep1_R2.fastq.gz,paired,treatment
treatment_rep2,treatment_rep2_R1.fastq.gz,treatment_rep2_R2.fastq.gz,paired,treatment
```

#### CSV Columns Explained

| Column | Required? | Description |
|--------|-----------|-------------|
| `sample` | ✅ **Yes** | A unique name for each sample (no spaces, use underscores) |
| `fastq_1` | ✅ **Yes** | Name of your R1 FASTQ file (or single-end file) |
| `fastq_2` | For paired-end | Name of your R2 FASTQ file |
| `type` | ✅ **Yes** | Must be: `single`, `paired`, or `technical_replicate` |
| `group` | Optional | Group name for comparisons (e.g., "control", "treatment") |

> **💡 Tips for creating your CSV:**
> - You can create this in Excel and save as "CSV (Comma delimited)"
> - Make sure column names are exactly as shown (lowercase!)
> - File names should match your actual FASTQ files exactly
> - Don't include the full path in `fastq_1`/`fastq_2` - just the file names

---

### Step 2: Upload Your Data

1. Upload your **FASTQ files** to a folder on the cluster (e.g., `/users/yourname/projects/my_project/data/`)
2. Upload your **CSV sample sheet** to the same project folder

---

### Step 3: Run the Pipeline

Connect to the cluster and navigate to **this pipeline folder** (`vast-tools_nextflow/`), then run:

```bash
srun submit_nf.sh main.nf \
    --sample_csv /path/to/your/samples.csv \
    --data_dir /path/to/your/fastq_files/ \
    --vastdb_path /users/mirimia/projects/vast-tools/VASTDB/ \
    --species hg19 \
    --project_name my_project_name \
    -c nextflow.config \
    -work-dir /nfs/scratch01/yourlab/yourname/
```

#### 🔧 Replace the placeholders (YOU MUST CHANGE THESE):
- `/path/to/your/samples.csv` → Full path to your CSV file
- `/path/to/your/fastq_files/` → Folder containing your FASTQ files
- `/nfs/scratch01/yourlab/yourname/` → Your scratch folder for temporary files (use scratch!)
- `my_project_name` → A name for your project (no spaces!)
- `hg19` → Your species (see available options below)

---

## 📌 Parameters Reference

### Mandatory Parameters (You Must Specify)

| Parameter | Description | Example |
|-----------|-------------|--------|
| `--sample_csv` | Path to your CSV sample sheet | `/users/me/project/samples.csv` |
| `--data_dir` | Folder containing your FASTQ files | `/users/me/project/data/` |
| `-work-dir` | Scratch folder for temporary files (**use scratch!**) | `/nfs/scratch01/yourlab/yourname/` |
| `--vastdb_path` | `/users/aaljord/agordo/VASTDB/` | VAST-DB database (already downloaded on cluster) |
| `-c` | `nextflow.config` | Pipeline config file (included in this repo) |

### Optional Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--species` | `mm10` | Species code (see list below) |
| `--project_name` | `oocyte_splicing_analysis` | Name for your analysis |
| `--outdir` | `nextflow_results` | Where results will be saved |
| `--skip_fastqc` | `false` | Add this flag to skip quality control |
| `--skip_trimming` | `false` | Add this flag to skip read trimming |
| `--skip_compare` | `false` | Add this flag to skip differential splicing analysis |
| `--min_dPSI` | `10` | Minimum delta PSI threshold for differential splicing |
| `--min_range` | `5` | Minimum range threshold for differential splicing |
| `--skip_matt` | `false` | Add this flag to skip MATT feature analysis |
| `--matt_intron_length` | `150` | Intronic region length for SF1 binding site search |
| `--skip_betas` | `false` | Add this flag to skip betAS simulation analysis |
| `--betas_filter_n` | `10` | Minimum N for filtering events in betAS |
| `--betas_nsim` | `1000` | Number of simulations for betAS analysis |
| `--betas_npoints` | `500` | Number of points for betAS volcano plot |

### Available Species

| Code | Species |
|------|---------|
| `hg19` | Human (GRCh37) |
| `hg38` | Human (GRCh38) |
| `mm9` | Mouse (older assembly) |
| `mm10` | Mouse (GRCm38) |
| `rn6` | Rat |
| `dm6` | Drosophila |
| `ce11` | C. elegans |

> **📍 VAST-DB Location on CRG Cluster:**
> The database is already available at: `/users/aaljord/agordo/VASTDB/`
> Currently it olny contains the Mm10 genome build (Mm2). You may downlaod the rest from the official [vast-tools repo](https://github.com/vastgroup/vast-tools?tab=readme-ov-file#vastdb-libraries).

---

## 🔬 Differential Splicing Analysis (vast-tools compare)

**Automatic group comparison!** If your sample sheet contains a `group` column with 2 or more different groups, the pipeline will automatically run `vast-tools compare` for all pairwise group comparisons.

### How It Works

1. The pipeline detects groups from your CSV file
2. For each pair of groups (e.g., "control" vs "treatment"), it runs `vast-tools compare`
3. If your samples are **paired-end**, the `--paired` flag is automatically added

### Example with 2 Groups

With this sample sheet:
```csv
sample,fastq_1,fastq_2,type,group
ctrl_1,ctrl_1_R1.fastq.gz,ctrl_1_R2.fastq.gz,paired,control
ctrl_2,ctrl_2_R1.fastq.gz,ctrl_2_R2.fastq.gz,paired,control
treat_1,treat_1_R1.fastq.gz,treat_1_R2.fastq.gz,paired,treatment
treat_2,treat_2_R1.fastq.gz,treat_2_R2.fastq.gz,paired,treatment
```

The pipeline will automatically run:
```
vast-tools compare INCLUSION_TABLE.tab -a ctrl_1,ctrl_2 -b treat_1,treat_2 --min_dPSI 10 --min_range 5 --paired
```

### Example with 3 Groups

With 3 groups (control, treatment_A, treatment_B), the pipeline will run **3 comparisons**:
- control vs treatment_A
- control vs treatment_B
- treatment_A vs treatment_B

### Customizing Compare Parameters

You can adjust the thresholds:
```bash
sbatch submit_nf.sh main.nf \
    --sample_csv /path/to/samples.csv \
    --data_dir /path/to/data/ \
    --vastdb_path /users/aaljord/agordo/VASTDB/ \
    --species mm10 \
    --project_name my_project \
    --min_dPSI 15 \
    --min_range 10 \
    -c nextflow.config \
    -work-dir /nfs/scratch01/yourlab/yourname/
```

### Skipping Compare

If you don't want differential analysis, add `--skip_compare`:
```bash
srun submit_nf.sh main.nf \
    ... \
    --skip_compare
```

---


## 📊 Output Files

After the pipeline finishes, you'll find results in the `nextflow_results/` folder (or your custom `--outdir`):

```
nextflow_results/
├── pipeline_summary_*.txt       # ⭐ Pipeline execution summary (text)
├── pipeline_summary_*.json      # ⭐ Pipeline execution summary (JSON)
├── qc/
│   ├── fastqc/                  # Quality reports for each sample
│   ├── fastqc_trimmed/          # Quality reports after trimming
│   ├── trimming_reports/        # Trim Galore reports
│   └── multiqc_report.html      # Summary QC report (open in browser)
├── trimmed_reads/               # Cleaned FASTQ files
├── vast_alignment/              # VAST-tools alignment outputs
├── inclusion_tables/            # ⭐ Main results: splicing quantification tables
├── compare_results/             # 🔬 Differential splicing results (if groups defined)
│   └── compare_groupA_vs_groupB/
├── betas_analysis/              # 📈 betAS simulation-based analysis
│   ├── betas_filtered_events.RData    # Filtered splicing events
│   ├── betas_filtering_summary.txt    # Filtering statistics
│   └── groupA_vs_groupB/              # Per-comparison results
│       ├── betas_*_results.csv        # Full results table
│       └── betas_*_summary.txt        # Comparison summary
├── matt_references/             # Downloaded GTF and FASTA files for MATT
└── matt_analysis/               # 🧬 MATT feature analysis results
    └── groupA_vs_groupB/
        ├── input/               # Prepared input tables for MATT
        ├── exons/               # Exon feature analysis (PDF report + tables)
        └── introns/             # Intron feature analysis (PDF report + tables)
```

### Key Output Files

| File | Description |
|------|-------------|
| `pipeline_summary_*.txt` | Human-readable summary of all pipeline stages |
| `pipeline_summary_*.json` | Machine-readable summary (for automated processing) |
| `*_INCLUSION_LEVELS_FULL-*.tab` | Main splicing quantification table (PSI values for all events) |
| `compare_*/DiffAS-*.tab` | Differentially spliced events between groups |
| `betas_*_results.csv` | betAS simulation results with FDR-corrected p-values |
| `matt_analysis/*/exons/summary.pdf` | MATT exon feature comparison report with box plots |
| `matt_analysis/*/introns/summary.pdf` | MATT intron feature comparison report with box plots |
| `multiqc_report.html` | Summary quality control report |

---

## 📈 betAS Simulation-Based Analysis

The pipeline includes **betAS** (beta-distribution-based Alternative Splicing) analysis, which performs simulation-based statistical testing of differential splicing. This runs automatically after VAST-tools combine when groups are defined.

### What betAS Does

1. **Filters events**: Only keeps events with sufficient coverage (N ≥ 10 by default)
2. **Simulates distributions**: Uses beta distributions to model PSI uncertainty
3. **Statistical testing**: Calculates FDR-corrected p-values for each event
4. **Generates results**: CSV files with dPSI, p-values, and significance calls

### betAS Output Files

| File | Description |
|------|-------------|
| `betas_filtered_events.RData` | R data object with filtered splicing events |
| `betas_filtering_summary.txt` | Events before/after filtering, sample counts |
| `betas_*_results.csv` | Full results table with dPSI, p-values, FDR |
| `betas_*_summary.txt` | Summary statistics per comparison |

### Customizing betAS Parameters

```bash
srun submit_nf.sh main.nf \
    ... \
    --betas_filter_n 15 \     # Require N ≥ 15 (stricter filtering)
    --betas_nsim 2000 \       # More simulations (more accurate, slower)
    --betas_npoints 1000      # More points for volcano plot
```

### Skipping betAS

If you don't need simulation-based analysis, add `--skip_betas`:
```bash
srun submit_nf.sh main.nf \
    ... \
    --skip_betas
```

> **⚠️ Resource Note:** betAS simulations are memory-intensive (60GB RAM per comparison). The pipeline handles this automatically.

---

## 🧬 MATT Feature Analysis

MATT (Motif Analysis using Transcript features) automatically analyzes sequence and structural features of differentially spliced exons and introns. It runs by default when you have groups defined in your sample CSV.

### What MATT Analyzes

**For exons (`cmpr_exons`):**
- Exon/intron lengths
- Splice site strength
- RNA secondary structure
- Binding motifs (SF1, etc.)
- GC content and more

**For introns (`cmpr_introns`):**
- Intron lengths
- Branch point features
- Polypyrimidine tract
- Splice site sequences

### MATT Output

MATT generates:
- **`summary.pdf`** - PDF report with box plots comparing feature distributions
- **Feature tables** - Detailed tables with all extracted features
- **Statistical tests** - P-values from Mann-Whitney U tests

### Skipping MATT

If you don't need feature analysis, add `--skip_matt`:
```bash
srun submit_nf.sh main.nf \
    ... \
    --skip_matt
```

---

## 📋 Pipeline Summary Report

At the end of each run, the pipeline generates a comprehensive summary report with statistics from all stages:

### Summary Contents

1. **VAST-tools Alignment Statistics**
   - Number of samples aligned
   - Output files per sample

2. **VAST-tools Combine Statistics**
   - Total splicing events in inclusion table
   - Event type breakdown (exon skipping, intron retention, microexons)

3. **VAST-tools Compare Statistics** (if groups defined)
   - Differential events per comparison
   - Up/down-regulated event counts

4. **betAS Analysis Statistics** (if enabled)
   - Events before/after filtering
   - Significant events (FDR < 0.05) per comparison
   - Processing time

### Output Files

- `pipeline_summary_{project_name}.txt` - Human-readable text summary
- `pipeline_summary_{project_name}.json` - Machine-readable JSON (for scripts/automation)

---

## 🗂️ How the Pipeline Organizes Results

The pipeline uses the `--outdir` parameter (default: `nextflow_results/`) to determine where all results are saved. Each process specifies its output location using `publishDir`:

```groovy
publishDir "${params.outdir}/inclusion_tables", mode: 'copy'
```

### Output Directory Structure

| Parameter | Controls | Default |
|-----------|----------|---------|
| `--outdir` | Base output directory | `nextflow_results/` |
| `--project_name` | Names output files | `oocyte_splicing_analysis` |

### How Results Are Organized

1. **QC reports** → `{outdir}/qc/`
2. **Trimmed reads** → `{outdir}/trimmed_reads/`
3. **Alignments** → `{outdir}/vast_alignment/`
4. **Inclusion tables** → `{outdir}/inclusion_tables/`
5. **Compare results** → `{outdir}/compare_results/`
6. **betAS results** → `{outdir}/betas_analysis/`
7. **MATT results** → `{outdir}/matt_analysis/`
8. **Summary** → `{outdir}/pipeline_summary_*.{txt,json}`

### Customizing Output Location

```bash
srun submit_nf.sh main.nf \
    ... \
    --outdir /users/me/projects/my_analysis/results/ \
    --project_name taxol_splicing_2026
```

This would create:
- `/users/me/projects/my_analysis/results/inclusion_tables/taxol_splicing_2026_INCLUSION_LEVELS_FULL-hg38.tab`
- `/users/me/projects/my_analysis/results/pipeline_summary_taxol_splicing_2026.txt`

---

## ❓ Troubleshooting

| Problem | Solution |
|---------|----------|
| "sample CSV file does not exist" | Check the path to your CSV file is correct |
| "data directory does not exist" | Check the path to your FASTQ folder |
| "VASTDB directory does not exist" | Use the path: `/users/mirimia/projects/vast-tools/VASTDB/` |
| Pipeline stuck or failed | Check the `.nextflow.log` file for error messages |

---

## 📚 Citation

If you use this pipeline in your research, please cite:

- **VAST-tools**: Tapial et al. (2017). *Genome Research*. [DOI](https://doi.org/10.1101/gr.220962.117)
- **Nextflow**: Di Tommaso et al. (2017). *Nature Biotechnology*. [DOI](https://doi.org/10.1038/nbt.3820)

---

## 👤 Author

**Andrés Gordo** - CRG Barcelona

For questions or issues, please open a GitHub issue or contact the author.
