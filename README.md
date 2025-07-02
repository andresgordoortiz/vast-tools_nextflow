# VAST-TOOLS Splicing Analysis Pipeline

A Nextflow pipeline for analyzing alternative splicing patterns in RNA-Seq data using VAST-tools, specifically designed for FMN2/Spire knockout studies in oocytes.

## Overview

This pipeline processes RNA-Seq data to identify and quantify alternative splicing events using the VAST-tools suite. It supports both single-end and paired-end sequencing data, handles technical replicates, and provides comprehensive quality control and reporting.

## Features

- **Flexible Input**: Supports single-end, paired-end, and technical replicate data
- **Quality Control**: Integrated FastQC and MultiQC reporting
- **Read Processing**: Optional trimming with Trim Galore
- **Splicing Analysis**: VAST-tools alignment and splicing quantification
- **Automated Reporting**: R Markdown-based analysis reports
- **Containerized**: All processes run in Docker/Singularity containers

## Quick Start

```bash
# Clone the repository
git clone https://github.com/your-username/vast-tools_nextflow.git
cd vast-tools_nextflow

# Run the pipeline
nextflow run main.nf \
  --sample_csv samples.csv \
  --vastdb_path /path/to/vastdb \
  --data_dir /path/to/fastq_files \
  --outdir results
```

## Requirements

### Software Dependencies

- [Nextflow](https://nextflow.io/) (≥22.04.0)
- [Docker](https://www.docker.com/) or [Singularity](https://sylabs.io/singularity/)

### Reference Data

- **VAST-DB**: Download the appropriate species database from [VAST-tools](https://github.com/vastgroup/vast-tools)
  - For mouse (mm10): Use the Mm2 database
  - For human (hg38): Use the Hs2 database

## Input Files

### Sample Sheet (CSV)

Create a CSV file with the following required columns:

| Column | Description | Required |
|--------|-------------|----------|
| `sample` | Unique sample identifier | ✓ |
| `fastq_1` | Path to R1 FASTQ file (or single-end file) | ✓ |
| `fastq_2` | Path to R2 FASTQ file (paired-end only) | For paired-end |
| `type` | Data type: `single`, `paired`, or `technical_replicate` | ✓ |
| `group` | Group identifier for the sample | Optional |

#### Example Sample Sheet

```csv
sample,fastq_1,fastq_2,type,group
sample1,reads/sample1_R1.fastq.gz,reads/sample1_R2.fastq.gz,paired,control
sample2,reads/sample2_R1.fastq.gz,reads/sample2_R2.fastq.gz,paired,knockout
tech_rep1,reads/tech1_1.fastq.gz,reads/tech1_2.fastq.gz,technical_replicate,control
```

For technical replicates, you can add additional columns (`fastq_3`, `fastq_4`, etc.) to specify more files to concatenate.

## Parameters

### Mandatory Parameters

| Parameter | Description |
|-----------|-------------|
| `--sample_csv` | Path to the sample sheet CSV file |
| `--vastdb_path` | Path to the VAST-DB directory |
| `--data_dir` | Directory containing input FASTQ files |

### Optional Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--outdir` | `nextflow_results` | Output directory |
| `--species` | `mm10` | Species for VAST-tools alignment |
| `--project_name` | `oocyte_splicing_analysis` | Custom project name |
| `--skip_fastqc` | `false` | Skip FastQC quality control |
| `--skip_trimming` | `false` | Skip read trimming |
| `--skip_fastqc_in_trimming` | `false` | Skip FastQC within trim_galore |
| `--skip_rmarkdown` | `false` | Skip R Markdown report generation |
| `--rmd_file` | `scripts/R/notebooks/Oocyte_fmndko_spireko_complete.Rmd` | Path to R Markdown file |
| `--multiqc_config` | `null` | Path to MultiQC config file |

### Supported Species

The pipeline supports the following species (use the species code as `--species` parameter):

- `hg19`, `hg38` - Human
- `mm9`, `mm10` - Mouse
- `rn6` - Rat
- `dm6` - Drosophila
- `ce11` - C. elegans
- And many others (see VAST-tools documentation)

## Output Structure

```
results/
├── qc/
│   ├── fastqc/              # FastQC reports
│   ├── trimming_reports/    # Trim Galore reports
│   └── multiqc_report.html  # Combined QC report
├── trimmed_reads/           # Trimmed FASTQ files
├── vast_alignment/          # VAST-tools alignment outputs
├── inclusion_tables/        # Splicing inclusion tables
└── report/                  # R Markdown analysis report
```

### Key Output Files

- **`*_INCLUSION_LEVELS_FULL-*.tab`**: Main splicing quantification table
- **`multiqc_report.html`**: Comprehensive quality control report
- **`*.html`**: R Markdown analysis report (if enabled)

## Usage Examples

### Basic Usage

```bash
nextflow run main.nf \
  --sample_csv samples.csv \
  --vastdb_path /data/vastdb \
  --data_dir /data/fastq
```

### Skip Quality Control Steps

```bash
nextflow run main.nf \
  --sample_csv samples.csv \
  --vastdb_path /data/vastdb \
  --data_dir /data/fastq \
  --skip_fastqc \
  --skip_trimming
```

### Human Data Analysis

```bash
nextflow run main.nf \
  --sample_csv human_samples.csv \
  --vastdb_path /data/vastdb \
  --data_dir /data/human_fastq \
  --species hg38 \
  --project_name human_splicing_study
```

### Custom Resource Configuration

```bash
nextflow run main.nf \
  --sample_csv samples.csv \
  --vastdb_path /data/vastdb \
  --data_dir /data/fastq \
  -profile cluster \
  --max_cpus 32 \
  --max_memory 128.GB
```

## Container Information

The pipeline uses the following containers:

- **VAST-tools**: `andresgordoortiz/vast-tools:latest`
- **FastQC**: `quay.io/biocontainers/fastqc:0.11.9--0`
- **Trim Galore**: `https://depot.galaxyproject.org/singularity/trim-galore:0.6.9--hdfd78af_0`
- **MultiQC**: `multiqc/multiqc:latest`
- **R Analysis**: `andresgordoortiz/splicing_analysis_r_crg:v1.5`

## Resource Requirements

### Minimum Requirements

- **CPU**: 8 cores
- **Memory**: 16 GB RAM
- **Storage**: 100 GB free space

### Recommended for Large Datasets

- **CPU**: 16+ cores
- **Memory**: 64+ GB RAM
- **Storage**: 500+ GB free space

## Troubleshooting

### Common Issues

1. **VAST-DB Path Error**

   ```text
   ERROR: VASTDB directory /path/to/vastdb does not exist!
   ```

   **Solution**: Ensure the VAST-DB path is correct and accessible.

2. **Missing Sample Files**

   ```text
   ERROR: The specified sample CSV file does not exist
   ```

   **Solution**: Check the path to your sample CSV file and ensure all FASTQ files listed exist.

3. **Memory Issues**

   ```text
   Process exceeded available memory
   ```

   **Solution**: Increase memory allocation or use fewer parallel processes.

4. **Container Issues**

   ```text
   Unable to pull Docker image
   ```

   **Solution**: Ensure Docker/Singularity is properly installed and configured.

### Getting Help

- Check the [Nextflow documentation](https://nextflow.io/docs/latest/)
- Review [VAST-tools documentation](https://github.com/vastgroup/vast-tools)
- Open an issue on the project GitHub repository

## Citation

If you use this pipeline in your research, please cite:

1. **VAST-tools**: Tapial et al. (2017). An atlas of alternative splicing profiles and functional associations reveals new regulatory programs and genes that simultaneously express multiple major isoforms. Genome Research.

2. **Nextflow**: Di Tommaso et al. (2017). Nextflow enables reproducible computational workflows. Nature Biotechnology.

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Authors

- **Andrés Gordo** - *Initial work*

## Acknowledgments

- VAST-tools development team
- Nextflow community
- nf-core project for best practices

## Version History

- **v1.0.0** - Initial release
  - Support for single-end, paired-end, and technical replicate data
  - Integrated quality control and reporting
  - R Markdown analysis reports
