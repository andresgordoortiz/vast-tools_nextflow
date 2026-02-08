#!/usr/bin/env nextflow

/*
==========================================
 VAST-TOOLS Splicing Analysis Pipeline
==========================================
 Authors: Andrés Gordo
 Date: April 2025

 Description:
 This pipeline processes RNA-Seq data to analyze alternative splicing patterns in FMN2/Spire knockouts using VAST-tools.
*/

nextflow.enable.dsl=2

// Default parameters
params.outdir = "$projectDir/nextflow_results"
params.vastdb_path = "/path/to/vastdb"
params.sample_csv = null   // CSV file with sample info (replaces reads_dir)
params.help = false
params.skip_fastqc = false  // Option to skip FastQC step
params.skip_trimming = false  // Option to skip trimming step
params.skip_fastqc_in_trimming = false  // Option to skip FastQC within trim_galore
params.skip_rmarkdown = true  // DEPRECATED - RMarkdown report generation is disabled
params.rmd_file = "$projectDir/scripts/R/notebooks/Oocyte_fmndko_spireko_complete.Rmd"  // DEPRECATED
params.prot_impact_url = "https://vastdb.crg.eu/downloads/mm10/PROT_IMPACT-mm10-v3.tab.gz"
params.species = "mm10"  // Default species for VAST-tools alignment
params.multiqc_config = null  // Path to multiqc config file (optional)
params.project_name = "oocyte_splicing_analysis"  // Custom project name for output files

// Change from default value to null to make it mandatory
params.data_dir = null  // Now mandatory - directory containing input FASTQ files

// VAST-tools compare parameters
params.skip_compare = false  // Skip differential splicing analysis (vast-tools compare)
params.min_dPSI = 10  // Minimum delta PSI for vast-tools compare
params.min_range = 5  // Minimum range for vast-tools compare

// MATT analysis parameters
params.skip_matt = false  // Skip MATT feature analysis (requires compare to run)
params.matt_intron_length = 150  // Length of intronic region to search for SF1 hits

// betAS analysis parameters
params.skip_betas = false  // Skip betAS simulation-based splicing analysis
params.betas_filter_n = 10  // Minimum N for filtering events in betAS
params.betas_nsim = 1000  // Number of simulations for betAS analysis
params.betas_npoints = 500  // Number of points for volcano plot preparation

// Display help message
def helpMessage() {
    log.info"""
    ============================================================
     FMN2/SPIRE ALTERNATIVE SPLICING ANALYSIS PIPELINE
    ============================================================
    Usage:

    nextflow run main.nf --sample_csv sample_sheet.csv --vastdb_path /path/to/vastdb --data_dir /path/to/data

    Mandatory Arguments:
      --sample_csv          Path to CSV file defining samples to analyze
      --vastdb_path         Path to VASTDB directory
      --data_dir            Path to directory containing FASTQ files referenced in sample_csv

    Sample CSV Format:
      The CSV file should contain the following columns:
      - sample: Unique sample identifier (required)
      - fastq_1: Path to FASTQ file (R1 for paired-end data) (required)
      - fastq_2: Path to FASTQ file for R2 (required for paired-end data)
      - type: One of 'single', 'paired', 'technical_replicate' (required)
      - group: Group identifier for the sample (optional)

      For technical replicates, additional columns (fastq_3, fastq_4, etc.) can
      be added to define more files to concatenate.

    Optional Arguments:
      --outdir              Output directory (default: ${params.outdir})
      --species             Species for VAST-tools alignment (default: ${params.species})
      --project_name        Custom project name for output files (default: ${params.project_name})
      --skip_fastqc         Skip FastQC quality control (default: ${params.skip_fastqc})
      --skip_trimming       Skip trimming of reads (default: ${params.skip_trimming})
      --skip_fastqc_in_trimming  Skip FastQC within trim_galore (default: ${params.skip_fastqc_in_trimming})
      --multiqc_config      Path to MultiQC config file (default: none)

    Differential Splicing (vast-tools compare) Arguments:
      --skip_compare        Skip differential splicing analysis (default: ${params.skip_compare})
      --min_dPSI            Minimum delta PSI threshold (default: ${params.min_dPSI})
      --min_range           Minimum range threshold (default: ${params.min_range})

    MATT Feature Analysis Arguments:
      --skip_matt           Skip MATT feature analysis (default: ${params.skip_matt})
      --matt_intron_length  Intronic region length for SF1 search (default: ${params.matt_intron_length})

    betAS Simulation Analysis Arguments:
      --skip_betas          Skip betAS simulation-based splicing analysis (default: ${params.skip_betas})
      --betas_filter_n      Minimum N for filtering events (default: ${params.betas_filter_n})
      --betas_nsim          Number of simulations (default: ${params.betas_nsim})
      --betas_npoints       Number of points for volcano plot (default: ${params.betas_npoints})

    Note: If your sample CSV contains a 'group' column with 2+ groups, the pipeline
    will automatically run vast-tools compare for all pairwise group comparisons.
    For paired-end samples, the --paired flag is automatically added.
    MATT analysis runs automatically after compare (if not skipped).

    Example Command:
      nextflow run main.nf --sample_csv samples.csv --vastdb_path /path/to/vastdb --data_dir /path/to/data --species mm10 --outdir results
    ============================================================
    """.stripIndent()
}

// Help message function is defined above, will be called from the workflow

// Define parameter validation function
def validateParameters() {
    if (params.vastdb_path == "/path/to/vastdb") {
        log.error "ERROR: No path to VAST-DB has been specified. Please set the --vastdb_path parameter to the location of your VAST-DB directory (e.g., --vastdb_path /path/to/vastdb)."
        exit 1
    }

    if (params.sample_csv == null) {
        log.error "ERROR: No sample CSV file has been specified. Please set the --sample_csv parameter to the location of your sample sheet CSV file (e.g., --sample_csv /path/to/samples.csv)."
        exit 1
    }

    if (params.data_dir == null) {
        log.error "ERROR: No data directory has been specified. Please set the --data_dir parameter to the location containing your FASTQ files (e.g., --data_dir /path/to/fastq_files)."
        exit 1
    }

    def dataDir = file(params.data_dir)
    if (!dataDir.exists()) {
        log.error "ERROR: The specified data directory '${params.data_dir}' does not exist. Please verify the path."
        exit 1
    }

    if (!dataDir.isDirectory()) {
        log.error "ERROR: The specified path '${params.data_dir}' is not a directory."
        exit 1
    }

    def sampleCsv = file(params.sample_csv)
    if (!sampleCsv.exists()) {
        log.error "ERROR: The specified sample CSV file '${params.sample_csv}' does not exist. Please verify the path."
        exit 1
    }

    // Validate RMarkdown file if not skipping
    if (!params.skip_rmarkdown) {
        def rmdFile = file(params.rmd_file)
        if (!rmdFile.exists()) {
            log.warn "WARNING: The specified RMarkdown file '${params.rmd_file}' does not exist. RMarkdown report will be skipped."
            params.skip_rmarkdown = true
        }
    }

    // Check if the CSV file is properly formatted (header check)
    def firstLine = sampleCsv.withReader { it.readLine() }
    def requiredColumns = ['sample', 'fastq_1', 'type']
    def headers = firstLine.split(',').collect { it.trim().toLowerCase() }

    def missingColumns = requiredColumns.findAll { !headers.contains(it) }
    if (missingColumns) {
        log.error "ERROR: The sample CSV is missing required columns: ${missingColumns.join(', ')}."
        log.error "Required format: sample,fastq_1,fastq_2,type,group"
        log.error "Where type must be 'single', 'paired', or 'technical_replicate'"
        exit 1
    }

    // Additional validation: Check that paired samples have fastq_2
    def csvData = sampleCsv.withReader { reader ->
        reader.readLines().drop(1) // Skip header
    }

    csvData.each { line ->
        def values = line.split(',')
        if (values.size() >= 4 && values[3].trim().toLowerCase() == 'paired') {
            if (values.size() < 3 || !values[2] || values[2].trim().isEmpty()) {
                log.error "ERROR: Sample '${values[0]}' is marked as 'paired' but missing fastq_2 column or value."
                exit 1
            }
        }
    }
}


process concatenate_technical_replicates {
    tag "Concatenating technical replicates: ${sample_id}"
    label 'process_medium'

    // Resource requirements
    cpus 4
    memory { 4.GB }
    time { 30.min }

    input:
    tuple val(sample_id), val(sample_type), path(fastq_files), val(group)

    output:
    tuple val(sample_id), val('single'), path("${sample_id}.fastq.gz"), val(group), emit: sample_data

    when:
    sample_type == 'technical_replicate'

    script:
    """
    echo "Concatenating ${fastq_files.size()} technical replicate files for sample ${sample_id}..."
    cat ${fastq_files.join(' ')} > ${sample_id}.fastq.gz
    echo "✓ Created merged file for ${sample_id}"
    """
}

// Process to handle paired-end reads
process prepare_paired_reads {
    tag "Preparing paired-end reads: ${sample_id}"

    // Resource requirements
    cpus 1
    memory { 2.GB }
    time { 10.min }

    input:
    tuple val(sample_id), val(sample_type), path(fastq_files), val(group)

    output:
    tuple val(sample_id), val('paired'), path("${sample_id}_R{1,2}.fastq.gz"), val(group), emit: sample_data

    when:
    sample_type == 'paired'

    script:
    """
    # Link files with standard names for downstream processing
    ln -s ${fastq_files[0]} ${sample_id}_R1.fastq.gz
    ln -s ${fastq_files[1]} ${sample_id}_R2.fastq.gz
    """
}

// Process to handle single-end reads
process prepare_single_reads {
    tag "Preparing single-end reads: ${sample_id}"

    // Resource requirements
    cpus 1
    memory { 2.GB }
    time { 10.min }

    input:
    tuple val(sample_id), val(sample_type), path(fastq_files), val(group)

    output:
    tuple val(sample_id), val('single'), path("${sample_id}.fastq.gz"), val(group), emit: sample_data

    when:
    sample_type == 'single'

    script:
    """
    # Link file with standard name for downstream processing
    ln -s ${fastq_files[0]} ${sample_id}.fastq.gz
    """
}

process verify_files {
    debug true

    // Resource requirements
    cpus 1
    memory { 1.GB }
    time { 5.min }

    input:
    path dir

    script:
    """
    echo "Verifying directory contents: ${dir}"
    echo "File count: \$(find ${dir} -type f | wc -l)"
    find ${dir} -type f | head -n 5
    """
}

process run_fastqc {
    tag "Quality control: ${sample_id}"
    label 'process_medium'
    publishDir "${params.outdir}/qc/fastqc", mode: 'copy'
    container 'quay.io/biocontainers/fastqc:0.11.9--0' // Using biocontainer for FastQC

    // Resource requirements
    cpus 4
    memory { 8.GB }
    time { 1.hour }

    input:
    tuple val(sample_id), val(sample_type), path(fastq_files), val(group)

    output:
    path "*_fastqc.{zip,html}", emit: fastqc_results

    script:
    """
    echo "Running FastQC quality control on sample ${sample_id}..."
    fastqc -t ${task.cpus} ${fastq_files}
    echo "FastQC analysis complete for ${sample_id}."
    """
}

process run_trim_galore {
    tag "Trimming: ${sample_id}"
    label 'process_medium'
    publishDir "${params.outdir}/trimmed_reads", mode: 'copy', pattern: "*.{fq.gz,fastq.gz}"
    publishDir "${params.outdir}/qc/trimming_reports", mode: 'copy', pattern: "*_trimming_report.txt"
    publishDir "${params.outdir}/qc/fastqc_trimmed", mode: 'copy', pattern: "*_fastqc.{html,zip}"
    container 'https://depot.galaxyproject.org/singularity/trim-galore:0.6.9--hdfd78af_0' // Using galaxyproject container for trim_galore

    // Retry configuration for handling large files that may need more resources
    maxRetries 3
    errorStrategy { task.attempt <= maxRetries ? 'retry' : 'terminate' }

    // Resource requirements - scale memory with retry attempts for large files
    // trim_galore uses cutadapt which scales well up to 4 cores
    cpus 4
    memory { 32.GB * task.attempt }
    time { 2.hours * task.attempt }

    input:
    tuple val(sample_id), val(sample_type), path(fastq_files), val(group)

    output:
    tuple val(sample_id), val(sample_type), path("*_trimmed.{fq.gz,fastq.gz}"), val(group), emit: trimmed_reads
    path "*_trimming_report.txt", emit: trim_log
    path "*_fastqc.{zip,html}", optional: true, emit: fastqc_results

    script:
    def fastqc_option = params.skip_fastqc_in_trimming ? "" : "--fastqc"
    def fastqc_args = params.skip_fastqc_in_trimming ? "" : "--fastqc_args '-t ${task.cpus}'"

    if (sample_type == 'paired') {
        // For paired-end data
        """
        echo "Trimming paired-end reads for sample ${sample_id}..."
        trim_galore ${fastqc_option} -j ${task.cpus} ${fastqc_args} \\
            --paired --quality 20 ${fastq_files[0]} ${fastq_files[1]} \\
            --basename ${sample_id}

        # More robust renaming with error checking
        if [ -f "${sample_id}_val_1.fq.gz" ]; then
            mv "${sample_id}_val_1.fq.gz" "${sample_id}_R1_trimmed.fq.gz"
        else
            echo "ERROR: Expected output file ${sample_id}_val_1.fq.gz not found"
            ls -la
            exit 1
        fi

        if [ -f "${sample_id}_val_2.fq.gz" ]; then
            mv "${sample_id}_val_2.fq.gz" "${sample_id}_R2_trimmed.fq.gz"
        else
            echo "ERROR: Expected output file ${sample_id}_val_2.fq.gz not found"
            ls -la
            exit 1
        fi

        echo "Trimming complete for ${sample_id}."
        """
    } else {
        // For single-end data - fix the naming logic
        """
        echo "Trimming single-end reads for sample ${sample_id}..."
        trim_galore ${fastqc_option} -j ${task.cpus} ${fastqc_args} \\
            --quality 20 ${fastq_files}

        # More robust file renaming
        TRIMMED_FILE=\$(find . -name "*_trimmed.fq.gz" -o -name "*_trimmed.fastq.gz" | head -n 1)
        if [ -n "\$TRIMMED_FILE" ]; then
            mv "\$TRIMMED_FILE" "${sample_id}_trimmed.fq.gz"
        else
            echo "ERROR: No trimmed output file found"
            ls -la
            exit 1
        fi

        echo "Trimming complete for ${sample_id}."
        """
    }
}

process run_multiqc {
    tag "MultiQC Report"
    label 'process_low'
    publishDir "${params.outdir}/qc", mode: 'copy'
    container 'multiqc/multiqc:latest' // Using the multiqc container

    // Resource requirements
    cpus 2
    memory { 4.GB }
    time { 20.min }

    input:
    path ('fastqc/*')
    path ('trim_logs/*')
    path ('vast_dirs/*'), stageAs: 'vast_*'

    output:
    path "multiqc_report.html", emit: report
    path "multiqc_report_data"

    script:
    def config_arg = params.multiqc_config ? "--config ${params.multiqc_config}" : ''
    """
    echo "Generating MultiQC report..."
    # Make sure the directory structure is preserved for MultiQC to recognize file types
    mkdir -p trim_galore_reports
    mkdir -p vast_out

    # Move files to appropriate directories based on naming patterns, but avoid moving files to themselves
    find . -name "*trimming_report.txt" -not -path "./trim_galore_reports/*" -exec mv {} ./trim_galore_reports/ \\;
    find . -name "*.tab" -o -name "*.log" | grep -v trimming | grep -v "./vast_out/" | xargs -I{} mv {} ./vast_out/ 2>/dev/null || true

    multiqc . ${config_arg} -f -n multiqc_report.html
    echo "MultiQC report generated."
    """
}

process prepare_vastdb {
    tag "Validate VASTDB"

    // Resource requirements
    cpus 1
    memory { 1.GB }
    time { 5.min }

    input:
    val vastdb_path

    output:
    val vastdb_path, emit: vastdb_path

    script:
    def species_dir = getVastdbDirName(params.species)
    """
    echo "Validating VASTDB structure at ${vastdb_path}"
    echo "Using species directory: ${species_dir} for ${params.species}"

    # Verify the VASTDB path exists
    if [ ! -d "${vastdb_path}" ]; then
        echo "ERROR: VASTDB directory ${vastdb_path} does not exist!"
        exit 1
    fi

    # Check if the species directory exists in VASTDB
    if [ ! -d "${vastdb_path}/${species_dir}" ]; then
        echo "ERROR: Species directory ${species_dir} not found in VASTDB ${vastdb_path}"
        echo "Available directories in VASTDB:"
        ls -la ${vastdb_path}/
        exit 1
    fi

    # Verify key files and directories exist
    echo "Verifying VASTDB structure..."
    ls -la ${vastdb_path}/${species_dir}/ || {
        echo "ERROR: Cannot list contents of species directory"
        exit 1
    }

    echo "VASTDB structure validated successfully"
    """
}

process align_reads {
    tag "VAST-tools alignment: ${sample_id}"
    label 'process_high'
    publishDir "${params.outdir}/vast_alignment", mode: 'copy', pattern: "vast_out/**"
    container 'andresgordoortiz/vast-tools:latest'

    // Retry configuration - retry twice if the process fails
    maxRetries 2
    errorStrategy { task.attempt <= maxRetries ? 'retry' : 'terminate' }

    // Resource requirements - scale both memory AND time with retries for large files
    // For 200M+ paired-end reads, expect 10-15 hours
    cpus 8
    memory { 30.GB * task.attempt }
    time { 25.hours + (3.hours * (task.attempt - 1)) }  // 15h -> 18h -> 21h on retries

    input:
    tuple val(sample_id), val(sample_type), path(fastq_files), val(group)
    val vastdb_path

    output:
    tuple val(sample_id), val(group), path("vast_out"), emit: alignment_output
    path "vast_out", emit: vast_out_dir

    script:
    def vast_options = "--IR_version 2 -c ${task.cpus} -n ${sample_id} -sp ${params.species} --verbose"

    if (sample_type == 'paired') {
        // Paired-end alignment
        """
        echo "VAST-tools alignment attempt ${task.attempt} for paired-end sample ${sample_id}..."
        mkdir -p vast_out/to_combine
        echo "Starting VAST-tools alignment for paired-end sample ${sample_id}..."
        echo "Using VASTDB path: ${vastdb_path}"

        # Debug info
        echo "Container environment:"
        env | grep SINGULARITY || true

        # Ensure VASTDB is correctly set, with multiple fallback options
        export VASTDB=/usr/local/vast-tools/VASTDB
        echo "VASTDB set to \$VASTDB"

        # Multiple command options for robustness
        if vast-tools align ${fastq_files[0]} ${fastq_files[1]} -o vast_out ${vast_options}; then
            echo "VAST-tools alignment completed successfully"
        else
            echo "First attempt failed, trying with explicit VASTDB path..."
            export VASTDB=${vastdb_path}
            echo "VASTDB now set to \$VASTDB"

            if vast-tools align ${fastq_files[0]} ${fastq_files[1]} -o vast_out ${vast_options}; then
                echo "VAST-tools alignment completed successfully with explicit path"
            else
                echo "Both alignment attempts failed - debugging information:"
                echo "Container filesystem:"
                ls -la /usr/local/vast-tools/ || true
                echo "VASTDB location:"
                ls -la \$VASTDB || true
                vast-tools --version || true
                exit 1
            fi
        fi

        # Verify alignment outputs and debug
        echo "Checking alignment outputs..."
        find vast_out -type f | sort | head -20
        echo "Directory structure:"
        find vast_out -type d | sort

        # Ensure all necessary files are in to_combine directory
        echo "Ensuring all necessary files are in to_combine directory..."
        mkdir -p vast_out/to_combine

        # Find all the different file types VAST-tools creates
        find vast_out -type f -name "*.eej*" -o -name "*.exskX" -o -name "*.info" -o -name "*.IR*" -o -name "*.mic*" -o -name "*.MULTI*" -o -name "*.tab" | while read file; do
            echo "Found file: \$file"
            cp "\$file" vast_out/to_combine/ || echo "Failed to copy \$file"
        done

        # Also look for files directly in cRPKM, to_combine and other subdirectories
        for subdir in vast_out/*/; do
            if [ -d "\$subdir" ] && [ "\$subdir" != "vast_out/to_combine/" ]; then
                echo "Checking subdirectory: \$subdir"
                find "\$subdir" -maxdepth 1 -type f | while read file; do
                    echo "Found file in subdir: \$file"
                    cp "\$file" vast_out/to_combine/ || echo "Failed to copy \$file"
                done
            fi
        done

        # List to_combine contents for verification
        echo "Contents of to_combine directory:"
        ls -la vast_out/to_combine/
        echo "File count in to_combine: \$(find vast_out/to_combine/ -type f | wc -l)"

        echo "VAST-tools alignment complete for ${sample_id}"
        """
    } else {
        // Single-end alignment - same changes as paired-end
        """
        echo "VAST-tools alignment attempt ${task.attempt} for single-end sample ${sample_id}..."
        mkdir -p vast_out/to_combine
        echo "Starting VAST-tools alignment for single-end sample ${sample_id}..."
        echo "Using VASTDB path: ${vastdb_path}"

        # Debug info
        echo "Container environment:"
        env | grep SINGULARITY || true

        # Ensure VASTDB is correctly set, with multiple fallback options
        export VASTDB=/usr/local/vast-tools/VASTDB
        echo "VASTDB set to \$VASTDB"

        # Multiple command options for robustness
        if vast-tools align ${fastq_files} -o vast_out ${vast_options}; then
            echo "VAST-tools alignment completed successfully"
        else
            echo "First attempt failed, trying with explicit VASTDB path..."
            export VASTDB=${vastdb_path}
            echo "VASTDB now set to \$VASTDB"

            if vast-tools align ${fastq_files} -o vast_out ${vast_options}; then
                echo "VAST-tools alignment completed successfully with explicit path"
            else
                echo "Both alignment attempts failed - debugging information:"
                echo "Container filesystem:"
                ls -la /usr/local/vast-tools/ || true
                echo "VASTDB location:"
                ls -la \$VASTDB || true
                vast-tools --version || true
                exit 1
            fi
        fi

        # Verify alignment outputs and debug
        echo "Checking alignment outputs..."
        find vast_out -type f | sort | head -20
        echo "Directory structure:"
        find vast_out -type d | sort

        # Ensure all necessary files are in to_combine directory
        echo "Ensuring all necessary files are in to_combine directory..."
        mkdir -p vast_out/to_combine

        # Find all the different file types VAST-tools creates
        find vast_out -type f -name "*.eej*" -o -name "*.exskX" -o -name "*.info" -o -name "*.IR*" -o -name "*.mic*" -o -name "*.MULTI*" -o -name "*.tab" | while read file; do
            echo "Found file: \$file"
            cp "\$file" vast_out/to_combine/ || echo "Failed to copy \$file"
        done

        # Also look for files directly in cRPKM, to_combine and other subdirectories
        for subdir in vast_out/*/; do
            if [ -d "\$subdir" ] && [ "\$subdir" != "vast_out/to_combine/" ]; then
                echo "Checking subdirectory: \$subdir"
                find "\$subdir" -maxdepth 1 -type f | while read file; do
                    echo "Found file in subdir: \$file"
                    cp "\$file" vast_out/to_combine/ || echo "Failed to copy \$file"
                done
            fi
        done

        # List to_combine contents for verification
        echo "Contents of to_combine directory:"
        ls -la vast_out/to_combine/
        echo "File count in to_combine: \$(find vast_out/to_combine/ -type f | wc -l)"

        echo "VAST-tools alignment complete for ${sample_id}"
        """
    }
}

process combine_results {
    tag "VAST-tools combine"
    label 'process_high'
    publishDir "${params.outdir}/inclusion_tables", mode: 'copy', pattern: '*INCLUSION_LEVELS_FULL*.tab'
    container 'andresgordoortiz/vast-tools:latest'

    // Container options:
    // DO NOT use --cleanenv or --no-home: vast-tools combine internally calls
    // R scripts (RI_MakeTablePIR.R) that need a functional R environment.
    // --cleanenv strips env vars causing R to hang waiting for stdin.
    // --no-home prevents R from finding user libraries.
    containerOptions {
        if (workflow.containerEngine == 'singularity') {
            return "--bind ${params.vastdb_path}:/usr/local/vast-tools/VASTDB"
        } else {
            return ''
        }
    }

    // Retry configuration
    maxRetries 3
    errorStrategy { task.exitStatus in [140, 139, 137, 143, 1] ? 'retry' : 'terminate' }

    // CRITICAL: vast-tools combine MUST use single core (--cores 1)
    // Multiple cores fork children that each load VASTDB into memory,
    // causing OOM kills (exit 140). See GitHub issue vastgroup/vast-tools#131
    // Also: the default in RunDBS_2.pl is already Ncores=1.
    cpus 1
    memory { 16.GB * task.attempt }  // 16GB -> 32GB -> 48GB on retries
    time { 2.hours * task.attempt }

    input:
    path vast_out_dirs, stageAs: "vast_*"
    val vastdb_path
    val output_name

    output:
    path "${output_name}_INCLUSION_LEVELS_FULL-${params.species}.tab", emit: inclusion_table

    script:
    """
    #!/bin/bash
    set -euo pipefail

    # ── Helper: timestamp logger ──
    ts() { echo "[\$(date '+%Y-%m-%d %H:%M:%S')] \$*"; }

    ts "=== VAST-tools combine START ==="
    ts "Attempt: ${task.attempt} | Memory: ${task.memory} | CPUs: ${task.cpus}"
    ts "Species: ${params.species} | Working dir: \$(pwd)"

    # Show system info for debugging
    ts "--- System Info ---"
    free -h 2>/dev/null || true
    df -h . 2>/dev/null | head -2 || true
    echo ""

    # ── Clean up leftover temp dirs from any prior failed attempts ──
    rm -rf raw_incl/ raw_reads/ tmp/ 2>/dev/null || true

    # ── Create to_combine/ directory ──
    mkdir -p to_combine

    # ── Collect VAST-tools align output files ──
    ts "--- Collecting align output files ---"
    for dir in vast_*; do
        [ -d "\$dir" ] || continue
        if [ -d "\$dir/to_combine" ]; then
            sample_files=\$(ls "\$dir/to_combine/" 2>/dev/null | head -1)
            ts "  \$dir -> \$(ls "\$dir/to_combine/" 2>/dev/null | wc -l) files"
            cp "\$dir/to_combine"/*.{eej2,exskX,IR2,IR.summary_v2.txt,micX,MULTI3X,info} to_combine/ 2>/dev/null || true
        else
            ts "  WARNING: \$dir/to_combine/ not found!"
        fi
    done

    # ── Verify collected files ──
    file_count=\$(find to_combine -type f | wc -l)
    ts "Total files in to_combine/: \$file_count"
    ts "File types breakdown:"
    for ext in eej2 exskX IR2 IR.summary_v2.txt micX MULTI3X info; do
        count=\$(ls to_combine/*.\$ext 2>/dev/null | wc -l)
        ts "  *.\$ext: \$count files"
    done

    if [ \$file_count -eq 0 ]; then
        ts "ERROR: No files found to combine!"
        exit 1
    fi

    # ── Verify VASTDB is accessible ──
    ts "--- VASTDB check ---"
    if [ -d "/usr/local/vast-tools/VASTDB" ]; then
        ts "VASTDB found at /usr/local/vast-tools/VASTDB"
        ls /usr/local/vast-tools/VASTDB/ 2>/dev/null | head -5
    else
        ts "WARNING: VASTDB not found at /usr/local/vast-tools/VASTDB!"
        ts "Checking VASTDB env: \${VASTDB:-not set}"
    fi

    # ── Background progress monitor ──
    # Monitors file creation in raw_incl/ and raw_reads/ to track stage progress
    (
        while true; do
            sleep 30
            now=\$(date '+%H:%M:%S')
            # Count output files being generated
            ri_count=\$(find raw_incl -type f 2>/dev/null | wc -l)
            rr_count=\$(find raw_reads -type f 2>/dev/null | wc -l)
            # Check what processes are running
            procs=\$(ps aux 2>/dev/null | grep -E 'perl|Rscript|vast-tools|Add_to|RI_Make|GetPSI|MakeTable' | grep -v grep | awk '{print \$11}' | xargs -r basename -a 2>/dev/null | sort -u | tr '\\n' ', ' || true)
            mem_used=\$(free -m 2>/dev/null | awk '/Mem:/{print \$3"/"\$2"MB"}' || echo "N/A")
            echo "[\$now] MONITOR: raw_incl=\$ri_count files, raw_reads=\$rr_count files | mem=\$mem_used | running: \${procs:-none}"
        done
    ) &
    MONITOR_PID=\$!
    trap "kill \$MONITOR_PID 2>/dev/null || true" EXIT

    # ── Run vast-tools combine ──
    ts "--- Running vast-tools combine ---"
    ts "Command: vast-tools combine -o . -sp ${params.species} --cores 1 --IR_version 2 --verbose"
    ts "NOTE: vast-tools buffers verbose output per sub-job (8 sequential jobs with --cores 1)"
    ts "       Progress monitor above tracks file creation every 30s"

    combine_start=\$(date +%s)

    # Run combine and tee stderr so we see output in real-time
    # vast-tools combine sends verbose output to stderr
    vast-tools combine -o . -sp ${params.species} --cores 1 --IR_version 2 --verbose 2>&1 | \
        while IFS= read -r line; do
            echo "[\$(date '+%H:%M:%S')] \$line"
        done
    combine_exit=\${PIPESTATUS[0]}

    combine_end=\$(date +%s)
    combine_duration=\$(( combine_end - combine_start ))
    ts "vast-tools combine finished in \${combine_duration}s (exit code: \$combine_exit)"

    if [ \$combine_exit -ne 0 ]; then
        ts "ERROR: vast-tools combine failed with exit code \$combine_exit"
        ts "--- Directory listing ---"
        ls -la
        ts "--- raw_incl/ contents ---"
        ls -la raw_incl/ 2>/dev/null || echo "  (not found)"
        ts "--- raw_reads/ contents ---"
        ls -la raw_reads/ 2>/dev/null || echo "  (not found)"
        exit \$combine_exit
    fi

    ts "--- Post-combine directory listing ---"
    ls -la *.tab 2>/dev/null || echo "  No .tab files found!"
    ls -la raw_incl/ 2>/dev/null | head -5 || true
    ls -la raw_reads/ 2>/dev/null | head -5 || true

    # ── Find and rename output file ──
    output_file=\$(find . -maxdepth 1 -name "INCLUSION_LEVELS_FULL-${params.species}*.tab" -type f | head -n 1)
    if [ -n "\$output_file" ]; then
        ts "Found output: \$output_file (\$(wc -c < "\$output_file") bytes)"
        cp "\$output_file" "${output_name}_INCLUSION_LEVELS_FULL-${params.species}.tab"
        event_count=\$(tail -n +2 "${output_name}_INCLUSION_LEVELS_FULL-${params.species}.tab" | wc -l)
        sample_count=\$(head -1 "${output_name}_INCLUSION_LEVELS_FULL-${params.species}.tab" | awk -F'\\t' '{print NF}')
        ts "Combined table: \$event_count events, \$sample_count columns"
    else
        ts "ERROR: Output file INCLUSION_LEVELS_FULL-${params.species}*.tab not found"
        ts "All files in working directory:"
        find . -maxdepth 2 -type f | sort
        exit 1
    fi

    ts "=== VAST-tools combine DONE (total: \${combine_duration}s) ==="
    """
}


process compare_groups {
    tag "VAST-tools compare: ${group_a} vs ${group_b}"
    label 'process_medium'
    publishDir "${params.outdir}/compare_results", mode: 'copy'
    container 'andresgordoortiz/vast-tools:latest'

    // Resource requirements
    cpus 4
    memory { 16.GB }
    time { 2.hours }

    input:
    path inclusion_table
    val group_a
    val group_b
    val samples_a
    val samples_b
    val is_paired
    val vastdb_path

    output:
    tuple val(group_a), val(group_b), path("compare_${group_a}_vs_${group_b}"), emit: compare_output
    path "compare_${group_a}_vs_${group_b}/*", emit: compare_results

    script:
    def paired_flag = is_paired ? "--paired" : ""
    def samples_a_str = samples_a.join(',')
    def samples_b_str = samples_b.join(',')
    """
    echo "Running VAST-tools compare: ${group_a} vs ${group_b}"
    echo "Samples in ${group_a}: ${samples_a_str}"
    echo "Samples in ${group_b}: ${samples_b_str}"
    echo "Paired-end data: ${is_paired}"
    echo "Parameters: --min_dPSI ${params.min_dPSI} --min_range ${params.min_range}"

    # Set up VASTDB
    export VASTDB=/usr/local/vast-tools/VASTDB

    # Create output directory
    mkdir -p compare_${group_a}_vs_${group_b}

    # Run vast-tools compare
    # vast-tools compare outputs files to the same directory as the input file
    # So we copy the input file to the output directory and run from there
    cp ${inclusion_table} compare_${group_a}_vs_${group_b}/
    cd compare_${group_a}_vs_${group_b}

    vast-tools compare ${inclusion_table} \\
        -a ${samples_a_str} \\
        -b ${samples_b_str} \\
        --min_dPSI ${params.min_dPSI} \\
        --min_range ${params.min_range} \\
        ${paired_flag} \\
        -sp ${params.species}

    # Remove the copied input file to keep only comparison results
    rm -f ${inclusion_table}

    cd ..

    echo "✓ Comparison complete: ${group_a} vs ${group_b}"

    # List output files
    echo "Output files:"
    ls -la compare_${group_a}_vs_${group_b}/
    """
}

// Process to download reference files (GTF and FASTA) for MATT
process download_matt_references {
    tag "Download references for ${species}"
    label 'process_low'
    storeDir "${params.outdir}/matt_references"

    cpus 2
    memory { 4.GB }
    time { 2.hours }

    input:
    val species

    output:
    tuple val(species), path("*.gtf"), path("*.fa"), emit: references

    shell:
    '''
    #!/bin/bash
    set -e

    # Function to download with retry and validation
    download_file() {
        local url="$1"
        local output="$2"
        local max_retries=3
        local retry=0

        while [ $retry -lt $max_retries ]; do
            echo "Attempting download (try $((retry+1))/${max_retries}): ${url}"
            rm -f "${output}"

            if wget --timeout=60 -q "${url}" -O "${output}" 2>/dev/null; then
                # Check if file exists and has content
                if [ -s "${output}" ]; then
                    # Verify it's actually gzip format
                    if file "${output}" | grep -q "gzip"; then
                        echo "Download successful: ${output}"
                        return 0
                    else
                        echo "Downloaded file is not gzip format, retrying..."
                    fi
                fi
            fi

            # Try curl as fallback
            rm -f "${output}"
            if curl --connect-timeout 60 -sL "${url}" -o "${output}" 2>/dev/null; then
                if [ -s "${output}" ] && file "${output}" | grep -q "gzip"; then
                    echo "Download successful (curl): ${output}"
                    return 0
                fi
            fi

            retry=$((retry + 1))
            [ $retry -lt $max_retries ] && sleep 5
        done

        return 1
    }

    # Species-specific configuration
    # Note: hg19/GRCh37 uses the dedicated GRCh37 Ensembl archive
    case "!{species}" in
        hg19)
            RELEASE="87"
            ASSEMBLY="GRCh37"
            SPECIES_NAME="homo_sapiens"
            SPECIES_CAP="Homo_sapiens"
            # GRCh37 has its own dedicated Ensembl server
            BASE_URL="https://ftp.ensembl.org/pub/grch37/release-${RELEASE}"
            ;;
        hg38)
            RELEASE="110"
            ASSEMBLY="GRCh38"
            SPECIES_NAME="homo_sapiens"
            SPECIES_CAP="Homo_sapiens"
            BASE_URL="https://ftp.ensembl.org/pub/release-${RELEASE}"
            ;;
        mm9)
            RELEASE="67"
            ASSEMBLY="NCBIM37"
            SPECIES_NAME="mus_musculus"
            SPECIES_CAP="Mus_musculus"
            BASE_URL="https://ftp.ensembl.org/pub/release-${RELEASE}"
            ;;
        mm10)
            RELEASE="102"
            ASSEMBLY="GRCm38"
            SPECIES_NAME="mus_musculus"
            SPECIES_CAP="Mus_musculus"
            BASE_URL="https://ftp.ensembl.org/pub/release-${RELEASE}"
            ;;
        rn6)
            RELEASE="104"
            ASSEMBLY="Rnor_6.0"
            SPECIES_NAME="rattus_norvegicus"
            SPECIES_CAP="Rattus_norvegicus"
            BASE_URL="https://ftp.ensembl.org/pub/release-${RELEASE}"
            ;;
        dm6)
            RELEASE="104"
            ASSEMBLY="BDGP6.32"
            SPECIES_NAME="drosophila_melanogaster"
            SPECIES_CAP="Drosophila_melanogaster"
            BASE_URL="https://ftp.ensembl.org/pub/release-${RELEASE}"
            ;;
        *)
            echo "ERROR: Species !{species} not supported"
            exit 1
            ;;
    esac

    GTF_URL="${BASE_URL}/gtf/${SPECIES_NAME}/${SPECIES_CAP}.${ASSEMBLY}.${RELEASE}.gtf.gz"
    FASTA_URL="${BASE_URL}/fasta/${SPECIES_NAME}/dna/${SPECIES_CAP}.${ASSEMBLY}.dna.primary_assembly.fa.gz"
    FASTA_URL_TOPLEVEL="${BASE_URL}/fasta/${SPECIES_NAME}/dna/${SPECIES_CAP}.${ASSEMBLY}.dna.toplevel.fa.gz"

    echo "Downloading GTF and FASTA for !{species} (Ensembl release ${RELEASE})..."
    echo "Base URL: ${BASE_URL}"

    # Download GTF
    echo "Downloading GTF from: ${GTF_URL}"
    if ! download_file "${GTF_URL}" "!{species}.gtf.gz"; then
        echo "ERROR: Failed to download GTF file"
        exit 1
    fi
    gunzip !{species}.gtf.gz

    # Download FASTA - try primary_assembly first, fall back to toplevel
    echo "Downloading FASTA from: ${FASTA_URL}"
    if ! download_file "${FASTA_URL}" "!{species}.fa.gz"; then
        echo "Primary assembly not found, trying toplevel..."
        echo "Downloading FASTA from: ${FASTA_URL_TOPLEVEL}"
        if ! download_file "${FASTA_URL_TOPLEVEL}" "!{species}.fa.gz"; then
            echo "ERROR: Failed to download FASTA file (tried both primary_assembly and toplevel)"
            exit 1
        fi
    fi
    gunzip !{species}.fa.gz

    # Final validation
    if [ ! -s "!{species}.gtf" ]; then
        echo "ERROR: GTF file is empty or missing"
        exit 1
    fi
    if [ ! -s "!{species}.fa" ]; then
        echo "ERROR: FASTA file is empty or missing"
        exit 1
    fi

    echo "Downloaded reference files for !{species}"
    ls -la
    '''
}

// Process to prepare MATT input from vast-tools compare output
process prepare_matt_input {
    tag "Prepare MATT input: ${group_a} vs ${group_b}"
    label 'process_low'
    publishDir "${params.outdir}/matt_analysis/${group_a}_vs_${group_b}/input", mode: 'copy'

    // Resource requirements
    cpus 1
    memory { 2.GB }
    time { 30.min }

    input:
    path compare_dir
    val group_a
    val group_b
    val species

    output:
    tuple val(group_a), val(group_b), path("exons_for_matt.tab"), path("introns_for_matt.tab"), emit: matt_input

    script:
    """
    #!/bin/bash
    set -e

    echo "Preparing MATT input for ${group_a} vs ${group_b}"
    echo "Looking for files in ${compare_dir}"

    # Find the diff file
    DIFF_FILE=\$(find ${compare_dir} -name "*.tab" -type f | head -1)

    if [ -z "\$DIFF_FILE" ]; then
        echo "WARNING: No differential splicing files found. Creating empty MATT input files."
        echo -e "START\\tEND\\tSCAFFOLD\\tSTRAND\\tGENEID\\tEVENT_ID\\tDATASET" > exons_for_matt.tab
        echo -e "START\\tEND\\tSCAFFOLD\\tSTRAND\\tGENEID\\tEVENT_ID\\tDATASET" > introns_for_matt.tab
        exit 0
    fi

    echo "Processing file: \$DIFF_FILE"

    # Create header for output files
    echo -e "START\\tEND\\tSCAFFOLD\\tSTRAND\\tGENEID\\tEVENT_ID\\tDATASET\\tdPSI" > exons_for_matt.tab
    echo -e "START\\tEND\\tSCAFFOLD\\tSTRAND\\tGENEID\\tEVENT_ID\\tDATASET\\tdPSI" > introns_for_matt.tab

    # Process the vast-tools compare output
    # Typical columns: GENE, EVENT, COORD, LENGTH, FullCO, COMPLEX, dPSI, etc.
    # COORD format: chrX:start-end or more complex patterns

    awk -F'\\t' '
    BEGIN { OFS="\\t" }
    NR==1 {
        # Find column indices
        for (i=1; i<=NF; i++) {
            if (\$i == "GENE" || \$i == "GeneID") gene_col = i
            if (\$i == "EVENT" || \$i == "EventType") event_col = i
            if (\$i == "COORD" || \$i == "CO" || \$i == "FullCO") coord_col = i
            if (\$i == "dPSI" || \$i == "deltaPSI") dpsi_col = i
            if (\$i == "STRAND" || \$i == "Strand") strand_col = i
        }
        # Set defaults if not found
        if (!gene_col) gene_col = 1
        if (!event_col) event_col = 2
        if (!coord_col) coord_col = 3
        if (!dpsi_col) dpsi_col = 0
        if (!strand_col) strand_col = 0
        next
    }
    NR>1 {
        gene = \$gene_col
        event = \$event_col
        coord = \$coord_col
        dpsi = (dpsi_col > 0) ? \$dpsi_col : 0
        strand = (strand_col > 0) ? \$strand_col : "+"

        # Parse coordinate: chr:start-end
        if (match(coord, /([^:]+):([0-9]+)-([0-9]+)/, arr)) {
            chrom = arr[1]
            start = arr[2]
            end = arr[3]
        } else {
            next
        }

        # Determine dataset based on dPSI
        if (dpsi == "" || dpsi == "NA") {
            dataset = "ndiff"
        } else if (dpsi + 0 > ${params.min_dPSI}) {
            dataset = "up"
        } else if (dpsi + 0 < -${params.min_dPSI}) {
            dataset = "down"
        } else {
            dataset = "ndiff"
        }

        event_id = gene "_" event "_" chrom ":" start "-" end

        # Output line
        line = start "\\t" end "\\t" chrom "\\t" strand "\\t" gene "\\t" event_id "\\t" dataset "\\t" dpsi

        # Classify as exon or intron
        event_lower = tolower(event)
        if (index(event_lower, "ir") > 0 || index(event_lower, "intron") > 0) {
            print line >> "introns_for_matt.tab"
        } else {
            print line >> "exons_for_matt.tab"
        }
    }
    ' "\$DIFF_FILE"

    echo "Exon events:"
    wc -l exons_for_matt.tab

    echo "Intron events:"
    wc -l introns_for_matt.tab

    echo "✓ MATT input files created"
    """
}
// Process to run MATT cmpr_exons
process run_matt_exons {
    tag "MATT exons: ${group_a} vs ${group_b}"
    label 'process_medium'
    publishDir "${params.outdir}/matt_analysis/${group_a}_vs_${group_b}/exons", mode: 'copy'
    container 'andresgordoortiz/matt-container:latest'

    // Resource requirements
    cpus 4
    memory { 16.GB }
    time { 4.hours }

    input:
    tuple val(group_a), val(group_b), path(exons_tab), path(introns_tab)
    tuple val(species), path(gtf), path(fasta)

    output:
    path "matt_exons_${group_a}_vs_${group_b}/*", emit: exon_results, optional: true

    script:
    def matt_map = ['hg19': 'Hsap', 'hg38': 'Hsap', 'mm9': 'Mmus', 'mm10': 'Mmus', 'rn6': 'Rnor', 'dm6': 'Dmel']
    def matt_species = matt_map.containsKey(species) ? matt_map[species] : 'Hsap'

    """
    echo "Running MATT cmpr_exons for ${group_a} vs ${group_b}"
    echo "Species: ${species} (MATT code: ${matt_species})"
    echo "Intron length parameter: ${params.matt_intron_length}"

    # Check if we have exons to analyze
    exon_count=\$(wc -l < ${exons_tab})
    echo "Number of exon events: \$((exon_count - 1))"

    if [ \$exon_count -le 1 ]; then
        echo "No exon events to analyze. Skipping MATT cmpr_exons."
        mkdir -p matt_exons_${group_a}_vs_${group_b}
        echo "No exon events found for analysis" > matt_exons_${group_a}_vs_${group_b}/no_events.txt
        exit 0
    fi

    # Create output directory
    mkdir -p matt_exons_${group_a}_vs_${group_b}

    # Run MATT cmpr_exons
    # Format: matt cmpr_exons TABLE START END SCAFFOLD STRAND GENEID GTF FASTA SPECIES INTRON_LENGTH DATASET[groups] OUTPUT_DIR
    matt cmpr_exons ${exons_tab} \\
        START END SCAFFOLD STRAND GENEID \\
        ${gtf} ${fasta} ${matt_species} ${params.matt_intron_length} \\
        DATASET[up,down,ndiff] \\
        matt_exons_${group_a}_vs_${group_b} || {
            echo "MATT cmpr_exons failed, but continuing..."
            echo "MATT analysis failed" > matt_exons_${group_a}_vs_${group_b}/error.txt
        }

    echo "✓ MATT exon analysis complete"
    ls -la matt_exons_${group_a}_vs_${group_b}/ || true
    """
}

// Process to run MATT cmpr_introns
process run_matt_introns {
    tag "MATT introns: ${group_a} vs ${group_b}"
    label 'process_medium'
    publishDir "${params.outdir}/matt_analysis/${group_a}_vs_${group_b}/introns", mode: 'copy'
    container 'andresgordoortiz/matt-container:latest'

    // Resource requirements
    cpus 4
    memory { 16.GB }
    time { 4.hours }

    input:
    tuple val(group_a), val(group_b), path(exons_tab), path(introns_tab)
    tuple val(species), path(gtf), path(fasta)

    output:
    path "matt_introns_${group_a}_vs_${group_b}/*", emit: intron_results, optional: true

    script:
    def matt_map = ['hg19': 'Hsap', 'hg38': 'Hsap', 'mm9': 'Mmus', 'mm10': 'Mmus', 'rn6': 'Rnor', 'dm6': 'Dmel']
    def matt_species = matt_map.containsKey(species) ? matt_map[species] : 'Hsap'

    """
    echo "Running MATT cmpr_introns for ${group_a} vs ${group_b}"
    echo "Species: ${species} (MATT code: ${matt_species})"
    echo "Intron length parameter: ${params.matt_intron_length}"

    # Check if we have introns to analyze
    intron_count=\$(wc -l < ${introns_tab})
    echo "Number of intron events: \$((intron_count - 1))"

    if [ \$intron_count -le 1 ]; then
        echo "No intron events to analyze. Skipping MATT cmpr_introns."
        mkdir -p matt_introns_${group_a}_vs_${group_b}
        echo "No intron events found for analysis" > matt_introns_${group_a}_vs_${group_b}/no_events.txt
        exit 0
    fi

    # Create output directory
    mkdir -p matt_introns_${group_a}_vs_${group_b}

    # Run MATT cmpr_introns
    matt cmpr_introns ${introns_tab} \\
        START END SCAFFOLD STRAND GENEID \\
        ${gtf} ${fasta} ${matt_species} ${params.matt_intron_length} \\
        DATASET[up,down,ndiff] \\
        matt_introns_${group_a}_vs_${group_b} || {
            echo "MATT cmpr_introns failed, but continuing..."
            echo "MATT analysis failed" > matt_introns_${group_a}_vs_${group_b}/error.txt
        }

    echo "✓ MATT intron analysis complete"
    ls -la matt_introns_${group_a}_vs_${group_b}/ || true
    """
}

// Process to filter inclusion table and prepare betAS data object
process prepare_betas_data {
    tag "betAS data preparation"
    label 'process_high'
    publishDir "${params.outdir}/betas_analysis", mode: 'copy'
    container 'andresgordoortiz/splicing_analysis_r_crg:v1.5'

    // Resource requirements - betAS filtering is memory intensive
    cpus 2
    memory { 16.GB }
    time { 1.hours }

    input:
    path inclusion_table
    val species
    val filter_n

    output:
    path "betas_filtered_events.RData", emit: betas_data
    path "betas_filtering_summary.txt", emit: summary

    script:
    """
    #!/usr/bin/env Rscript

    # Load required libraries
    library(betAS)

    # Log start
    cat("Starting betAS data preparation\\n")
    cat("Inclusion table:", "${inclusion_table}", "\\n")
    cat("Species:", "${species}", "\\n")
    cat("Filter N:", ${filter_n}, "\\n")

    # Load and filter data using betAS
    cat("Loading splicing data...\\n")
    splicing_data <- getDataset(
        pathTables = "${inclusion_table}",
        tool = "vast-tools"
    )

    cat("Extracting events...\\n")
    splicing_events_raw <- getEvents(splicing_data, tool = "vast-tools")

    cat("Filtering events with N >= ${filter_n}...\\n")
    splicing_events <- filterEvents(splicing_events_raw, N = ${filter_n})

    # Summary statistics
    n_events_raw <- nrow(splicing_events_raw\$PSI)
    n_events_filtered <- nrow(splicing_events\$PSI)
    n_samples <- ncol(splicing_events\$PSI)

    cat("\\n=== Filtering Summary ===\\n")
    cat("Total events before filtering:", n_events_raw, "\\n")
    cat("Events after filtering (N >=", ${filter_n}, "):", n_events_filtered, "\\n")
    cat("Number of samples:", n_samples, "\\n")
    cat("Sample names:", paste(colnames(splicing_events\$PSI), collapse = ", "), "\\n")

    # Save filtered data
    cat("\\nSaving filtered betAS data...\\n")
    save(splicing_events, file = "betas_filtered_events.RData")

    # Write summary file
    summary_text <- paste0(
        "betAS Filtering Summary\\n",
        "=======================\\n",
        "Input file: ${inclusion_table}\\n",
        "Species: ${species}\\n",
        "Filter threshold (N): ${filter_n}\\n",
        "\\n",
        "Results:\\n",
        "- Events before filtering: ", n_events_raw, "\\n",
        "- Events after filtering: ", n_events_filtered, "\\n",
        "- Number of samples: ", n_samples, "\\n",
        "- Sample names: ", paste(colnames(splicing_events\$PSI), collapse = ", "), "\\n"
    )
    writeLines(summary_text, "betas_filtering_summary.txt")

    cat("✓ betAS data preparation complete\\n")
    """
}

// Process to run betAS simulation for a group comparison
process run_betas_simulation {
    tag "betAS simulation: ${group_a} vs ${group_b}"
    label 'process_high'
    publishDir "${params.outdir}/betas_analysis/${group_a}_vs_${group_b}", mode: 'copy'
    container 'andresgordoortiz/splicing_analysis_r_crg:v1.5'

    // Resource requirements - simulations are very memory intensive
    cpus 2
    memory { 60.GB }
    time { 72.hours }

    input:
    path betas_data
    val group_a
    val group_b
    val samples_a
    val samples_b
    val nsim
    val npoints

    output:
    path "betas_${group_a}_vs_${group_b}_results.csv", emit: results
    path "betas_${group_a}_vs_${group_b}_summary.txt", emit: summary

    script:
    def samples_a_r = samples_a.collect { "\"${it}\"" }.join(', ')
    def samples_b_r = samples_b.collect { "\"${it}\"" }.join(', ')
    """
    #!/usr/bin/env Rscript

    # Load required libraries
    library(betAS)

    # Log start
    cat("Starting betAS simulation analysis\\n")
    cat("Comparison: ${group_a} vs ${group_b}\\n")
    cat("Samples in ${group_a}:", "${samples_a.join(', ')}", "\\n")
    cat("Samples in ${group_b}:", "${samples_b.join(', ')}", "\\n")
    cat("Number of simulations:", ${nsim}, "\\n")
    cat("Number of points:", ${npoints}, "\\n")

    # Load filtered betAS data
    cat("\\nLoading filtered betAS data...\\n")
    load("${betas_data}")

    # Get sample names from the data
    all_samples <- colnames(splicing_events\$PSI)
    cat("Available samples in data:", paste(all_samples, collapse = ", "), "\\n")

    # Define sample groups
    samples_group_a <- c(${samples_a_r})
    samples_group_b <- c(${samples_b_r})

    cat("\\nLooking for group A samples:", paste(samples_group_a, collapse = ", "), "\\n")
    cat("Looking for group B samples:", paste(samples_group_b, collapse = ", "), "\\n")

    # Find matching columns (handle potential naming differences)
    find_matching_cols <- function(sample_names, all_cols) {
        matched <- c()
        for (s in sample_names) {
            # Try exact match first
            if (s %in% all_cols) {
                matched <- c(matched, which(all_cols == s))
            } else {
                # Try partial match (sample name might be part of column name)
                partial <- grep(s, all_cols, fixed = TRUE)
                if (length(partial) > 0) {
                    matched <- c(matched, partial[1])
                }
            }
        }
        return(unique(matched))
    }

    colsA <- find_matching_cols(samples_group_a, all_samples)
    colsB <- find_matching_cols(samples_group_b, all_samples)

    cat("\\nMatched column indices for ${group_a}:", paste(colsA, collapse = ", "), "\\n")
    cat("Matched column indices for ${group_b}:", paste(colsB, collapse = ", "), "\\n")
    cat("Matched samples for ${group_a}:", paste(all_samples[colsA], collapse = ", "), "\\n")
    cat("Matched samples for ${group_b}:", paste(all_samples[colsB], collapse = ", "), "\\n")

    # Check if we found samples
    if (length(colsA) == 0 || length(colsB) == 0) {
        cat("\\nERROR: Could not find matching samples for one or both groups\\n")
        cat("Creating empty results file...\\n")

        # Create empty results
        empty_df <- data.frame(
            Error = "No matching samples found",
            Group_A = "${group_a}",
            Group_B = "${group_b}",
            Samples_A_requested = paste(samples_group_a, collapse = ", "),
            Samples_B_requested = paste(samples_group_b, collapse = ", "),
            Available_samples = paste(all_samples, collapse = ", ")
        )
        write.csv(empty_df, "betas_${group_a}_vs_${group_b}_results.csv", row.names = FALSE)

        summary_text <- paste0(
            "betAS Simulation Summary - ERROR\\n",
            "================================\\n",
            "Comparison: ${group_a} vs ${group_b}\\n",
            "Error: Could not find matching samples\\n",
            "Requested ${group_a} samples: ", paste(samples_group_a, collapse = ", "), "\\n",
            "Requested ${group_b} samples: ", paste(samples_group_b, collapse = ", "), "\\n",
            "Available samples: ", paste(all_samples, collapse = ", "), "\\n"
        )
        writeLines(summary_text, "betas_${group_a}_vs_${group_b}_summary.txt")
        quit(status = 0)
    }

    # Set seed for reproducibility
    set.seed(42)
    cat("\\nSeed set to 42 for reproducible simulations\\n")

    # Run betAS simulation
    cat("\\nRunning betAS simulation (this may take a while)...\\n")
    cat("Using", ${nsim}, "simulations and", ${npoints}, "points\\n")

    start_time <- Sys.time()

    splicing_betAS_df <- prepareTableVolcanoFDR(
        psitable = splicing_events\$PSI,
        qualtable = splicing_events\$Qual,
        npoints = ${npoints},
        colsA = colsA,
        colsB = colsB,
        labA = "${group_a}",
        labB = "${group_b}",
        basalColor = "#89C0AE",
        interestColor = "#E69A9C",
        maxDevTable = maxDevSimulationN100,
        nsim = ${nsim},
        seed = TRUE,
        CoverageWeight = FALSE
    )

    end_time <- Sys.time()
    elapsed <- difftime(end_time, start_time, units = "hours")

    cat("\\nSimulation completed in", round(as.numeric(elapsed), 2), "hours\\n")

    # Save results
    cat("\\nSaving results...\\n")
    write.csv(splicing_betAS_df, "betas_${group_a}_vs_${group_b}_results.csv", row.names = TRUE)

    # Calculate summary statistics
    n_events <- nrow(splicing_betAS_df)
    n_significant <- if ("FDR" %in% colnames(splicing_betAS_df)) {
        sum(splicing_betAS_df\$FDR < 0.05, na.rm = TRUE)
    } else if ("pval" %in% colnames(splicing_betAS_df)) {
        sum(splicing_betAS_df\$pval < 0.05, na.rm = TRUE)
    } else {
        NA
    }

    # Write summary
    summary_text <- paste0(
        "betAS Simulation Summary\\n",
        "========================\\n",
        "Comparison: ${group_a} vs ${group_b}\\n",
        "\\n",
        "Parameters:\\n",
        "- Number of simulations: ${nsim}\\n",
        "- Number of points: ${npoints}\\n",
        "- Seed: 42\\n",
        "\\n",
        "Samples:\\n",
        "- ${group_a} (n=", length(colsA), "): ", paste(all_samples[colsA], collapse = ", "), "\\n",
        "- ${group_b} (n=", length(colsB), "): ", paste(all_samples[colsB], collapse = ", "), "\\n",
        "\\n",
        "Results:\\n",
        "- Total events analyzed: ", n_events, "\\n",
        "- Significant events (FDR < 0.05): ", n_significant, "\\n",
        "- Processing time: ", round(as.numeric(elapsed), 2), " hours\\n"
    )
    writeLines(summary_text, "betas_${group_a}_vs_${group_b}_summary.txt")

    cat("\\n✓ betAS simulation complete for ${group_a} vs ${group_b}\\n")
    """
}

// Process to generate pipeline execution summary
process generate_pipeline_summary {
    tag "Pipeline Summary"
    label 'process_low'
    publishDir "${params.outdir}", mode: 'copy'

    // Resource requirements
    cpus 1
    memory { 2.GB }
    time { 30.min }

    input:
    path inclusion_table
    path vast_dirs, stageAs: "vast_out_*"
    val project_name
    val species
    val outdir

    output:
    path "pipeline_summary_${project_name}.txt", emit: summary_txt
    path "pipeline_summary_${project_name}.json", emit: summary_json

    script:
    """
    #!/bin/bash
    set -e

    echo "Generating pipeline execution summary..."

    # Initialize summary file
    SUMMARY_FILE="pipeline_summary_${project_name}.txt"
    JSON_FILE="pipeline_summary_${project_name}.json"

    # Header
    cat > \$SUMMARY_FILE << 'HEADER'
================================================================================
                    VAST-TOOLS SPLICING ANALYSIS PIPELINE SUMMARY
================================================================================
HEADER

    echo "Project: ${project_name}" >> \$SUMMARY_FILE
    echo "Species: ${species}" >> \$SUMMARY_FILE
    echo "Generated: \$(date)" >> \$SUMMARY_FILE
    echo "" >> \$SUMMARY_FILE

    # Initialize JSON
    echo "{" > \$JSON_FILE
    echo "  \"project\": \"${project_name}\"," >> \$JSON_FILE
    echo "  \"species\": \"${species}\"," >> \$JSON_FILE
    echo "  \"timestamp\": \"\$(date -Iseconds)\"," >> \$JSON_FILE

    # ============================================================
    # 1. VAST-TOOLS ALIGN STATISTICS
    # ============================================================
    echo "--------------------------------------------------------------------------------" >> \$SUMMARY_FILE
    echo "1. VAST-TOOLS ALIGNMENT STATISTICS" >> \$SUMMARY_FILE
    echo "--------------------------------------------------------------------------------" >> \$SUMMARY_FILE

    echo "  \"alignment\": {" >> \$JSON_FILE
    echo "    \"samples\": [" >> \$JSON_FILE

    SAMPLE_COUNT=0
    FIRST_SAMPLE=true
    for vast_dir in vast_out_*; do
        if [ -d "\$vast_dir" ]; then
            SAMPLE_COUNT=\$((SAMPLE_COUNT + 1))
            SAMPLE_NAME=\$(basename "\$vast_dir" | sed 's/vast_out_//')

            # Count files in to_combine as a proxy for successful alignment
            TO_COMBINE_FILES=0
            if [ -d "\$vast_dir/to_combine" ]; then
                TO_COMBINE_FILES=\$(find "\$vast_dir/to_combine" -type f | wc -l)
            fi

            echo "  Sample: \$SAMPLE_NAME" >> \$SUMMARY_FILE
            echo "    - Output files in to_combine: \$TO_COMBINE_FILES" >> \$SUMMARY_FILE

            # JSON entry
            if [ "\$FIRST_SAMPLE" = "false" ]; then
                echo "," >> \$JSON_FILE
            fi
            FIRST_SAMPLE=false
            echo -n "      {\"name\": \"\$SAMPLE_NAME\", \"output_files\": \$TO_COMBINE_FILES}" >> \$JSON_FILE
        fi
    done

    echo "" >> \$JSON_FILE
    echo "    ]," >> \$JSON_FILE
    echo "    \"total_samples\": \$SAMPLE_COUNT" >> \$JSON_FILE
    echo "  }," >> \$JSON_FILE

    echo "" >> \$SUMMARY_FILE
    echo "  Total samples aligned: \$SAMPLE_COUNT" >> \$SUMMARY_FILE
    echo "" >> \$SUMMARY_FILE

    # ============================================================
    # 2. VAST-TOOLS COMBINE STATISTICS
    # ============================================================
    echo "--------------------------------------------------------------------------------" >> \$SUMMARY_FILE
    echo "2. VAST-TOOLS COMBINE STATISTICS" >> \$SUMMARY_FILE
    echo "--------------------------------------------------------------------------------" >> \$SUMMARY_FILE

    if [ -f "${inclusion_table}" ]; then
        # Count total events (excluding header)
        TOTAL_EVENTS=\$(tail -n +2 "${inclusion_table}" | wc -l)

        # Count samples (columns after the first few metadata columns)
        HEADER_LINE=\$(head -1 "${inclusion_table}")
        TOTAL_COLUMNS=\$(echo "\$HEADER_LINE" | awk -F'\\t' '{print NF}')

        # Count event types if possible
        EXON_EVENTS=\$(grep -c "EX\\|Alt\\|ALTA\\|ALTD" "${inclusion_table}" 2>/dev/null || echo "N/A")
        IR_EVENTS=\$(grep -c "^IR\\|\\tIR" "${inclusion_table}" 2>/dev/null || echo "N/A")
        MIC_EVENTS=\$(grep -c "MIC" "${inclusion_table}" 2>/dev/null || echo "N/A")

        echo "  Inclusion table: ${inclusion_table}" >> \$SUMMARY_FILE
        echo "  Total splicing events: \$TOTAL_EVENTS" >> \$SUMMARY_FILE
        echo "  Total columns: \$TOTAL_COLUMNS" >> \$SUMMARY_FILE
        echo "" >> \$SUMMARY_FILE
        echo "  Event type estimates:" >> \$SUMMARY_FILE
        echo "    - Exon skipping (EX/Alt): ~\$EXON_EVENTS" >> \$SUMMARY_FILE
        echo "    - Intron retention (IR): ~\$IR_EVENTS" >> \$SUMMARY_FILE
        echo "    - Microexons (MIC): ~\$MIC_EVENTS" >> \$SUMMARY_FILE

        echo "  \"combine\": {" >> \$JSON_FILE
        echo "    \"inclusion_table\": \"${inclusion_table}\"," >> \$JSON_FILE
        echo "    \"total_events\": \$TOTAL_EVENTS," >> \$JSON_FILE
        echo "    \"total_columns\": \$TOTAL_COLUMNS" >> \$JSON_FILE
        echo "  }," >> \$JSON_FILE
    else
        echo "  WARNING: Inclusion table not found!" >> \$SUMMARY_FILE
        echo "  \"combine\": {\"error\": \"inclusion table not found\"}," >> \$JSON_FILE
    fi
    echo "" >> \$SUMMARY_FILE

    # ============================================================
    # 3. VAST-TOOLS COMPARE STATISTICS (scan output directory)
    # ============================================================
    echo "--------------------------------------------------------------------------------" >> \$SUMMARY_FILE
    echo "3. VAST-TOOLS COMPARE STATISTICS (Differential Splicing)" >> \$SUMMARY_FILE
    echo "--------------------------------------------------------------------------------" >> \$SUMMARY_FILE

    echo "  \"compare\": [" >> \$JSON_FILE
    FIRST_COMPARE=true
    COMPARE_COUNT=0

    # Look for compare results in the output directory
    if [ -d "${outdir}/compare_results" ]; then
        for compare_dir in ${outdir}/compare_results/compare_*; do
            if [ -d "\$compare_dir" ]; then
                COMPARE_COUNT=\$((COMPARE_COUNT + 1))
                COMPARISON_NAME=\$(basename "\$compare_dir" | sed 's/compare_//')

                # Find diff files and count significant events
                DIFF_FILE=\$(find "\$compare_dir" -name "*.tab" -type f 2>/dev/null | head -1)
                DIFF_EVENTS="N/A"
                UP_EVENTS="N/A"
                DOWN_EVENTS="N/A"

                if [ -n "\$DIFF_FILE" ] && [ -f "\$DIFF_FILE" ]; then
                    DIFF_EVENTS=\$(tail -n +2 "\$DIFF_FILE" | wc -l)
                    UP_EVENTS=\$(awk -F'\\t' 'NR>1 && \$NF > 0 {count++} END {print count+0}' "\$DIFF_FILE" 2>/dev/null || echo "N/A")
                    DOWN_EVENTS=\$(awk -F'\\t' 'NR>1 && \$NF < 0 {count++} END {print count+0}' "\$DIFF_FILE" 2>/dev/null || echo "N/A")
                fi

                echo "  Comparison: \$COMPARISON_NAME" >> \$SUMMARY_FILE
                echo "    - Differential events: \$DIFF_EVENTS" >> \$SUMMARY_FILE
                echo "    - Up-regulated (dPSI > 0): \$UP_EVENTS" >> \$SUMMARY_FILE
                echo "    - Down-regulated (dPSI < 0): \$DOWN_EVENTS" >> \$SUMMARY_FILE
                echo "" >> \$SUMMARY_FILE

                if [ "\$FIRST_COMPARE" = "false" ]; then
                    echo "," >> \$JSON_FILE
                fi
                FIRST_COMPARE=false
                echo -n "    {\"comparison\": \"\$COMPARISON_NAME\", \"diff_events\": \"\$DIFF_EVENTS\", \"up\": \"\$UP_EVENTS\", \"down\": \"\$DOWN_EVENTS\"}" >> \$JSON_FILE
            fi
        done
    fi

    if [ \$COMPARE_COUNT -eq 0 ]; then
        echo "  No comparisons performed (skipped or no groups defined)" >> \$SUMMARY_FILE
    fi

    echo "" >> \$JSON_FILE
    echo "  ]," >> \$JSON_FILE

    # ============================================================
    # 4. betAS ANALYSIS STATISTICS (scan output directory)
    # ============================================================
    echo "--------------------------------------------------------------------------------" >> \$SUMMARY_FILE
    echo "4. betAS SIMULATION-BASED ANALYSIS STATISTICS" >> \$SUMMARY_FILE
    echo "--------------------------------------------------------------------------------" >> \$SUMMARY_FILE

    echo "  \"betas\": {" >> \$JSON_FILE

    # Look for betAS filtering summary
    FILTER_SUMMARY="${outdir}/betas_analysis/betas_filtering_summary.txt"
    if [ -f "\$FILTER_SUMMARY" ]; then
        echo "  Data Filtering:" >> \$SUMMARY_FILE
        grep -E "Events|samples|threshold" "\$FILTER_SUMMARY" 2>/dev/null | sed 's/^/    /' >> \$SUMMARY_FILE || true
        echo "" >> \$SUMMARY_FILE

        EVENTS_BEFORE=\$(grep "before filtering" "\$FILTER_SUMMARY" 2>/dev/null | grep -oE '[0-9]+' | tail -1 || echo "N/A")
        EVENTS_AFTER=\$(grep "after filtering" "\$FILTER_SUMMARY" 2>/dev/null | grep -oE '[0-9]+' | tail -1 || echo "N/A")
        N_SAMPLES=\$(grep "Number of samples" "\$FILTER_SUMMARY" 2>/dev/null | grep -oE '[0-9]+' | tail -1 || echo "N/A")

        echo "    \"filtering\": {" >> \$JSON_FILE
        echo "      \"events_before\": \"\$EVENTS_BEFORE\"," >> \$JSON_FILE
        echo "      \"events_after\": \"\$EVENTS_AFTER\"," >> \$JSON_FILE
        echo "      \"samples\": \"\$N_SAMPLES\"" >> \$JSON_FILE
        echo "    }," >> \$JSON_FILE
    else
        echo "  No betAS filtering data found (analysis may be skipped or pending)" >> \$SUMMARY_FILE
        echo "    \"filtering\": null," >> \$JSON_FILE
    fi

    echo "    \"simulations\": [" >> \$JSON_FILE
    FIRST_BETAS=true
    BETAS_COUNT=0

    # Look for betAS simulation summaries in output directory
    if [ -d "${outdir}/betas_analysis" ]; then
        for betas_dir in ${outdir}/betas_analysis/*/; do
            if [ -d "\$betas_dir" ]; then
                BETAS_SUMMARY="\${betas_dir}betas_*_summary.txt"
                for summary_file in \$BETAS_SUMMARY; do
                    if [ -f "\$summary_file" ]; then
                        BETAS_COUNT=\$((BETAS_COUNT + 1))
                        COMPARISON_NAME=\$(basename "\$(dirname "\$summary_file")")

                        echo "  Simulation: \$COMPARISON_NAME" >> \$SUMMARY_FILE

                        if grep -q "ERROR" "\$summary_file" 2>/dev/null; then
                            echo "    - Status: ERROR (see details in output file)" >> \$SUMMARY_FILE
                            STATUS="error"
                            TOTAL_EVENTS="N/A"
                            SIG_EVENTS="N/A"
                        else
                            STATUS="completed"
                            TOTAL_EVENTS=\$(grep "Total events" "\$summary_file" 2>/dev/null | grep -oE '[0-9]+' | tail -1 || echo "N/A")
                            SIG_EVENTS=\$(grep "Significant events" "\$summary_file" 2>/dev/null | grep -oE '[0-9]+' | tail -1 || echo "N/A")
                            PROC_TIME=\$(grep "Processing time" "\$summary_file" 2>/dev/null | grep -oE '[0-9.]+' | tail -1 || echo "N/A")

                            echo "    - Status: Completed" >> \$SUMMARY_FILE
                            echo "    - Total events analyzed: \$TOTAL_EVENTS" >> \$SUMMARY_FILE
                            echo "    - Significant events (FDR < 0.05): \$SIG_EVENTS" >> \$SUMMARY_FILE
                            echo "    - Processing time: \$PROC_TIME hours" >> \$SUMMARY_FILE
                        fi
                        echo "" >> \$SUMMARY_FILE

                        if [ "\$FIRST_BETAS" = "false" ]; then
                            echo "," >> \$JSON_FILE
                        fi
                        FIRST_BETAS=false
                        echo -n "      {\"comparison\": \"\$COMPARISON_NAME\", \"status\": \"\$STATUS\", \"total_events\": \"\$TOTAL_EVENTS\", \"significant_events\": \"\$SIG_EVENTS\"}" >> \$JSON_FILE
                    fi
                done
            fi
        done
    fi

    if [ \$BETAS_COUNT -eq 0 ]; then
        echo "  No betAS simulations found (analysis may be skipped or pending)" >> \$SUMMARY_FILE
    fi

    echo "" >> \$JSON_FILE
    echo "    ]" >> \$JSON_FILE
    echo "  }" >> \$JSON_FILE

    # Close JSON
    echo "}" >> \$JSON_FILE

    # ============================================================
    # FOOTER
    # ============================================================
    echo "" >> \$SUMMARY_FILE
    echo "================================================================================" >> \$SUMMARY_FILE
    echo "                              END OF SUMMARY" >> \$SUMMARY_FILE
    echo "================================================================================" >> \$SUMMARY_FILE

    echo ""
    echo "✓ Pipeline summary generated:"
    echo "  - Text: \$SUMMARY_FILE"
    echo "  - JSON: \$JSON_FILE"
    cat \$SUMMARY_FILE
    """
}

process run_rmarkdown_report {
    tag "Generate R analysis report"
    label 'process_high'
    publishDir "${params.outdir}/report", mode: 'copy', pattern: '*.html'
    container 'andresgordoortiz/splicing_analysis_r_crg:v1.5'

    // Resource requirements
    cpus 2
    memory { 30.GB }
    time { 72.hours } // 3 days

    input:
    path inclusion_table
    val rmd_file

    output:
    path "*.html", optional: true, emit: report_html

    when:
    !params.skip_rmarkdown

    script:
    def rmd_basename = file(rmd_file).name
    def rmd_dir = file(rmd_file).parent
    def html_basename = rmd_basename.replace('.Rmd', '.html')

    """
    # Verify the RMarkdown file exists
    if [ ! -f "${rmd_file}" ]; then
        echo "ERROR: RMarkdown file ${rmd_file} not found. Cannot generate the report."
        exit 1
    fi

    # Get absolute path to the RMarkdown file directory
    RMD_DIR="\$(dirname ${rmd_file})"
    echo "RMarkdown directory: \$RMD_DIR"

    # Store the original working directory before changing to RMarkdown directory
    ORIGINAL_WD="\$PWD"
    echo "Original working directory: \$ORIGINAL_WD"

    # Copy the inclusion table to the RMarkdown directory to make it accessible
    cp ${inclusion_table} "\$RMD_DIR/"

    # List files in the RMarkdown directory for debugging
    echo "Files in RMarkdown directory before rendering:"
    ls -la "\$RMD_DIR"

    # Navigate to the RMarkdown directory to render the file in its original context
    cd "\$RMD_DIR"

    # Render the RMarkdown file in its original location
    Rscript -e "rmarkdown::render('${rmd_basename}', output_dir=getwd())"

    # Check if rendering was successful
    if [ -f "${html_basename}" ]; then
        echo "✓ R Markdown report generated successfully in original directory"
        # Copy the HTML report back to the original working directory for publishing
        cp "${html_basename}" "\$ORIGINAL_WD/"
        echo "HTML report copied to working directory for publishing"
    else
        echo "ERROR: R Markdown rendering failed. Check for errors in the RMarkdown file."
        echo "Current directory contents:"
        ls -la
        exit 1
    fi
    """
}

// Function to get VASTDB directory from species code
def getVastdbDirName(species) {
    // Map species names to VASTDB folder names
    def speciesDirectoryMap = [
        'hg19': 'Hsa',
        'hg38': 'Hs2',
        'mm9': 'Mmu',
        'mm10': 'Mm2',
        'rn6': 'Rno',
        'bosTau6': 'Bta',
        'galGal3': 'Gg3',
        'galGal4': 'Gg4',
        'xenTro3': 'Xt1',
        'danRer10': 'Dre',
        'braLan2': 'Bl1',
        'strPur4': 'Spu',
        'dm6': 'Dme',
        'strMar1': 'Sma',
        'ce11': 'Cel',
        'schMed31': 'Sme',
        'nemVec1': 'Nve',
        'araTha10': 'Ath'
    ]

    return speciesDirectoryMap.containsKey(species) ?
           speciesDirectoryMap[species] :
           species
}

// Function to parse sample CSV and extract group information for compare
def parseGroupsFromCsv(sampleCsv) {
    def groups = [:]
    def has_paired = false

    def csvFile = file(sampleCsv)
    def lines = csvFile.readLines()
    def headers = lines[0].split(',').collect { it.trim().toLowerCase() }

    def sampleIdx = headers.indexOf('sample')
    def typeIdx = headers.indexOf('type')
    def groupIdx = headers.indexOf('group')

    if (groupIdx == -1) {
        return [groups: [:], has_paired: false]
    }

    lines.drop(1).each { line ->
        def values = line.split(',').collect { it.trim() }
        if (values.size() > groupIdx && values[groupIdx]) {
            def sample = values[sampleIdx]
            def type = values[typeIdx].toLowerCase()
            def group = values[groupIdx]

            if (!groups.containsKey(group)) {
                groups[group] = []
            }
            groups[group] << sample

            if (type == 'paired') {
                has_paired = true
            }
        }
    }

    return [groups: groups, has_paired: has_paired]
}

// Function to generate pairwise comparisons
def generatePairwiseComparisons(groupsMap, has_paired) {
    def comparisons = []
    def groupNames = groupsMap.keySet().toList()

    for (int i = 0; i < groupNames.size(); i++) {
        for (int j = i + 1; j < groupNames.size(); j++) {
            def group_a = groupNames[i]
            def group_b = groupNames[j]
            comparisons << [
                group_a: group_a,
                group_b: group_b,
                samples_a: groupsMap[group_a],
                samples_b: groupsMap[group_b],
                is_paired: has_paired
            ]
        }
    }

    return comparisons
}

// Define a function to parse the sample CSV and create channels
def parseSamplesCsv(sampleCsv) {
    Channel.fromPath(sampleCsv)
        .splitCsv(header: true, strip: true)
        .map { row ->
            // Check if the required fields are present
            if (row.sample == null || row.fastq_1 == null || row.type == null) {
                log.error "ERROR: Missing required fields in sample CSV. Each row must have 'sample', 'fastq_1', and 'type'."
                exit 1
            }

            // Check if the type is valid (single, paired, technical_replicate)
            if (!['single', 'paired', 'technical_replicate'].contains(row.type.toLowerCase())) {
                log.error "ERROR: Invalid value for 'type' in sample CSV. Must be 'single', 'paired', or 'technical_replicate'."
                exit 1
            }

            def sampleId = row.sample

            // FIX: Handle relative paths by prepending the data directory if needed
            def fastq1Path = row.fastq_1
            if (!fastq1Path.startsWith("/") && !fastq1Path.startsWith("./") && !fastq1Path.startsWith("../")) {
                fastq1Path = "${params.data_dir}/${fastq1Path}"
            }
            def fastq1 = file(fastq1Path, checkIfExists: true)

            // For paired-end reads, check if fastq_2 is provided
            def fastq2 = null
            if (row.fastq_2) {
                def fastq2Path = row.fastq_2
                if (!fastq2Path.startsWith("/") && !fastq2Path.startsWith("./") && !fastq2Path.startsWith("../")) {
                    fastq2Path = "${params.data_dir}/${fastq2Path}"
                }
                fastq2 = file(fastq2Path, checkIfExists: true)
            }

            if (row.type.toLowerCase() == 'paired' && fastq2 == null) {
                log.error "ERROR: Type is set to 'paired' for sample '${sampleId}' but 'fastq_2' is missing."
                exit 1
            }

            // For technical replicates, get all fastq files
            def fastq_files = []
            if (row.type.toLowerCase() == 'technical_replicate') {
                fastq_files.add(fastq1)
                if (fastq2) fastq_files.add(fastq2)

                // Check for additional technical replicates (handle relative paths)
                def additionalFastqs = (3..10).collect { i ->
                    def key = "fastq_${i}"
                    if (row.containsKey(key) && row[key]) {
                        def fastqPath = row[key]
                        if (!fastqPath.startsWith("/") && !fastqPath.startsWith("./") && !fastqPath.startsWith("../")) {
                            fastqPath = "${params.data_dir}/${fastqPath}"
                        }
                        return file(fastqPath, checkIfExists: true)
                    }
                    return null
                }
                additionalFastqs.removeAll([null])
                fastq_files.addAll(additionalFastqs)
            } else if (row.type.toLowerCase() == 'paired') {
                fastq_files = [fastq1, fastq2]
            } else {
                fastq_files = [fastq1]
            }

            // Return a tuple with sample info
            return [sampleId, row.type.toLowerCase(), fastq_files, row.group ?: sampleId]
        }
}

// Define the workflow
workflow {
    // Show help message if help flag is specified
    if (params.help) {
        helpMessage()
        exit 0
    }

    // Validate parameters first
    validateParameters()

    // Execute processes in order with proper logging
    log.info """
    Starting pipeline with parameters:
      - Sample CSV:        ${params.sample_csv}
      - VASTDB path:       ${params.vastdb_path}
      - Data directory:    ${params.data_dir} (MANDATORY)
      - Output directory:  ${params.outdir}
      - Species:           ${params.species}
      - Skip FastQC:       ${params.skip_fastqc}
      - Skip Trimming:     ${params.skip_trimming}
      - Skip FastQC in Trimming: ${params.skip_fastqc_in_trimming}
      - Skip Compare:      ${params.skip_compare}
      - Min dPSI:          ${params.min_dPSI}
      - Min Range:         ${params.min_range}
      - Skip MATT:         ${params.skip_matt}
      - MATT intron length: ${params.matt_intron_length}
      - Skip betAS:        ${params.skip_betas}
      - betAS filter N:    ${params.betas_filter_n}
      - betAS simulations: ${params.betas_nsim}
      - betAS npoints:     ${params.betas_npoints}
    """

    // Prepare VASTDB - now just validates paths
    vastdb_path = prepare_vastdb(params.vastdb_path)

    // Parse sample CSV and create sample channels
    sample_channel = parseSamplesCsv(params.sample_csv)

    // Split samples by type to process differently
    tech_rep_samples = sample_channel.filter { it[1] == 'technical_replicate' }
    paired_samples = sample_channel.filter { it[1] == 'paired' }
    single_samples = sample_channel.filter { it[1] == 'single' }

    // Process samples according to their type
    concatenated_tech_rep = concatenate_technical_replicates(tech_rep_samples)
    prepared_paired = prepare_paired_reads(paired_samples)
    prepared_single = prepare_single_reads(single_samples)

    // Combine all prepared samples
    all_prepared_samples = concatenated_tech_rep
        .mix(prepared_paired)
        .mix(prepared_single)

    // Process flow depends on whether trimming is enabled
    // If trimming is enabled, trim_galore can run FastQC automatically
    // If trimming is disabled, we might run standalone FastQC

    if (!params.skip_trimming) {
        // Run trim_galore on all samples
        trimmed_results = run_trim_galore(all_prepared_samples)

        // Use the trimmed reads for downstream processes
        samples_for_alignment = trimmed_results.trimmed_reads

        // Get FastQC results from trim_galore if fastqc in trimming isn't skipped
        fastqc_results = params.skip_fastqc_in_trimming ?
            Channel.empty() :
            trimmed_results.fastqc_results
    } else {
        // Skip trimming - use the raw reads
        samples_for_alignment = all_prepared_samples

        // Only run standalone FastQC if trimming is skipped (to avoid duplicate FastQC runs)
        fastqc_results = params.skip_fastqc ? Channel.empty() : run_fastqc(all_prepared_samples)
    }

    // Align all samples with VAST-tools
    align_results = align_reads(samples_for_alignment, vastdb_path)

    // Use the named output channels directly
    vast_out_dirs = align_results.vast_out_dir.collect()

    // Get a unique name for the output based on the project name parameter
    output_name = params.project_name

    // Combine all alignment results
    inclusion_table = combine_results(vast_out_dirs, vastdb_path, output_name)

    // Run MultiQC if FastQC was run (either standalone or via trim_galore)
    if ((!params.skip_trimming && !params.skip_fastqc_in_trimming) || (params.skip_trimming && !params.skip_fastqc)) {
        // Use the named output channel directly - no need to map again
        vast_dirs_for_multiqc = align_results.vast_out_dir.collect()

        // Handle optional channels properly
        trim_logs = params.skip_trimming ?
            Channel.empty() :
            run_trim_galore.out.trim_log.collect()

        fastqc_for_multiqc = fastqc_results.collect()

        // Run MultiQC with properly separated inputs
        multiqc_report = run_multiqc(
            fastqc_for_multiqc,
            trim_logs,
            vast_dirs_for_multiqc
        )
    }

    // RMarkdown report generation is deprecated/disabled
    // To re-enable in the future, set params.skip_rmarkdown = false

    // Run betAS simulation-based splicing analysis if not skipped
    // This runs independently after vast-tools combine, doesn't require compare
    if (!params.skip_betas) {
        // Parse groups from CSV file for betAS
        def groupInfoBetas = parseGroupsFromCsv(params.sample_csv)
        def groupsBetas = groupInfoBetas.groups
        def has_paired_betas = groupInfoBetas.has_paired

        if (groupsBetas.size() >= 2) {
            log.info "betAS analysis enabled - preparing filtered data..."
            log.info "Found ${groupsBetas.size()} groups for betAS analysis: ${groupsBetas.keySet().join(', ')}"

            // Generate pairwise comparisons for betAS
            def comparisonsBetas = generatePairwiseComparisons(groupsBetas, has_paired_betas)
            log.info "Will run ${comparisonsBetas.size()} betAS simulation(s)"

            // First, prepare the betAS filtered data object from the inclusion table
            betas_data = prepare_betas_data(
                inclusion_table,
                params.species,
                params.betas_filter_n
            )

            // Create channel for betAS simulations from comparisons
            betas_channel = Channel.from(comparisonsBetas)
                .combine(betas_data.betas_data)

            // Run betAS simulation for each pairwise comparison
            betas_results = run_betas_simulation(
                betas_channel.map { it[1] },           // betas_data RData file
                betas_channel.map { it[0].group_a },
                betas_channel.map { it[0].group_b },
                betas_channel.map { it[0].samples_a },
                betas_channel.map { it[0].samples_b },
                params.betas_nsim,
                params.betas_npoints
            )

            log.info "betAS simulation jobs submitted"
        } else {
            log.info "Skipping betAS analysis: Need at least 2 groups, found ${groupsBetas.size()}"
        }
    } else {
        log.info "Skipping betAS analysis as per user request."
    }

    // Run vast-tools compare for pairwise group comparisons if groups are defined
    if (!params.skip_compare) {
        // Parse groups from CSV file
        def groupInfo = parseGroupsFromCsv(params.sample_csv)
        def groups = groupInfo.groups
        def has_paired = groupInfo.has_paired

        if (groups.size() >= 2) {
            log.info "Found ${groups.size()} groups for comparison: ${groups.keySet().join(', ')}"
            log.info "Paired-end data detected: ${has_paired}"

            // Generate pairwise comparisons
            def comparisons = generatePairwiseComparisons(groups, has_paired)
            log.info "Will run ${comparisons.size()} pairwise comparison(s)"

            // Create channel from comparisons and combine with inclusion table
            // Each item will be: [comparison_map, inclusion_table_path]
            compare_channel = Channel.from(comparisons)
                .combine(inclusion_table)

            // Run compare for each pair
            compare_results = compare_groups(
                compare_channel.map { it[1] },         // inclusion_table (from combine)
                compare_channel.map { it[0].group_a },
                compare_channel.map { it[0].group_b },
                compare_channel.map { it[0].samples_a },
                compare_channel.map { it[0].samples_b },
                compare_channel.map { it[0].is_paired },
                params.vastdb_path
            )

            // Run MATT analysis if not skipped
            if (!params.skip_matt) {
                log.info "MATT analysis enabled - downloading reference files..."

                // Download reference files for MATT (GTF and FASTA)
                matt_refs = download_matt_references(params.species)

                // Prepare MATT input from compare output and run MATT
                // Use the compare_output channel which emits [group_a, group_b, compare_dir]
                matt_input = prepare_matt_input(
                    compare_results.compare_output.map { it[2] },  // compare directory
                    compare_results.compare_output.map { it[0] },  // group_a
                    compare_results.compare_output.map { it[1] },  // group_b
                    params.species
                )

                // Run MATT cmpr_exons
                exon_analysis = run_matt_exons(
                    matt_input,
                    matt_refs
                )

                // Run MATT cmpr_introns
                intron_analysis = run_matt_introns(
                    matt_input,
                    matt_refs
                )

                log.info "MATT analysis jobs submitted"
            } else {
                log.info "Skipping MATT analysis as per user request."
            }

        } else {
            log.info "Skipping vast-tools compare: Need at least 2 groups, found ${groups.size()}"
        }
    } else {
        log.info "Skipping vast-tools compare as per user request."
    }

    // ============================================================
    // PIPELINE SUMMARY - Collect all outputs and generate report
    // ============================================================

    // Collect vast alignment directories (always available)
    vast_dirs_for_summary = align_results.vast_out_dir.collect()

    // Generate pipeline summary - scans output directory for compare and betAS results
    generate_pipeline_summary(
        inclusion_table,
        vast_dirs_for_summary,
        params.project_name,
        params.species,
        params.outdir
    )
}