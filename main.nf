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
    cpus 2
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

    // Retry configuration - retry once if the process fails
    maxRetries 2
    errorStrategy { task.attempt <= maxRetries ? 'retry' : 'terminate' }

    // Resource requirements
    cpus 8
    memory { 30.GB * task.attempt}
    time { 10.hours }

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

    // HPC-optimized container options
    containerOptions {
        def baseOptions = '--writable-tmpfs --no-home --cleanenv'
        if (workflow.containerEngine == 'singularity') {
            return "${baseOptions} --bind /tmp:/tmp --bind \$PWD:\$PWD"
        } else {
            return '--ulimit stack=unlimited --ulimit memlock=unlimited --shm-size=16g'
        }
    }

    // Increase resources to handle larger datasets
    cpus 8
    memory { 64.GB }
    time { 6.hours }  // Increase time limit to 6 hours

    input:
    path vast_out_dirs, stageAs: "vast_*"
    val vastdb_path
    val output_name

    output:
    path "${output_name}_INCLUSION_LEVELS_FULL-${params.species}.tab", emit: inclusion_table

    script:
    """
    echo "VAST-tools combine attempt ${task.attempt} with ${task.memory} memory"
    echo "Using VASTDB path: ${vastdb_path}"

    # Create a clean working directory for combining
    mkdir -p to_combine

    # Copy files with improved error handling
    echo "Collecting VAST-tools output files..."

    # Debug: List directories
    echo "Available vast_* directories:"
    ls -la vast_* || true

    # Improved file copying - avoid copying files to themselves
    for dir in vast_*; do
        if [ -d "\$dir" ]; then
            echo "Processing directory: \$dir"

            # Look for to_combine subdirectory first
            if [ -d "\$dir/to_combine" ]; then
                echo "Found to_combine directory in \$dir"
                # Copy files with explicit error handling for each type
                for pattern in "*.eej*" "*.exskX" "*.info" "*.IR*" "*.mic*" "*.MULTI*" "*.tab"; do
                    echo "Copying \$pattern files from \$dir/to_combine/"
                    find "\$dir/to_combine" -type f -name "\$pattern" -exec cp {} to_combine/ \\; 2>/dev/null || true
                done
            else
                echo "No to_combine directory in \$dir, searching for files directly"
                # Search the entire directory tree for relevant files
                find "\$dir" -type f \\( -name "*.eej*" -o -name "*.exskX" -o -name "*.info" -o -name "*.IR*" -o -name "*.mic*" -o -name "*.MULTI*" -o -name "*.tab" \\) -exec cp {} to_combine/ \\; 2>/dev/null || true
            fi
        fi
    done

    # Count files and check
    file_count=\$(find to_combine -type f | wc -l)
    echo "Found \$file_count files to combine"
    echo "Sample of files:"
    ls -la to_combine | head -20

    if [ \$file_count -gt 0 ]; then
        # HPC-specific optimizations for Singularity
        echo "Setting environment for HPC Singularity execution..."
        export TMPDIR=\${PWD}/tmp
        mkdir -p \$TMPDIR

        # Set resource limits more conservatively for HPC
        ulimit -s 8192  # 8MB stack instead of unlimited
        ulimit -n 4096  # Increase file descriptors

        # Clean environment for Singularity
        unset LD_LIBRARY_PATH
        unset PERL5LIB

        # Set VASTDB path
        export VASTDB=${vastdb_path}
        echo "VASTDB set to \$VASTDB"
        echo "TMPDIR set to \$TMPDIR"

        # Check if running in Singularity vs Docker
        if [ -n "\${SINGULARITY_CONTAINER:-}" ]; then
            echo "Running in Singularity container: \$SINGULARITY_CONTAINER"

            # Singularity-specific environment
            export PERL5OPT=""  # Clear potentially problematic Perl options
            export PERL_HASH_SEED=0

            # Run with much longer timeout (24 hours) and periodic progress monitoring
            echo "Starting vast-tools combine (Singularity mode)..."

            # Start progress monitoring in background
            (
                while true; do
                    echo "\$(date): Combine operation in progress... Files: \$(find to_combine -type f | wc -l)"
                    if [ -d "tmp" ]; then
                        echo "Temp directory size: \$(du -sh tmp 2>/dev/null || echo 'N/A')"
                    fi
                    echo "Memory usage: \$(free -h 2>/dev/null || echo 'N/A')"
                    sleep 900  # Check every 15 minutes
                done
            ) > combine_progress.log 2>&1 &
            PROGRESS_PID=\$!

            # Use a much longer timeout (86400 seconds = 24 hours)
            timeout 86400 vast-tools combine to_combine/ -sp ${params.species} -o . --verbose > combine.log 2>&1
            combine_exit_code=\$?

            # Kill the progress monitor
            kill \$PROGRESS_PID 2>/dev/null || true
        else
            echo "Running in Docker container"
            # Docker execution with longer timeout
            timeout -k 2h 24h bash -c "vast-tools combine to_combine/ -sp ${params.species} -o . --verbose" 2> combine_error.log
            combine_exit_code=\$?
        fi

        if [ \$combine_exit_code -eq 0 ]; then
            echo "VAST-tools combine completed successfully"

            # Look for output file
            inclusion_file=\$(find . -name "INCLUSION_LEVELS_FULL*${params.species}*" -type f | head -n 1)
            if [ -n "\$inclusion_file" ]; then
                echo "Found inclusion file: \$inclusion_file"
                cp "\$inclusion_file" "${output_name}_INCLUSION_LEVELS_FULL-${params.species}.tab"
                echo "✓ Results successfully combined"
            else
                echo "No inclusion table found, creating placeholder"
                echo "# VAST-tools combine completed but no output file found" > "${output_name}_INCLUSION_LEVELS_FULL-${params.species}.tab"
                echo "# Files processed: \$file_count" >> "${output_name}_INCLUSION_LEVELS_FULL-${params.species}.tab"
                echo "# Created: \$(date)" >> "${output_name}_INCLUSION_LEVELS_FULL-${params.species}.tab"
            fi
        else
            echo "Combine failed with exit code \$combine_exit_code"

            # Comprehensive debug output
            echo "=== DEBUG INFORMATION ==="
            echo "Container environment: \${SINGULARITY_CONTAINER:-docker}"
            echo "Working directory: \$(pwd)"
            echo "Files in to_combine: \$(ls -la to_combine/ | wc -l)"
            echo "Available memory: \$(cat /proc/meminfo | grep MemAvailable || echo 'N/A')"
            echo "Disk space: \$(df -h . || echo 'N/A')"

            # Check logs
            if [ -f "combine.log" ]; then
                echo "=== COMBINE LOG ==="
                tail -100 combine.log
            fi

            if [ -f "combine_error.log" ]; then
                echo "=== ERROR LOG ==="
                cat combine_error.log
            fi

            if [ -f "combine_progress.log" ]; then
                echo "=== PROGRESS LOG ==="
                tail -50 combine_progress.log
            fi

            echo "Creating error report file"
            echo "# VAST-tools combine failed with exit code \$combine_exit_code" > "${output_name}_INCLUSION_LEVELS_FULL-${params.species}.tab"
            echo "# Container: \${SINGULARITY_CONTAINER:-docker}" >> "${output_name}_INCLUSION_LEVELS_FULL-${params.species}.tab"
            echo "# Input files: \$file_count" >> "${output_name}_INCLUSION_LEVELS_FULL-${params.species}.tab"
            echo "# Error occurred at: \$(date)" >> "${output_name}_INCLUSION_LEVELS_FULL-${params.species}.tab"

            # Include log content if available
            if [ -f "combine.log" ]; then
                echo "# Log content:" >> "${output_name}_INCLUSION_LEVELS_FULL-${params.species}.tab"
                tail -50 combine.log >> "${output_name}_INCLUSION_LEVELS_FULL-${params.species}.tab" 2>/dev/null || true
            fi

            # If timeout occurred, make that clear
            if [ \$combine_exit_code -eq 124 ]; then
                echo "# ERROR: Process timed out after 24 hours" >> "${output_name}_INCLUSION_LEVELS_FULL-${params.species}.tab"
                echo "# Consider splitting your dataset or increasing the timeout further" >> "${output_name}_INCLUSION_LEVELS_FULL-${params.species}.tab"
            fi
        fi
    else
        echo "No files found to combine, creating empty output file"
        echo "# No VAST-tools output files found to combine" > "ewing_splicing_PRJNA407215_INCLUSION_LEVELS_FULL-hg19.tab"
        echo "# Created: \$(date)" >> "ewing_splicing_PRJNA407215_INCLUSION_LEVELS_FULL-hg19.tab"
    fi
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
    # Note: vast-tools compare outputs to current directory by default
    # We run from inside the output directory to capture all output files
    cd compare_${group_a}_vs_${group_b}

    vast-tools compare ../${inclusion_table} \\
        -a ${samples_a_str} \\
        -b ${samples_b_str} \\
        --min_dPSI ${params.min_dPSI} \\
        --min_range ${params.min_range} \\
        ${paired_flag} \\
        -sp ${params.species}

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

    # Species-specific configuration
    case "!{species}" in
        hg19)
            RELEASE="75"
            ASSEMBLY="GRCh37"
            SPECIES_NAME="homo_sapiens"
            SPECIES_CAP="Homo_sapiens"
            ;;
        hg38)
            RELEASE="110"
            ASSEMBLY="GRCh38"
            SPECIES_NAME="homo_sapiens"
            SPECIES_CAP="Homo_sapiens"
            ;;
        mm9)
            RELEASE="67"
            ASSEMBLY="NCBIM37"
            SPECIES_NAME="mus_musculus"
            SPECIES_CAP="Mus_musculus"
            ;;
        mm10)
            RELEASE="102"
            ASSEMBLY="GRCm38"
            SPECIES_NAME="mus_musculus"
            SPECIES_CAP="Mus_musculus"
            ;;
        rn6)
            RELEASE="104"
            ASSEMBLY="Rnor_6.0"
            SPECIES_NAME="rattus_norvegicus"
            SPECIES_CAP="Rattus_norvegicus"
            ;;
        dm6)
            RELEASE="104"
            ASSEMBLY="BDGP6.32"
            SPECIES_NAME="drosophila_melanogaster"
            SPECIES_CAP="Drosophila_melanogaster"
            ;;
        *)
            echo "ERROR: Species !{species} not supported"
            exit 1
            ;;
    esac

    GTF_URL="https://ftp.ensembl.org/pub/release-${RELEASE}/gtf/${SPECIES_NAME}/${SPECIES_CAP}.${ASSEMBLY}.${RELEASE}.gtf.gz"
    FASTA_URL="https://ftp.ensembl.org/pub/release-${RELEASE}/fasta/${SPECIES_NAME}/dna/${SPECIES_CAP}.${ASSEMBLY}.dna.primary_assembly.fa.gz"

    echo "Downloading GTF and FASTA for !{species} (Ensembl release ${RELEASE})..."

    # Download GTF
    echo "Downloading GTF from: ${GTF_URL}"
    wget -q "${GTF_URL}" -O !{species}.gtf.gz || curl -sL "${GTF_URL}" -o !{species}.gtf.gz
    gunzip !{species}.gtf.gz

    # Download FASTA - try primary_assembly first, fall back to toplevel
    echo "Downloading FASTA..."
    if ! wget -q "${FASTA_URL}" -O !{species}.fa.gz 2>/dev/null; then
        echo "Primary assembly not found, trying toplevel..."
        FASTA_TOPLEVEL=$(echo "${FASTA_URL}" | sed 's/primary_assembly/toplevel/')
        wget -q "${FASTA_TOPLEVEL}" -O !{species}.fa.gz || curl -sL "${FASTA_TOPLEVEL}" -o !{species}.fa.gz
    fi
    gunzip !{species}.fa.gz

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
}