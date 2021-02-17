#######################################
#### Workflow for FastQC, MultiQC, ####
## and trimming of FASTQ read files. ##
#######################################

## - Information about project & samples - ##
# Original fastq files were obtained from the Short Read Archive (SRA).
# All sample names (indicated by {sample}) include both the Sequence Read Archive (SRA) Sample ID (SRS#) and the treatment group - e.g. SRS5770694_C1abs.
# SRA Project ID: SRP234876; Bioproject ID: PRJNA593775; Gene Expression Omnibus Series ID: GSE141520

## - Using This Script - ##
# When run in combination with the other Carpenter Snakemake files and the corresponding JSON file in the same directory, 
# this script creates most of the necessary sub-directories within your starting directory. See carpenter_directory_structure.txt for 
# the complete directory structure that will result from running these scripts.
#
# To use this script, you need to create the following directory structure:
#
# <your starting directory>
# └── carpenter
# 	└── 01_fastq
# 		└── <all fastq files>
#
# After creating the above structure, change the following variables in this script:
# 1) All tool paths
# 2) adapter_fasta
# 3) main_dir (your starting directory)

# THINGS TO CHANGE
# 1. Get appropriate adapter file
# 2. Change trimmomatic parameters
# 3. Change to final directory paths

#####################
## Start of Script ##
#####################

import os

## - USERS SHOULD CHANGE THE FOLLOWING VARIABLES - ##
# Tool Paths
fastqc_path = "/home/avannan/miniconda3/envs/rodent_addiction/opt/fastqc-0.11.9/fastqc" # Version 0.11.9
multiqc_path = "/home/avannan/miniconda3/envs/rodent_addiction/bin/multiqc" # Version 1.9
trimmomatic_path = "/home/avannan/miniconda3/envs/rodent_addiction/share/trimmomatic-0.39-1/trimmomatic.jar" # Version 0.39-1
# References & Reference Directories
adapter_fasta = "/home/avannan/miniconda3/envs/rodent_addiction/share/trimmomatic/adapters/TruSeq3-PE-2.fa"
# Starting Directory
start_dir = "/data/CEM/wilsonlab/projects/rodent_addiction"
## -- END -- ##

# Config file
configfile: "carpenter_config.json"

# Directory Variables
# FASTQs
fastq_dir = start_dir + config["fastq_dir"] # Initial FASTQs
fastqc_dir = start_dir + config["fastqc_dir"] # FastQC/MultiQC for initial FASTQs
# Trimmed FASTQs
trimmed_fastq_dir = start_dir + config["trimmed_fastq_dir"] # Trimmed FASTQs
trimmed_fastqc_dir = start_dir + config["trimmed_fastqc_dir"] # FASTQC/MultiQC for trimmed FASTQs

######################
## All Output Files ##
######################

rule all:
	input:
		# FASTQs
		# Initial quality control checks on individual FASTQ files (FastQC)
		# HTML files
		expand(fastqc_dir + "carpenter_{sample}_fq1_fastqc.html", sample = config["samples_ALL"]),
		expand(fastqc_dir + "carpenter_{sample}_fq2_fastqc.html", sample = config["samples_ALL"]),
		# Zip files
		expand(fastqc_dir + "carpenter_{sample}_fq1_fastqc.zip", sample = config["samples_ALL"]),
		expand(fastqc_dir + "carpenter_{sample}_fq2_fastqc.zip", sample = config["samples_ALL"]),
		# Initial quality control check on set of FASTQ files (MultiQC)
		(fastqc_dir + "carpenter_multiqc.html"),

		# TRIMMED FASTQs
		# Trimmed FASTQ files after quality control
		expand(trimmed_fastq_dir + "carpenter_{sample}_trim_fq1.fastq", sample = config["samples_ALL"]),
		expand(trimmed_fastq_dir + "carpenter_{sample}_trim_fq2.fastq", sample = config["samples_ALL"]),
		# Log files for trimming
		expand(trimmed_fastq_dir + "logfiles/carpenter_{sample}_trimmomatic.log", sample = config["samples_ALL"]),
		# Quality control checks on trimmed FASTQ files
		expand(trimmed_fastqc_dir + "carpenter_{sample}_trim_fq1_fastqc.html", sample = config["samples_ALL"]),
		expand(trimmed_fastqc_dir + "carpenter_{sample}_trim_fq2_fastqc.html", sample = config["samples_ALL"]),
		expand(trimmed_fastqc_dir + "carpenter_{sample}_trim_fq1_fastqc.zip", sample = config["samples_ALL"]),
		expand(trimmed_fastqc_dir + "carpenter_{sample}_trim_fq2_fastqc.zip", sample = config["samples_ALL"]),
		# Quality control check on set of trimmed FASTQ files
		(trimmed_fastqc_dir + "carpenter_trim_multiqc.html")

##############
## 01_fastq ##
##############

#####################
## Quality Control ##
#####################

rule initial_fastqc:
	input:
		FQ1 = fastq_dir + "carpenter_{sample}_fq1.fastq",
		FQ2 = fastq_dir + "carpenter_{sample}_fq2.fastq"
	output:
		FQ1_HTML = fastqc_dir + "carpenter_{sample}_fq1_fastqc.html",
		FQ2_HTML = fastqc_dir + "carpenter_{sample}_fq2_fastqc.html",
		FQ1_ZIP = fastqc_dir + "carpenter_{sample}_fq1_fastqc.zip",
		FQ2_ZIP = fastqc_dir + "carpenter_{sample}_fq2_fastqc.zip"
	params:
		fastqc = fastqc_path,
		fastqc_dir = fastqc_dir
	shell:
		"""
		{params.fastqc} {input.FQ1} {input.FQ2} -o {params.run_fastqc_dir}
		"""

rule initial_multiqc:
	input: 
		expand(fastqc_dir + "carpenter_{sample}_fq1_fastqc.zip", sample = config["samples_ALL"]),
		expand(fastqc_dir + "carpenter_{sample}_fq2_fastqc.zip", sample = config["samples_ALL"])
	output:
		MULTIQC_REPORT = fastqc_dir + "carpenter_multiqc.html"
	message: "Running MultiQC for FastQC reports on initial FASTQ files."
	params:
		multiqc = multiqc_path,
		fastqc_dir = fastqc_dir
	shell:
		"""
		export LC_ALL=en_US.UTF-8 && export LANG=en_US.UTF-8 && \
		{params.multiqc} -f {params.fastqc_dir} -n {output.MULTIQC_REPORT} --interactive --verbose
		"""

######################
## 02_trimmed_fastq ##
######################

##############
## Trimming ##
##############

rule trimmomatic:
	input:
		FQ1 = fastq_dir + "carpenter_{sample}_fq1.fastq",
		FQ2 = fastq_dir + "carpenter_{sample}_fq2.fastq",
		ADAPTER = adapter_fasta
	output:
		TRIMMED_FQ1 = trimmed_fastq_dir + "carpenter_{sample}_trim_fq1.fastq",
		TRIMMED_FQ2 = trimmed_fastq_dir + "carpenter_{sample}_trim_fq2.fastq",
		TRIMMED_UNPAIRED_FQ1 = trimmed_fastq_dir + "carpenter_{sample}_trim_unpair_fq1.fastq",
		TRIMMED_UNPAIRED_FQ2 = trimmed_fastq_dir + "carpenter_{sample}_trim_unpair_fq2.fastq",
		LOGFILE = trimmed_fastq_dir + "logfiles/carpenter_{sample}_trimmomatic.log"
	params:
		# CHANGE TRIMMOMATIC PARAMETERS. CHECK PHRED.
		trimmomatic_jar = trimmomatic_path,
		threads = 4,
		seed_mismatches = 2,
		palindrome_clip_threshold = 30,
		simple_clip_threshold = 10,
		leading = 3, # May not need
		trailing = 3, # May not need
		winsize = 4,
		winqual = 30,
		minlen = 63 # Half of read length (125 bp)
	shell:
		"""
		java -jar {params.trimmomatic_jar} PE -phred33 -threads {params.threads} -trimlog {output.LOGFILE} \
		{input.FQ1} {input.FQ2} \
		{output.TRIMMED_FQ1} {output.TRIMMED_UNPAIRED_FQ1} \
		{output.TRIMMED_FQ2} {output.TRIMMED_UNPAIRED_FQ2} \
		ILLUMINACLIP:{input.ADAPTER}:{params.seed_mismatches}:{params.palindrome_clip_threshold}:{params.simple_clip_threshold} \
		LEADING:{params.leading} TRAILING:{params.trailing} \
		SLIDINGWINDOW:{params.winsize}:{params.winqual} \
		MINLEN:{params.minlen}
		"""

#####################
## Quality Control ##
#####################

rule trimmed_fastqc:
	input:
		TRIMMED_FQ1 = trimmed_fastq_dir + "carpenter_{sample}_trim_fq1.fastq",
		TRIMMED_FQ2 = trimmed_fastq_dir + "carpenter_{sample}_trim_fq2.fastq"
	output:
		TRIMMED_FQ1_HTML = trimmed_fastqc_dir + "carpenter_{sample}_trim_fq1_fastqc.html",
		TRIMMED_FQ2_HTML = trimmed_fastqc_dir + "carpenter_{sample}_trim_fq2_fastqc.html",
		TRIMMED_FQ1_ZIP = trimmed_fastqc_dir + "carpenter_{sample}_trim_fq1_fastqc.zip",
		TRIMMED_FQ2_ZIP = trimmed_fastqc_dir + "carpenter_{sample}_trim_fq2_fastqc.zip"
	params:
		fastqc = fastqc_path,
		trimmed_fastqc_dir = trimmed_fastqc_dir
	shell:
		"""
		{params.fastqc} {input.TRIMMED_FQ1} {input.TRIMMED_FQ2} -o {params.trimmed_fastqc_dir}
		"""

rule trimmed_multiqc:
	input:
		expand(trimmed_fastqc_dir + "carpenter_{sample}_trim_fq1_fastqc.zip", sample = config["samples_ALL"]),
		expand(trimmed_fastqc_dir + "carpenter_{sample}_trim_fq2_fastqc.zip", sample = config["samples_ALL"])
	output:
		TRIMMED_MULTIQC_REPORT = trimmed_fastqc_dir + "carpenter_trim_multiqc.html"
	message: "Running MultiQC for FastQC reports on -trimmed- FASTQ files."
	params:
		multiqc = multiqc_path,
		trimmed_multiqc_dir = trimmed_fastqc_dir
	shell:
		"""
		export LC_ALL=en_US.UTF-8 && export LANG=en_US.UTF-8 && \
		{params.multiqc} -f {params.trimmed_multiqc_dir} -n {output.TRIMMED_MULTIQC_REPORT} --interactive --verbose
		"""