#######################################
#### Workflow for FastQC, MultiQC, ####
## and trimming of FASTQ read files. ##
#######################################

## - Information about project & samples - ##
# Original fastq files were obtained from the authors.
# SRA Project ID: SRP132477; Bioproject ID: PRJNA433508; Gene Expression Omnibus Series ID: GSE110344
# All sample names (indicated by {sample}) include both the SRA Sample ID (SRS#) and the treatment group - e.g. SRS2926577_C1no.
# Files were initially separated by run (8 runs per sample with the exception of SRS2926581_S1no); run # is listed for applicable files.

## - Using This Script - ##
# When run in combination with the other Walker Snakemake files and the corresponding JSON file in the same directory, 
# this script creates most of the necessary sub-directories within your starting directory. See walker_directory_structure.txt for 
# the complete directory structure that will result from running these scripts.
#
# To use this script, you need to create the following directory structure:
#
# <your starting directory>
# └── walker
# 	└── 01_fastq
# 		└── by_run
# 			└── <all fastq files except below>
# 		└── by_run
# 			└── walker_SRS2926581_S1no_cat.fastq
#
# After creating the above structure, change the following variables in this script:
# 1) All tool paths
# 2) adapter_fasta
# 3) main_dir (your starting directory)

# THINGS TO CHANGE
# 1. Get appropriate adapter file
# 2. Change trimmomatic parameters
# 3. Change to final directory paths
# 4. Fix concatenation of runs
# 5. Add all runs

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
configfile: "walker_config.json"

# Directory Variables
# FASTQs
run_fastq_dir = start_dir + config["run_fastq_dir"] # Initial FASTQs, separated by run
run_fastqc_dir = start_dir + config["run_fastqc_dir"] # FastQC for initial FASTQs separated by run + MultiQC by sample
cat_fastq_dir = start_dir + config["cat_fastq_dir"] # Initial FASTQs, concatenated by sample
cat_fastqc_dir = start_dir + config["cat_fastqc_dir"] # FastQC/MultiQC for initial FASTQs concatenated by sample
# Trimmed FASTQs
trimmed_fastq_dir = start_dir + config["trimmed_fastq_dir"] # Trimmed FASTQs
trimmed_fastqc_dir = start_dir + config["trimmed_fastqc_dir"] # FASTQC/MultiQC for trimmed FASTQs

######################
## All Output Files ##
######################

rule all:
	input:
		# FASTQs
		# Initial quality control checks on individual FASTQ files (FastQC), separated by run
		# HTML files
		expand(run_fastqc_dir + "{sample}/walker_{sample}_run1_fastqc.html", sample = config["samples_8runs"]),
		expand(run_fastqc_dir + "{sample}/walker_{sample}_run2_fastqc.html", sample = config["samples_8runs"]),
		expand(run_fastqc_dir + "{sample}/walker_{sample}_run3_fastqc.html", sample = config["samples_8runs"]),
		expand(run_fastqc_dir + "{sample}/walker_{sample}_run4_fastqc.html", sample = config["samples_8runs"]),
		expand(run_fastqc_dir + "{sample}/walker_{sample}_run5_fastqc.html", sample = config["samples_8runs"]),
		expand(run_fastqc_dir + "{sample}/walker_{sample}_run6_fastqc.html", sample = config["samples_8runs"]),
		expand(run_fastqc_dir + "{sample}/walker_{sample}_run7_fastqc.html", sample = config["samples_8runs"]),
		expand(run_fastqc_dir + "{sample}/walker_{sample}_run8_fastqc.html", sample = config["samples_8runs"]),
		# Zip files
		expand(run_fastqc_dir + "{sample}/walker_{sample}_run1_fastqc.zip", sample = config["samples_8runs"]),
		expand(run_fastqc_dir + "{sample}/walker_{sample}_run2_fastqc.zip", sample = config["samples_8runs"]),
		expand(run_fastqc_dir + "{sample}/walker_{sample}_run3_fastqc.zip", sample = config["samples_8runs"]),
		expand(run_fastqc_dir + "{sample}/walker_{sample}_run4_fastqc.zip", sample = config["samples_8runs"]),
		expand(run_fastqc_dir + "{sample}/walker_{sample}_run5_fastqc.zip", sample = config["samples_8runs"]),
		expand(run_fastqc_dir + "{sample}/walker_{sample}_run6_fastqc.zip", sample = config["samples_8runs"]),
		expand(run_fastqc_dir + "{sample}/walker_{sample}_run7_fastqc.zip", sample = config["samples_8runs"]),
		expand(run_fastqc_dir + "{sample}/walker_{sample}_run8_fastqc.zip", sample = config["samples_8runs"]),
		# Sample-level MultiQC files
		expand(run_fastqc_dir + "{sample}/walker_{sample}_multiqc.html", sample = config["samples_8runs"]),

		# Concatenated FASTQ files
		expand(cat_fastq_dir + "walker_{sample}_cat.fastq", sample = config["samples_ALL"]),
		# Quality control checks on individual FASTQ files (FastQC), concatenated by sample
		expand(cat_fastqc_dir + "walker_{sample}_cat_fastqc.html", sample = config["samples_ALL"]),
		expand(cat_fastqc_dir + "walker_{sample}_cat_fastqc.zip", sample = config["samples_ALL"]),
		# Quality control check on set of FASTQ files (MultiQC), concatenated by sample
		(cat_fastqc_dir + "walker_cat_multiqc.html"),

		# TRIMMED FASTQs
		# Trimmed FASTQ files after quality control and concatenation
		expand(trimmed_fastq_dir + "walker_{sample}_trim.fastq", sample = config["samples_ALL"]),
		# Log files for trimming
		expand(trimmed_fastq_dir + "logfiles/walker_{sample}_trimmomatic.log", sample = config["samples_ALL"]),
		# Quality control checks on individual FASTQ files (trimmed)
		expand(trimmed_fastqc_dir + "walker_{sample}_trim_fastqc.html", sample = config["samples_ALL"]),
		expand(trimmed_fastqc_dir + "walker_{sample}_trim_fastqc.zip", sample = config["samples_ALL"]),
		# Quality control check on set of FASTQ files (MultiQC), trimmed
		(trimmed_fastqc_dir + "walker_trim_multiqc_report.html")

##############
## 01_fastq ##
##############

##############################
## Quality Control - By Run ##
##############################

rule run_fastqc:
	input:
		RUN1 = run_fastq_dir + "walker_{sample}_run1.fastq",
		RUN2 = run_fastq_dir + "walker_{sample}_run2.fastq",
		RUN3 = run_fastq_dir + "walker_{sample}_run3.fastq",
		RUN4 = run_fastq_dir + "walker_{sample}_run4.fastq",
		RUN5 = run_fastq_dir + "walker_{sample}_run5.fastq",
		RUN6 = run_fastq_dir + "walker_{sample}_run6.fastq",
		RUN7 = run_fastq_dir + "walker_{sample}_run7.fastq",
		RUN8 = run_fastq_dir + "walker_{sample}_run8.fastq"
	output:
		# HTML files
		run_fastqc_dir + "{sample}/walker_{sample}_run1_fastqc.html",
		run_fastqc_dir + "{sample}/walker_{sample}_run2_fastqc.html",
		run_fastqc_dir + "{sample}/walker_{sample}_run3_fastqc.html",
		run_fastqc_dir + "{sample}/walker_{sample}_run4_fastqc.html",
		run_fastqc_dir + "{sample}/walker_{sample}_run5_fastqc.html",
		run_fastqc_dir + "{sample}/walker_{sample}_run6_fastqc.html",
		run_fastqc_dir + "{sample}/walker_{sample}_run7_fastqc.html",
		run_fastqc_dir + "{sample}/walker_{sample}_run8_fastqc.html",
		# Zip files
		run_fastqc_dir + "{sample}/walker_{sample}_run1_fastqc.zip",
		run_fastqc_dir + "{sample}/walker_{sample}_run2_fastqc.zip",
		run_fastqc_dir + "{sample}/walker_{sample}_run3_fastqc.zip",
		run_fastqc_dir + "{sample}/walker_{sample}_run4_fastqc.zip",
		run_fastqc_dir + "{sample}/walker_{sample}_run5_fastqc.zip",
		run_fastqc_dir + "{sample}/walker_{sample}_run6_fastqc.zip",
		run_fastqc_dir + "{sample}/walker_{sample}_run7_fastqc.zip",
		run_fastqc_dir + "{sample}/walker_{sample}_run8_fastqc.zip"
	params:
		fastqc = fastqc_path,
		run_fastqc_dir = run_fastqc_dir
	shell:
		"""
		{params.fastqc} {input.RUN1} {input.RUN2} {input.RUN3} {input.RUN4} \
		{input.RUN5} {input.RUN6} {input.RUN7} {input.RUN8} -o {params.run_fastqc_dir}
		"""

rule sample_multiqc:
	input: 
		run_fastqc_dir + "{sample}/walker_{sample}_run1_fastqc.zip",
		run_fastqc_dir + "{sample}/walker_{sample}_run2_fastqc.zip",
		run_fastqc_dir + "{sample}/walker_{sample}_run3_fastqc.zip",
		run_fastqc_dir + "{sample}/walker_{sample}_run4_fastqc.zip",
		run_fastqc_dir + "{sample}/walker_{sample}_run5_fastqc.zip",
		run_fastqc_dir + "{sample}/walker_{sample}_run6_fastqc.zip",
		run_fastqc_dir + "{sample}/walker_{sample}_run7_fastqc.zip",
		run_fastqc_dir + "{sample}/walker_{sample}_run8_fastqc.zip"
	output:
		MULTIQC_REPORT = run_fastqc_dir + "{sample}/walker_{sample}_multiqc.html"
	message: "Running MultiQC for FastQC reports on initial FASTQ files."
	params:
		multiqc = multiqc_path,
		multiqc_dir = run_fastqc_dir + "{sample}"
	shell:
		"""
		export LC_ALL=en_US.UTF-8 && export LANG=en_US.UTF-8 && \
		{params.multiqc} -f {params.multiqc_dir} -n {output.MULTIQC_REPORT} --interactive --verbose
		"""


#################################
## Quality Control - By Sample ##
#################################

rule concatenate_fastq:
	input:
		RUN1 = run_fastq_dir + "walker_{sample}_run1.fastq",
		RUN2 = run_fastq_dir + "walker_{sample}_run2.fastq",
		RUN3 = run_fastq_dir + "walker_{sample}_run3.fastq",
		RUN4 = run_fastq_dir + "walker_{sample}_run4.fastq",
		RUN5 = run_fastq_dir + "walker_{sample}_run5.fastq",
		RUN6 = run_fastq_dir + "walker_{sample}_run6.fastq",
		RUN7 = run_fastq_dir + "walker_{sample}_run7.fastq",
		RUN8 = run_fastq_dir + "walker_{sample}_run8.fastq"
	output:
		CAT_FQ = cat_fastq_dir + "walker_{sample}_cat.fastq"
	shell:
		"""
		cat {input.RUN1} {input.RUN2} {input.RUN3} {input.RUN4} \
		{input.RUN5} {input.RUN6} {input.RUN7} {input.RUN8} > {output.CAT_FQ}
		"""

rule cat_fastqc:
	input:
		CAT_FQ = cat_fastq_dir + "walker_{sample}_cat.fastq"
	output:
		FQ_HTML = cat_fastqc_dir + "walker_{sample}_cat_fastqc.html",
		FQ_ZIP = cat_fastqc_dir + "walker_{sample}_cat_fastqc.zip",
	params:
		fastqc = fastqc_path,
		cat_fastqc_dir = cat_fastqc_dir
	shell:
		"""
		{params.fastqc} -o {params.cat_fastqc_dir} {input.CAT_FQ}
		"""

rule cat_multiqc:
	input:
		expand(cat_fastqc_dir + "walker_{sample}_cat_fastqc.zip", sample = config["samples_ALL"])
	output:
		MULTIQC_REPORT = cat_fastqc_dir + "walker_cat_multiqc.html"
	message: "Running MultiQC for FastQC reports on -concatenated- FASTQ files."
	params:
		multiqc = multiqc_path,
		cat_fastqc_dir = cat_fastqc_dir
	shell:
		"""
		export LC_ALL=en_US.UTF-8 && export LANG=en_US.UTF-8 && \
		{params.multiqc} -f {params.cat_fastqc_dir} -n {output.MULTIQC_REPORT} --interactive --verbose
		"""

######################
## 02_trimmed_fastq ##
######################

rule trimmomatic:
	input:
		FQ = cat_fastq_dir + "walker_{sample}_cat.fastq",
		ADAPTER = adapter_fasta
	output:
		TRIMMED_FQ = trimmed_fastq_dir + "walker_{sample}_trim.fastq",
		LOGFILE = trimmed_fastq_dir + "logfiles/walker_{sample}_trimmomatic.log"
	params:
		trimmomatic_jar = trimmomatic_path,
		threads = 4,
		# LEADING & TRAILING
		leading = 3, # May not need
		trailing = 3, # May not need
		# SLIDINGWINDOW
		winsize = 4,
		winqual = 30,
		minlen = 50 # Half of read length
	shell:
		"""
		java -jar {params.trimmomatic_jar} SE -phred33 -threads {params.threads} -trimlog {output.LOGFILE} \
		{input.FQ} {output.TRIMMED_FQ} \
		LEADING:{params.leading} TRAILING:{params.trailing} \
		SLIDINGWINDOW:{params.winsize}:{params.winqual} \
		MINLEN:{params.minlen}
		"""

#####################
## Quality Control ##
#####################

rule trimmed_fastqc:
	input:
		TRIMMED_FQ = trimmed_fastq_dir + "walker_{sample}_trim.fastq",
	output:
		TRIMMED_HTML = trimmed_fastqc_dir + "walker_{sample}_trim_fastqc.html",
		TRIMMED_ZIP = trimmed_fastqc_dir + "walker_{sample}_trim_fastqc.zip",
	params:
		fastqc = fastqc_path,
		trimmed_fastqc_dir = trimmed_fastqc_dir
	shell:
		"""
		{params.fastqc} {input.TRIMMED_FQ} -o {params.trimmed_fastqc_dir}
		"""

rule trimmed_multiqc:
	input:
		expand(trimmed_fastqc_dir + "walker_{sample}_trim_fastqc.zip", sample = config["samples_ALL"])
	output:
		TRIMMED_MULTIQC_REPORT = trimmed_fastqc_dir + "walker_trim_multiqc_report.html"
	message: "Running MultiQC for FastQC reports on -trimmed- FASTQ files."
	params:
		multiqc = multiqc_path,
		trimmed_multiqc_dir = trimmed_fastqc_dir
	shell:
		"""
		export LC_ALL=en_US.UTF-8 && export LANG=en_US.UTF-8 && \
		{params.multiqc} -f {params.trimmed_multiqc_dir} -n {output.TRIMMED_MULTIQC_REPORT} --interactive --verbose
		"""