#######################################
#### Workflow for FastQC, MultiQC, ####
## and trimming of FASTQ read files. ##
#######################################

## - Information about project & samples - ##
# Original fastq files come from our laboratory and are separated into 4 runs per sample.
# All sample names (indicated by {sample}) include both the SRA sample ID (SRS#) and the treatment group - e.g. SRS6085261_CI1cue
# Run # is also listed for applicable files.
# Concatenated files can be downloaded from SRA (though FASTQ headers have been stripped).
# SRA Project ID: SRP246331; Bioproject ID: PRJNA604189; Gene Expression Omnibus Series ID: GSE144606

## - Using This Script - ##
# When run in combination with the other Powell Snakemake files and the corresponding JSON file in the same directory, 
# this script creates most of the necessary sub-directories within your starting directory. See powell_directory_structure.txt for 
# the complete directory structure that will result from running these scripts.
#
# To use this script, you need to create the following directory structure:
#
# <your starting directory>
# └── powell
# 	└── 01_fastq
# 		└── by_run
# 			└── <all fastq files>
#
# After creating the above structure, change the following variables in this script:
# 1) All tool paths
# 2) adapter_fasta
# 3) main_dir (your starting directory)

# THINGS TO CHANGE
# 1. Change trimmomatic parameters
# 2. Change to final directory paths

#####################
## Start of Script ##
#####################

import os

## - USERS SHOULD CHANGE THE FOLLOWING VARIABLES - ##
# Tool Paths
fastqc_path = "/home/avannan/miniconda3/envs/rodent_addiction/opt/fastqc-0.11.9/fastqc" # Version 0.11.9
multiqc_path = "/home/avannan/miniconda3/envs/rodent_addiction/bin/multiqc" # Version 1.9
bbduksh_path = "/home/avannan/miniconda3/envs/rodent_addiction/opt/bbmap-38.90-0/bbduk.sh" # Version 38.90
# References & Reference Directories
adapter_fasta = "/home/avannan/miniconda3/envs/rodent_addiction/opt/bbmap-38.90-0/resources/adapters.fa"
# Starting Directory
start_dir = "/data/CEM/wilsonlab/projects/rodent_addiction"
## -- END -- ##

# Config file
configfile: "powell_config.json"

# Directory Variables
# FASTQs
run_fastq_dir = start_dir + config["run_fastq_dir"] # Initial FASTQs, separated by run
run_fastqc_dir = start_dir + config["run_fastqc_dir"] # FastQC for initial FASTQs separated by run + MultiQC by sample
cat_fastq_dir = start_dir + config["cat_fastq_dir"] # Initial FASTQs, concatenated by sample
cat_fastqc_dir = start_dir + config["cat_fastqc_dir"] # FastQC/MultiQC for initial FASTQs concatenated by sample
# Trimmed FASTQs
trimmed_fastq_dir = start_dir + config["trimmed_fastq_dir"] # Trimmed, concatenated FASTQs
trimmed_fastqc_dir = start_dir + config["trimmed_fastqc_dir"] # FASTQC/MultiQC for trimmed FASTQs

######################
## All Output Files ##
######################

rule all:
	input:
		# FASTQs
		# Initial quality control checks on individual FASTQ files (FastQC), separated by run
		# HTML files
		expand(run_fastqc_dir + "{sample}/powell_{sample}_run1_fastqc.html", sample = config["samples_ALL"]),
		expand(run_fastqc_dir + "{sample}/powell_{sample}_run2_fastqc.html", sample = config["samples_ALL"]),
		expand(run_fastqc_dir + "{sample}/powell_{sample}_run3_fastqc.html", sample = config["samples_ALL"]),
		expand(run_fastqc_dir + "{sample}/powell_{sample}_run4_fastqc.html", sample = config["samples_ALL"]),
		# Zip files
		expand(run_fastqc_dir + "{sample}/powell_{sample}_run1_fastqc.zip", sample = config["samples_ALL"]),
		expand(run_fastqc_dir + "{sample}/powell_{sample}_run2_fastqc.zip", sample = config["samples_ALL"]),
		expand(run_fastqc_dir + "{sample}/powell_{sample}_run3_fastqc.zip", sample = config["samples_ALL"]),
		expand(run_fastqc_dir + "{sample}/powell_{sample}_run4_fastqc.zip", sample = config["samples_ALL"]),
		# Sample-level MultiQC files
		expand(run_fastqc_dir + "{sample}/powell_{sample}_multiqc.html", sample = config["samples_ALL"]),

		# Concatenated FASTQ files
		expand(cat_fastq_dir + "powell_{sample}_cat.fastq", sample = config["samples_ALL"]),
		# Quality control checks on individual FASTQ files (FastQC), concatenated by sample
		expand(cat_fastqc_dir + "powell_{sample}_cat_fastqc.html", sample = config["samples_ALL"]),
		expand(cat_fastqc_dir + "powell_{sample}_cat_fastqc.zip", sample = config["samples_ALL"]),
		# Quality control check on set of FASTQ files (MultiQC), concatenated by sample
		(cat_fastqc_dir + "powell_cat_multiqc.html"),

		# TRIMMED FASTQs
		# Trimmed FASTQ files after quality control and concatenation
		expand(trimmed_fastq_dir + "powell_{sample}_trim.fastq", sample = config["samples_ALL"]),
		# Stats files from trimming
		expand(trimmed_fastq_dir + "stats/powell_{sample}_bbduk_stats.txt", sample = config["samples_ALL"]),
		# Quality control checks on individual FASTQ files (trimmed)
		expand(trimmed_fastqc_dir + "powell_{sample}_trim_fastqc.html", sample = config["samples_ALL"]),
		expand(trimmed_fastqc_dir + "powell_{sample}_trim_fastqc.zip", sample = config["samples_ALL"]),
		# Quality control check on set of FASTQ files (MultiQC), trimmed
		(trimmed_fastqc_dir + "powell_trim_multiqc.html")

##############
## 01_fastq ##
##############

##############################
## Quality Reports - By Run ##
##############################

rule run_fastqc:
	input:
		RUN1 = run_fastq_dir + "powell_{sample}_run1.fastq",
		RUN2 = run_fastq_dir + "powell_{sample}_run2.fastq",
		RUN3 = run_fastq_dir + "powell_{sample}_run3.fastq",
		RUN4 = run_fastq_dir + "powell_{sample}_run4.fastq"
	output:
		# HTML files
		run_fastqc_dir + "{sample}/powell_{sample}_run1_fastqc.html",
		run_fastqc_dir + "{sample}/powell_{sample}_run2_fastqc.html",
		run_fastqc_dir + "{sample}/powell_{sample}_run3_fastqc.html",
		run_fastqc_dir + "{sample}/powell_{sample}_run4_fastqc.html",
		# Zip files
		run_fastqc_dir + "{sample}/powell_{sample}_run1_fastqc.zip",
		run_fastqc_dir + "{sample}/powell_{sample}_run2_fastqc.zip",
		run_fastqc_dir + "{sample}/powell_{sample}_run3_fastqc.zip",
		run_fastqc_dir + "{sample}/powell_{sample}_run4_fastqc.zip"
	params:
		fastqc = fastqc_path,
		run_fastqc_dir = run_fastqc_dir
	shell:
		"""
		{params.fastqc} {input.RUN1} {input.RUN2} {input.RUN3} {input.RUN4} -o {params.run_fastqc_dir}
		"""

rule sample_multiqc:
	input: 
		run_fastqc_dir + "{sample}/powell_{sample}_run1_fastqc.zip",
		run_fastqc_dir + "{sample}/powell_{sample}_run2_fastqc.zip",
		run_fastqc_dir + "{sample}/powell_{sample}_run3_fastqc.zip",
		run_fastqc_dir + "{sample}/powell_{sample}_run4_fastqc.zip"
	output:
		MULTIQC_REPORT = run_fastqc_dir + "{sample}/powell_{sample}_multiqc.html"
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
## Quality Reports - By Sample ##
#################################

rule concatenate_fastq:
	input:
		RUN1 = run_fastq_dir + "powell_{sample}_run1.fastq",
		RUN2 = run_fastq_dir + "powell_{sample}_run2.fastq",
		RUN3 = run_fastq_dir + "powell_{sample}_run3.fastq",
		RUN4 = run_fastq_dir + "powell_{sample}_run4.fastq"
	output:
		CAT_FQ = cat_fastq_dir + "powell_{sample}_cat.fastq"
	shell:
		"""
		cat {input.RUN1} {input.RUN2} {input.RUN3} {input.RUN4} > {output.CAT_FQ}
		"""

rule cat_fastqc:
	input:
		CAT_FQ = cat_fastq_dir + "powell_{sample}_cat.fastq"
	output:
		FQ_HTML = cat_fastqc_dir + "powell_{sample}_cat_fastqc.html",
		FQ_ZIP = cat_fastqc_dir + "powell_{sample}_cat_fastqc.zip",
	params:
		fastqc = fastqc_path,
		cat_fastqc_dir = cat_fastqc_dir
	shell:
		"""
		{params.fastqc} {input.CAT_FQ} -o {params.cat_fastqc_dir}
		"""

rule cat_multiqc:
	input:
		expand(cat_fastqc_dir + "powell_{sample}_cat_fastqc.zip", sample = config["samples_ALL"])
	output:
		MULTIQC_REPORT = cat_fastqc_dir + "powell_cat_multiqc.html"
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

##############
## Trimming ##
##############

rule trim_bbduk:
	input:
		FQ = cat_fastq_dir + "powell_{sample}_cat.fastq",
		ADAPTER = adapter_fasta
	output:
		TRIMMED_FQ = trimmed_fastq_dir + "powell_{sample}_trim.fastq",
		STATS = trimmed_fastq_dir + "stats/powell_{sample}_bbduk_stats.txt"
	params:
		bbduksh = bbduksh_path,
		trimq = 15,
		minlen = 50, # Half of read length
		maq = 20
	shell:
		"""
		{params.bbduksh} -Xmx3g in1={input.FQ} out={output.TRIMMED_FQ} \
		ref={input.ADAPTER} ktrim=r k=21 mink=11 hdist=2 stats={output.STATS} \
		trimpolya=10 trimpolyg=10 \
		qtrim=rl trimq={params.trimq} minlen={params.minlen} maq={params.maq}
		"""

#####################
## Quality Reports ##
#####################

rule trimmed_fastqc:
	input:
		TRIMMED_FQ = trimmed_fastq_dir + "powell_{sample}_trim.fastq",
	output:
		TRIMMED_HTML = trimmed_fastqc_dir + "powell_{sample}_trim_fastqc.html",
		TRIMMED_ZIP = trimmed_fastqc_dir + "powell_{sample}_trim_fastqc.zip",
	params:
		fastqc = fastqc_path,
		trimmed_fastqc_dir = trimmed_fastqc_dir
	shell:
		"""
		{params.fastqc} -o {params.trimmed_fastqc_dir} {input.TRIMMED_FQ}
		"""

rule trimmed_multiqc:
	input:
		expand(trimmed_fastqc_dir + "powell_{sample}_trim_fastqc.zip", sample = config["samples_ALL"])
	output:
		TRIMMED_MULTIQC_REPORT = trimmed_fastqc_dir + "powell_trim_multiqc.html"
	message: "Running MultiQC for FastQC reports on -trimmed- FASTQ files."
	params:
		multiqc = multiqc_path,
		trimmed_multiqc_dir = trimmed_fastqc_dir
	shell:
		"""
		export LC_ALL=en_US.UTF-8 && export LANG=en_US.UTF-8 && \
		{params.multiqc} -f {params.trimmed_multiqc_dir} -n {output.TRIMMED_MULTIQC_REPORT} --interactive --verbose
		"""