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
import regex

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
configfile: "carpenter_config.json"

# Directory Variables
# Misc
misc_dir = start_dir + config["misc_dir"] # Contains lists for MulitQC by sequencing cohort/abstinence group
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
		# Lists of 1abs and 28abs FastQC Zip files
		(misc_dir + "carpenter_1abs_multiqc_list.txt"),
		(misc_dir + "carpenter_1abs_multiqc_list.txt"),
		# Initial quality control check on set of FASTQ files (MultiQC)
		(fastqc_dir + "carpenter_1abs_multiqc.html"),
		(fastqc_dir + "carpenter_28abs_multiqc.html"),

		# TRIMMED FASTQs
		# Trimmed FASTQ files after quality control
		expand(trimmed_fastq_dir + "carpenter_{sample}_trim_fq1.fastq", sample = config["samples_ALL"]),
		expand(trimmed_fastq_dir + "carpenter_{sample}_trim_fq2.fastq", sample = config["samples_ALL"]),
		# Stats files from trimming
		expand(trimmed_fastq_dir + "stats/carpenter_{sample}_bbduk_stats.txt", sample = config["samples_ALL"]),
		# Quality control checks on trimmed FASTQ files
		expand(trimmed_fastqc_dir + "carpenter_{sample}_trim_fq1_fastqc.html", sample = config["samples_ALL"]),
		expand(trimmed_fastqc_dir + "carpenter_{sample}_trim_fq2_fastqc.html", sample = config["samples_ALL"]),
		expand(trimmed_fastqc_dir + "carpenter_{sample}_trim_fq1_fastqc.zip", sample = config["samples_ALL"]),
		expand(trimmed_fastqc_dir + "carpenter_{sample}_trim_fq2_fastqc.zip", sample = config["samples_ALL"]),
		# Lists of 1abs and 28abs FastQC Zip files
		(misc_dir + "carpenter_1abs_trim_multiqc_list.txt"),
		(misc_dir + "carpenter_1abs_trim_multiqc_list.txt"),
		# Quality control check on set of trimmed FASTQ files
		(trimmed_fastqc_dir + "carpenter_1abs_trim_multiqc.html"),
		(trimmed_fastqc_dir + "carpenter_28abs_trim_multiqc.html")

##############
## 01_fastq ##
##############

#####################
## Quality Reports ##
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
		{params.fastqc} {input.FQ1} {input.FQ2} -o {params.fastqc_dir}
		"""

rule initial_multiqc_setup:
	input:
		expand(fastqc_dir + "carpenter_{sample}_fq1_fastqc.zip", sample = config["samples_ALL"]),
		expand(fastqc_dir + "carpenter_{sample}_fq2_fastqc.zip", sample = config["samples_ALL"])
	output:
		LIST_1ABS = misc_dir + "carpenter_1abs_multiqc_list.txt",
		LIST_28ABS = misc_dir + "carpenter_28abs_multiqc_list.txt"
	priority: 1 # Needs to run before initial_multiqc_execute
	run:
		shell("echo -n > {output.LIST_1ABS}")
		shell("echo -n > {output.LIST_28ABS}")
		for i in input:
			if regex.search("1abs", i):
				shell("echo {} >> {{output.LIST_1ABS}}".format(i))
			elif regex.search("28abs", i):
				shell("echo {} >> {{output.LIST_28ABS}}".format(i))

rule initial_multiqc_execute:
	input:
		expand(fastqc_dir + "carpenter_{sample}_fq1_fastqc.zip", sample = config["samples_ALL"]),
		expand(fastqc_dir + "carpenter_{sample}_fq2_fastqc.zip", sample = config["samples_ALL"]),
		LIST_1ABS = misc_dir + "carpenter_1abs_multiqc_list.txt",
		LIST_28ABS = misc_dir + "carpenter_28abs_multiqc_list.txt"
	output:
		MULTIQC_1ABS_REPORT = fastqc_dir + "carpenter_1abs_multiqc.html",
		MULTIQC_28ABS_REPORT = fastqc_dir + "carpenter_28abs_multiqc.html"
	message: "Running MultiQC for FastQC reports on initial FASTQ files, separating by abstinence length."
	params:
		multiqc = multiqc_path
	shell:
		"""
		export LC_ALL=en_US.UTF-8 && export LANG=en_US.UTF-8 && {params.multiqc} -f \
		--file-list {input.LIST_1ABS} -n {output.MULTIQC_1ABS_REPORT} --interactive --verbose
		export LC_ALL=en_US.UTF-8 && export LANG=en_US.UTF-8 && {params.multiqc} -f \
		--file-list {input.LIST_28ABS} -n {output.MULTIQC_28ABS_REPORT} --interactive --verbose
		"""

######################
## 02_trimmed_fastq ##
######################

##############
## Trimming ##
##############

rule trim_bbduk:
	input:
		FQ1 = fastq_dir + "carpenter_{sample}_fq1.fastq",
		FQ2 = fastq_dir + "carpenter_{sample}_fq2.fastq",
		ADAPTER = adapter_fasta
	output:
		TRIMMED_FQ1 = trimmed_fastq_dir + "carpenter_{sample}_trim_fq1.fastq",
		TRIMMED_FQ2 = trimmed_fastq_dir + "carpenter_{sample}_trim_fq2.fastq",
		STATS = trimmed_fastq_dir + "stats/carpenter_{sample}_bbduk_stats.txt"
	params:
		bbduksh = bbduksh_path,
		trimq = 15,
		minlen = lambda wildcards: config[wildcards.sample]["Trim_Min_Length"], # Half of read length
		maq = 20
	shell:
		"""
		{params.bbduksh} -Xmx3g in1={input.FQ1} in2={input.FQ2} out1={output.TRIMMED_FQ1} out2={output.TRIMMED_FQ2} \
		ftm=5 ref={input.ADAPTER} ktrim=r k=21 mink=11 hdist=2 stats={output.STATS} tbo tpe \
		qtrim=rl trimq={params.trimq} minlen={params.minlen} maq={params.maq}
		"""

#####################
## Quality Reports ##
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

rule trimmed_multiqc_setup:
	input:
		expand(trimmed_fastqc_dir + "carpenter_{sample}_trim_fq1_fastqc.zip", sample = config["samples_ALL"]),
		expand(trimmed_fastqc_dir + "carpenter_{sample}_trim_fq2_fastqc.zip", sample = config["samples_ALL"])
	output:
		LIST_TRIM_1ABS = misc_dir + "carpenter_1abs_trim_multiqc_list.txt",
		LIST_TRIM_28ABS = misc_dir + "carpenter_28abs_trim_multiqc_list.txt"
	run:
		shell("echo -n > {output.LIST_TRIM_1ABS}")
		shell("echo -n > {output.LIST_TRIM_28ABS}")
		for i in input:
			if regex.search("1abs", i):
				shell("echo {} >> {{output.LIST_TRIM_1ABS}}".format(i))
			elif regex.search("28abs", i):
				shell("echo {} >> {{output.LIST_TRIM_28ABS}}".format(i))

rule trimmed_multiqc_execute:
	input:
		expand(trimmed_fastqc_dir + "carpenter_{sample}_trim_fq1_fastqc.zip", sample = config["samples_ALL"]),
		expand(trimmed_fastqc_dir + "carpenter_{sample}_trim_fq2_fastqc.zip", sample = config["samples_ALL"]),
		LIST_1ABS_TRIM = misc_dir + "carpenter_1abs_trim_multiqc_list.txt",
		LIST_28ABS_TRIM = misc_dir + "carpenter_28abs_trim_multiqc_list.txt"
	output:
		MULTIQC_1ABS_TRIM_REPORT = trimmed_fastqc_dir + "carpenter_1abs_trim_multiqc.html",
		MULTIQC_28ABS_TRIM_REPORT = trimmed_fastqc_dir + "carpenter_28abs_trim_multiqc.html"
	message: "Running MultiQC for FastQC reports on -trimmed- FASTQ files, separating by abstinence length."
	params:
		multiqc = multiqc_path
	shell:
		"""
		export LC_ALL=en_US.UTF-8 && export LANG=en_US.UTF-8 && {params.multiqc} -f \
		--file-list {input.LIST_1ABS_TRIM} -n {output.MULTIQC_1ABS_TRIM_REPORT} --interactive --verbose
		export LC_ALL=en_US.UTF-8 && export LANG=en_US.UTF-8 && {params.multiqc} -f \
		--file-list {input.LIST_28ABS_TRIM} -n {output.MULTIQC_28ABS_TRIM_REPORT} --interactive --verbose
		"""