#######################################
#### Workflow for FastQC, MultiQC, ####
## and trimming of FASTQ read files. ##
######## PART 1/2: SEE NOTES ##########
#######################################

## BEFORE RUNNING THIS SCRIPT ##
# This file, "walker_01_fastq_02_trimmed_fastq_1RUN.py", must be run before the file "walker_01_fastq_02_trimmed_fastq_8RUNS.py".
# This file assesses the initial, pre-trimming quality of a single sample, SRS2926581_S1no, that has 1 SRA run.
# In "walker_01_fastq_02_trimmed_fastq_8RUNS.py", samples with 8 SRA runs each are analyzed for initial quality
# & *ALL* samples are trimmed and re-assessed for quality.

## - Information about project & samples - ##
# Original fastq files were obtained from the Short Read Archive (SRA).
# SRA Project ID: SRP132477; Bioproject ID: PRJNA433508; Gene Expression Omnibus Series ID: GSE110344
# All sample names (indicated by {sample}) include both the Sequence Read Archive (SRA) Sample ID (SRS#) and the treatment group - e.g. SRS2926577_C1no.
# Files were initially separated by run (8 runs per sample with the exception of SRS2926581_S1no); run # is listed for applicable files.

# THINGS TO CHANGE
# 1. Get appropriate adapter file
# 2. Change trimmomatic parameters
# 3. Change to final directory paths

#####################
## Start of Script ##
#####################

import os

# Config file
configfile: "walker_config.json"

# Tool paths
fastqc_path = "/home/avannan/miniconda3/envs/rodent_addiction/opt/fastqc-0.11.9/fastqc" # Version 0.11.9
multiqc_path = "/home/avannan/miniconda3/envs/rodent_addiction/bin/multiqc" # Version 1.9
trimmomatic_path = "/home/avannan/miniconda3/envs/rodent_addiction/share/trimmomatic-0.39-1/trimmomatic.jar" # Version 0.39-1

# Directory Variables
# FASTQs
run_fastq_dir = config["run_fastq_dir"] # Initial FASTQ file
run_fastqc_dir = config["run_fastqc_dir"] # FastQC/MultiQC for initial FASTQ
cat_fastq_dir = config["cat_fastq_dir"] # FASTQs, concatenated by sample

######################
## All Output Files ##
######################

rule all:
	input:
		# FASTQs
		# Initial quality control check on individual FASTQ file (FastQC)
		expand(run_fastqc_dir + "walker_{sample}_fastqc.html", sample = config["samples_1run"]),
		expand(run_fastqc_dir + "walker_{sample}_fastqc.zip", sample = config["samples_1run"]),

		# Rename & move FASTQ file to match concatenated files for samples with 8 runs
		expand(cat_fastqc_dir + "walker_{sample}_fastqc.zip", sample = config["samples_1run"])

##############
## 01_fastq ##
##############

#####################
## Quality Control ##
#####################

rule run_fastqc:
	input:
		SINGLE = run_fastq_dir + "walker_{sample}.fastq"
	output:
		SINGLE_HTML = run_fastqc_dir + "walker_{sample}_fastqc.html",
		SINGLE_ZIP = run_fastqc_dir + "walker_{sample}_fastqc.zip"
	params:
		fastqc = fastqc_path,
		run_fastqc_dir = run_fastqc_dir
	shell:
		"""
		{params.fastqc} -o {params.run_fastqc_dir} {input.SINGLE}
		"""

###############################
## Homogenize FASTQ Filename ##
###############################

rule rename_single_fastq:
	input:
		SINGLE = run_fastq_dir + "walker_{sample}.fastq"
	output:
		SINGLE_RENAMED = cat_fastq_dir + "walker_{sample}_cat.fastq"
	shell:
		"""
		cp {input.SINGLE} {output.SINGLE_RENAMED}
		"""