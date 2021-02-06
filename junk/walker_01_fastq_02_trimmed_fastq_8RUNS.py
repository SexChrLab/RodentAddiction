#######################################
#### Workflow for FastQC, MultiQC, ####
## and trimming of FASTQ read files. ##
######## PART 2/2: SEE NOTES ##########
#######################################

## BEFORE RUNNING THIS SCRIPT ##
# This file, "walker_01_fastq_02_trimmed_fastq_8RUNS.py", must be run after the file "walker_01_fastq_02_trimmed_fastq_1RUN.py".
# In the "1RUN" file, a single sample with only 1 SRA run is analyzed for initial quality.
# In this file, *ALL* samples are trimmed and re-assessed for quality.

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
# References
adapter_dir = config["adapter_dir"] # Adapter
# FASTQs
run_fastq_dir = config["run_fastq_dir"] # Initial FASTQs, separated by run
run_fastqc_dir = config["run_fastqc_dir"] # FastQC/MultiQC for initial FASTQs separated by run
cat_fastq_dir = config["cat_fastq_dir"] # Initial FASTQs, concatenated by sample
cat_fastqc_dir = config["cat_fastqc_dir"] # FastQC/MultiQC for initial FASTQs concatenated by sample
# Trimmed FASTQs
trimmed_fastq_dir = config["trimmed_fastq_dir"] # Trimmed, concatenated FASTQs
trimmed_fastqc_dir = config["trimmed_fastqc_dir"] # FASTQC/MultiQC for trimmed FASTQs

######################
## All Output Files ##
######################

rule all:
	input:
		# FASTQs
		# Initial quality control checks on individual FASTQ files (FastQC), separated by run
		expand(run_fastqc_dir + "walker_{sample}_run1_fastqc.html", sample = config["samples_8runs"]),
		expand(run_fastqc_dir + "walker_{sample}_run2_fastqc.html", sample = config["samples_8runs"]),
#		expand(run_fastqc_dir + "walker_{sample}_run3_fastqc.html", sample = config["samples_8runs"]),
#		expand(run_fastqc_dir + "walker_{sample}_run4_fastqc.html", sample = config["samples_8runs"]),
#		expand(run_fastqc_dir + "walker_{sample}_run5_fastqc.html", sample = config["samples_8runs"]),
#		expand(run_fastqc_dir + "walker_{sample}_run6_fastqc.html", sample = config["samples_8runs"]),
#		expand(run_fastqc_dir + "walker_{sample}_run7_fastqc.html", sample = config["samples_8runs"]),
#		expand(run_fastqc_dir + "walker_{sample}_run8_fastqc.html", sample = config["samples_8runs"]),

		expand(run_fastqc_dir + "walker_{sample}_run1_fastqc.zip", sample = config["samples_8runs"]),
		expand(run_fastqc_dir + "walker_{sample}_run2_fastqc.zip", sample = config["samples_8runs"]),
#		expand(run_fastqc_dir + "walker_{sample}_run3_fastqc.zip", sample = config["samples_8runs"]),
#		expand(run_fastqc_dir + "walker_{sample}_run4_fastqc.zip", sample = config["samples_8runs"]),
#		expand(run_fastqc_dir + "walker_{sample}_run5_fastqc.zip", sample = config["samples_8runs"]),
#		expand(run_fastqc_dir + "walker_{sample}_run6_fastqc.zip", sample = config["samples_8runs"]),
#		expand(run_fastqc_dir + "walker_{sample}_run7_fastqc.zip", sample = config["samples_8runs"]),
#		expand(run_fastqc_dir + "walker_{sample}_run8_fastqc.zip", sample = config["samples_8runs"]),

		# Concatenated FASTQ files
		expand(cat_fastqc_dir + "walker_{sample}_cat.fastq", sample = config["samples_ALL"]),
		# Quality control checks on individual FASTQ files (FastQC), concatenated by sample
		expand(cat_fastqc_dir + "walker_{sample}_cat_fastqc.html", sample = config["samples_ALL"]),
		expand(cat_fastqc_dir + "walker_{sample}_cat_fastqc.zip", sample = config["samples_ALL"]),
		# Quality control check on set of FASTQ files (MultiQC), concatenated by sample
		(cat_fastqc_dir + "walker_multiqc_cat_report.html"),

		# TRIMMED FASTQs
		# Trimmed FASTQ files after quality control and concatenation
		expand(trimmed_fastq_dir + "walker_{sample}_trim.fastq", sample = config["samples_ALL"]),
		# Log files for trimming
		expand(trimmed_fastq_dir + "logfiles/walker_{sample}_trimmomatic.log", sample = config["samples_ALL"]),
		# Quality control checks on trimmed FASTQ files
		expand(trimmed_fastqc_dir + "walker_{sample}_trim_fastqc.html", sample = config["samples_ALL"]),
		expand(trimmed_fastqc_dir + "walker_{sample}_trim_fastqc.zip", sample = config["samples_ALL"]),
		# Quality control check on set of trimmed FASTQ files
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
		RUN2 = run_fastq_dir + "walker_{sample}_run2.fastq"
#		RUN3 = run_fastq_dir + "walker_{sample}_run3.fastq",
#		RUN4 = run_fastq_dir + "walker_{sample}_run4.fastq",
#		RUN5 = run_fastq_dir + "walker_{sample}_run5.fastq",
#		RUN6 = run_fastq_dir + "walker_{sample}_run6.fastq",
#		RUN7 = run_fastq_dir + "walker_{sample}_run7.fastq",
#		RUN8 = run_fastq_dir + "walker_{sample}_run8.fastq"
	output:
		# HTML files; also creates Zip files with same names but .zip extension
		run_fastqc_dir + "walker_{sample}_run1_fastqc.html",
		run_fastqc_dir + "walker_{sample}_run2.fastqc.html"
#		run_fastqc_dir + "walker_{sample}_run3.html",
#		run_fastqc_dir + "walker_{sample}_run4.html",
#		run_fastqc_dir + "walker_{sample}_run5.html",
#		run_fastqc_dir + "walker_{sample}_run6.html",
#		run_fastqc_dir + "walker_{sample}_run7.html",
#		run_fastqc_dir + "walker_{sample}_run8.html"
	params:
		fastqc = fastqc_path,
		run_fastqc_dir = run_fastqc_dir
	run:
		if wildcards.sample in samples_8runs:
			shell("{params.fastqc} -o {params.run_fastqc_dir} {input.RUN1} {input.RUN2}")


#################################
## Quality Control - By Sample ##
#################################

rule concatenate_fastq:
	input:
	output:
		CAT_FQ = cat_fastq_dir + "walker_{sample}_cat.fastq"
	params:
		run_fastq_dir = run_fastq_dir
	run:
		if wildcards.sample in samples_8runs:
			shell("cat {params.run_fastq_dir}*run1.fastq {params.run_fastq_dir}*run2.fastq > {output.CAT_FQ}")
		else:
			shell("cp {params.run_fastq_dir}*onerun.fastq {output.CAT_FQ}")

rule cat_fastqc:
	input:
		FQ = cat_fastq_dir + "walker_{sample}_cat.fastq"
	output:
		FQ_HTML = cat_fastqc_dir + "walker_{sample}_cat_fastqc.html",
		FQ_ZIP = cat_fastqc_dir + "walker_{sample}_cat_fastqc.zip",
	params:
		fastqc = fastqc_path,
		cat_fastqc_dir = cat_fastqc_dir
	shell:
		"""
		{params.fastqc} -o {params.cat_fastqc_dir} {input.FQ}
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

##############
## Trimming ##
##############

rule trimmomatic:
	input:
		FQ = cat_fastq_dir + "walker_{sample}_cat.fastq",
		ADAPTER = adapter_dir + "TruSeq3-SE.fa"
	output:
		TRIMMED_FQ = trimmed_fastq_dir + "walker_{sample}_trim.fastq",
		LOGFILE = trimmed_fastq_dir + "logfiles/walker_{sample}_trimmomatic.log"
	params:
		trimmomatic_jar = trimmomatic_path,
		threads = 4,
		# ILLUMINACLIP
		seed_mismatches = 2,
		palindrome_clip_threshold = 30, # Needs to be specified but has no effect in SE mode
		simple_clip_threshold = 10, # Applies to single-ended reads
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
		ILLUMINACLIP:{input.ADAPTER}:{params.seed_mismatches}:{params.palindrome_clip_threshold} \
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
		{params.fastqc} -o {params.trimmed_fastqc_dir} {input.TRIMMED_FQ}
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