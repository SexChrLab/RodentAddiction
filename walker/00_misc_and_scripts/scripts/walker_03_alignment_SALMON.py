###########################################
###### Workflow for Salmon alignment ######
###########################################

## - Information about project & samples - ##
# Original fastq files were obtained from the authors.
# SRA Project ID: SRP132477; Bioproject ID: PRJNA433508; Gene Expression Omnibus Series ID: GSE110344
# All sample names (indicated by {sample}) include both the SRA Sample ID (SRS#) and the treatment group - e.g. SRS2926577_C1no.
# Files were initially separated by run (8 runs per sample with the exception of SRS2926581_S1no); run # is listed for applicable files.

## - Using This Script - ##
# See README file in GitHub before running these scripts.
# Remember to change variables under the heading "USERS SHOULD CHANGE THE FOLLOWING VARIABLES"

#####################
## Start of Script ##
#####################

import os

## - USERS SHOULD CHANGE THE FOLLOWING VARIABLES - ##
# Config file
configfile: "walker_config.json"

# Reference Directories
salmon_index = "/scratch/avannan/refs/mouse/mm10/salmon_index"
# Starting Directory
start_dir = "/data/CEM/wilsonlab/projects/rodent_addiction/"
## -- END SECTION FOR CHANGING VARIABLES -- ##


# Directory Variables
# Trimmed FASTQs
trimmed_fastq_dir = os.path.join(start_dir, config["trimmed_fastq_dir"]) # Trimmed FASTQs
# Alignment
salmon_dir = os.path.join(start_dir, config["salmon_dir"]) # Checking library type with Salmon

######################
## All Output Files ##
######################

rule all:
	input:
		# CHECK LIBRARY TYPE (PLUS ALIGN/COUNT) WITH SALMON
		expand(directory(os.path.join(salmon_dir, "walker_{sample}_salmon_quant/")), sample = config["30sal_ALL"])

##################
## 03_alignment ##
##################

####################################
## Check library type with Salmon ##
####################################

rule salmon_quant:
	input:
		TRIMMED_FASTQ = os.path.join(trimmed_fastq_dir, "walker_{sample}_trim.fastq")
	output:
		OUTPUT = directory(os.path.join(salmon_dir, "walker_{sample}_salmon_quant/"))
	params:
		salmon_index = salmon_index,
		libtype = "A", # Automatic detection of library type
		threads = 8
	shell:
		"""
		source ~/miniconda3/etc/profile.d/conda.sh && conda deactivate && conda activate salmon &&
		salmon quant -i {params.salmon_index} -l {params.libtype} -r {input.TRIMMED_FASTQ} -o {output.OUTPUT}
		"""
