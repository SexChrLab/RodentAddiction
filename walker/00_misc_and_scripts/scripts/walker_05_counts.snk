##########################################
## Workflow for counting with StringTie ##
##########################################

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

# Tool Paths
stringtie_path = "/home/avannan/miniconda3/envs/rodent_addiction/bin/stringtie" # Version 2.1.4
# References
annotation_gtf = "/scratch/avannan/refs/mouse/GRCm39/Mus_musculus.GRCm39.105.gtf"
# Starting Directory
start_dir = "/data/CEM/wilsonlab/projects/rodent_addiction/"
# Script Paths
prepDE_path = "/data/CEM/wilsonlab/projects/rodent_addiction/prepDE.py3" # Last updated for StringTie version 2.0.4 release
## -- END SECTION FOR CHANGING VARIABLES -- ##

# Directory Variables
# Miscellaneous
misc_dir = os.path.join(start_dir, config["misc_dir"]) # Contains miscellaneous files, including text file for converting to raw counts
# Processing
hisat2_sorted_bams_dir = os.path.join(start_dir, config["hisat2_sorted_bams_dir"]) # Sorted HISAT2 BAM files
# Counting
hisat2_stringtie_dir = os.path.join(start_dir, config["hisat2_stringtie_dir"]) # Main HISAT2/StringTie directory
hisat2_stringtie_count_tran_dir = os.path.join(start_dir, config["hisat2_stringtie_count_tran_dir"]) # Transcript counts
hisat2_stringtie_cov_tran_dir = os.path.join(start_dir, config["hisat2_stringtie_cov_tran_dir"]) # Covered transcripts
hisat2_stringtie_gene_dir = os.path.join(start_dir, config["hisat2_stringtie_gene_dir"]) # Genes
hisat2_stringtie_final_count_dir = os.path.join(start_dir, config["hisat2_stringtie_final_count_dir"]) # Final raw counts from prepDE.py3
hisat2_featurecounts_gene_dir = os.path.join(start_dir, config["hisat2_featurecounts_gene_dir"]) # featureCounts

######################
## All Output Files ##
######################

rule all:
	input:
		# COUNTING WITH FEATURECOUNTS AFTER HISAT2 ALIGNMENT
		expand(os.path.join(hisat2_featurecounts_gene_dir, "walker_{sample}_hisat2_featurecounts_genes_s2.txt"), sample = config["30sal_ALL"])

###############
## 05_counts ##
###############

##########################################################
## Counting with featureCounts (after HISAT2 alignment) ##
##########################################################

rule hisat2_featurecounts_gene_count_per_sample_s2:
	priority: 1
	input:
		BAM = os.path.join(hisat2_sorted_bams_dir, "walker_{sample}_hisat2_trim_sort.bam"),
		GTF = annotation_gtf
	output:
		GENE_COUNTS = os.path.join(hisat2_featurecounts_gene_dir, "walker_{sample}_hisat2_featurecounts_genes_s2.txt")
	message: "Counting {wildcards.sample} reads (single end, reverse-stranded) using featureCounts for primary alignments only."
	params:
		s = 2, # Strandness; reverse
		threads = 5
	shell:
		"""
		featureCounts -T {params.threads} --primary -s {params.s} -t gene -g gene_id -a {input.GTF} -o {output.GENE_COUNTS} {input.BAM}
		"""

