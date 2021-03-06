##########################################
## Workflow for counting with StringTie ##
##########################################

## - Information about project & samples - ##
# Original fastq files come from our laboratory and are separated into 4 runs per sample.
# All sample names (indicated by {sample}) include both the SRA sample ID (SRS#) and the treatment group - e.g. SRS6085261_CI1cue
# Run # is also listed for applicable files.
# Concatenated files can be downloaded from SRA (though FASTQ headers have been stripped).
# SRA Project ID: SRP246331; Bioproject ID: PRJNA604189; Gene Expression Omnibus Series ID: GSE144606

## - Using This Script - ##
# See powell_01_fastq_02_trimmed_fastq.py for information on initial directory setup before running this script.
# Remember to change the following variables:
# 1) All tool paths
# 2) genome_fasta
# 3) annotation_gtf
# 3) main_dir (your starting directory)

# THINGS TO CHANGE
# 1. Change StringTie parameters
# 2. Change to final directory paths

#####################
## Start of Script ##
#####################

import os

## - USERS SHOULD CHANGE THE FOLLOWING VARIABLES - ##
# Tool Paths
stringtie_path = "/home/avannan/miniconda3/envs/rodent_addiction/bin/stringtie" # Version 2.1.4
# References & Reference Directories
genome_fasta = "/home/avannan/rodent_addiction/refs/rat/Rn6_genome.fa"
annotation_gtf = "/home/avannan/rodent_addiction/refs/rat/Rn6_annotation.gtf"
# Starting Directory
start_dir = "/data/CEM/wilsonlab/projects/rodent_addiction"
## -- END -- ##

# Config file
configfile: "powell_config.json"

# Directory Variables
# Processing
hisat2_sorted_bams_dir = start_dir + config["hisat2_sorted_bams_dir"] # Sorted HISAT2 BAM files
hisat2_sorted_bams_stats_dir = start_dir + config["hisat2_sorted_bams_stats_dir"] # Stats for sorted HISAT2 BAMs
# Counting
hisat2_stringtie_dir = start_dir + config["hisat2_stringtie_dir"] # HISAT2/StringTie directory; all StringTie files go in sample-specific subfolders

######################
## All Output Files ##
######################

rule all:
	input:
		# PROCESSING AFTER HISAT2 ALIGNMENT
		# HISAT2 paired, trimmed, sorted, BAM files and their indexes
		expand(hisat2_sorted_bams_dir + "powell_{sample}_hisat2_trim_sort.bam", sample = config["samples_ALL"]),
		expand(hisat2_sorted_bams_dir + "powell_{sample}_hisat2_trim_sort.bam.bai", sample = config["samples_ALL"]),
		# Stats for paired, trimmed, sorted, HISAT2 BAM files
		expand(hisat2_sorted_bams_stats_dir + "powell_{sample}_hisat2_trim_sort_stats.txt", sample = config["samples_ALL"]),

		# COUNTING WITH STRINGTIE AFTER HISAT2 ALIGNMENT
		# Tab files are also created in the same directory as the assembled transcripts, but not listed here
		expand(hisat2_stringtie_assem_tran_dir + "firstpass/powell_{sample}_hisat2_stringtie_assembled_transcripts_firstpass.gtf", sample = config["samples_ALL"]),
		(hisat2_stringtie_dir + "powell_hisat2_stringtie_list.txt"),
		(hisat2_stringtie_dir + "powell_hisat2_stringtie_merge.gtf"),
		expand(hisat2_stringtie_assem_tran_dir + "{sample}/powell_{sample}_hisat2_stringtie_assembled_transcripts_final.gtf", sample = config["samples_ALL"]),
		expand(hisat2_stringtie_cov_tran_dir + "powell_{sample}_hisat2_stringtie_covered_transcripts_final.gtf", sample = config["samples_ALL"]),
		expand(hisat2_stringtie_gene_dir + "powell_{sample}_hisat2_stringtie_genes_final.txt", sample = config["samples_ALL"])

###############
## 05_counts ##
###############

######################################################
## Counting with StringTie (after HISAT2 alignment) ##
######################################################

rule hisat2_stringtie_first_pass:
	input:
		BAM = hisat2_sorted_bams_dir + "powell_{sample}_hisat2_trim_sort.bam",
		GUIDE_GTF = annotation_gtf
	output: hisat2_stringtie_assem_tran_dir + "firstpass/powell_{sample}_hisat2_stringtie_assembled_transcripts_firstpass.gtf"
	params:
		stringtie = stringtie_path,
		p = 4 # Processing threads
	shell:
		"""
		{params.stringtie} {input.BAM} -o {output} -p {params.p} -G {input.GUIDE_GTF}
		"""

rule hisat2_stringtie_merge_list:
	input: expand(hisat2_stringtie_assem_tran_dir + "firstpass/powell_{sample}_hisat2_stringtie_assembled_transcripts_firstpass.gtf", sample = config["samples_ALL"])
	output: 
		STRINGTIE_MERGE_LIST = hisat2_stringtie_dir + "powell_hisat2_stringtie_list.txt"
	run:
		shell("echo -n > {output.STRINGTIE_MERGE_LIST}")
		for i in input:
			shell("echo {} >> {{output.STRINGTIE_MERGE_LIST}}".format(i))

rule hisat2_stringtie_merge:
	input: 
		STRINGTIE_MERGE_LIST = hisat2_stringtie_dir + "powell_hisat2_stringtie_list.txt",
		GUIDE_GTF = annotation_gtf
	output: 
		MERGED_GTF = hisat2_stringtie_dir + "powell_hisat2_stringtie_merge.gtf"
	params:
		stringtie = stringtie_path,
		p = 4 # Processing threads
	shell:
		"""
		{params.stringtie} --merge {input.STRINGTIE_MERGE_LIST} -o {output.MERGED_GTF} -p {params.p} -G {input.GUIDE_GTF}
		"""

rule hisat2_stringtie_final:
	input: 
		BAM = hisat2_sorted_bams_dir + "powell_{sample}_hisat2_trim_sort.bam",
		MERGED_GTF = hisat2_stringtie_dir + "powell_hisat2_stringtie_merge.gtf"
	output:
		ASSEMBLED_TRANSCRIPTS = hisat2_stringtie_assem_tran_dir + "{sample}/powell_{sample}_hisat2_stringtie_assembled_transcripts_final.gtf",
		COVERED_TRANSCRIPTS = hisat2_stringtie_cov_tran_dir + "powell_{sample}_hisat2_stringtie_covered_transcripts_final.gtf",
		GENE_ABUNDANCES = hisat2_stringtie_gene_dir + "powell_{sample}_hisat2_stringtie_genes_final.txt",
	message: "Counting transcripts from {input.BAM} using StringTie merge mode."
	params:
		stringtie = stringtie_path,
		f = 0.15, # Minimum isoform abundance of predicted transcripts, as a fraction of most abundant transcript assembled at a locus; Default = 0.01
		p = 8, # Processing threads; Default = 1
		m = 200, # Minimum length of predicted transcripts; Default = 50
		M = 0.95, # Maximum fraction of multiple-location-mapped reads allowed at a locus; Default = 0.95
		g = 50 # Minimum locus gap separation value; Default = 50
	shell:
		"""
		{params.stringtie} {input.BAM} -e -B -f {params.f} -p {params.p} -m {params.m} -M {params.M} -g {params.g} \
		-G {input.MERGED_GTF} -o {output.ASSEMBLED_TRANSCRIPTS} -C {output.COVERED_TRANSCRIPTS} -A {output.GENE_ABUNDANCES}
		"""