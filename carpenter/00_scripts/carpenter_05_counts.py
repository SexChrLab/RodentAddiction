##########################################
## Workflow for counting with StringTie ##
##########################################

## - Information about project & samples - ##
# Original fastq files were obtained from the Short Read Archive (SRA).
# All sample names (indicated by {sample}) include both the Sequence Read Archive (SRA) Sample ID (SRS#) and the treatment group - e.g. SRS5770694_C1abs.
# SRA Project ID: SRP234876; Bioproject ID: PRJNA593775; Gene Expression Omnibus Series ID: GSE141520

## - Using This Script - ##
# See README file in GitHub before running these scripts.
# Remember to change variables under the heading "USERS SHOULD CHANGE THE FOLLOWING VARIABLES"

#####################
## Start of Script ##
#####################

import os

## - USERS SHOULD CHANGE THE FOLLOWING VARIABLES - ##
# Config file
configfile: "carpenter_config.json"

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
		# COUNTING WITH STRINGTIE AFTER HISAT2 ALIGNMENT
		# Tab files are also created in the same directory as the counted transcripts, but are not listed here
		expand(os.path.join(hisat2_stringtie_count_tran_dir, "{sample}/carpenter_{sample}_hisat2_stringtie_count_transcripts.gtf"), sample = config["28abs_ALL"]),
		expand(os.path.join(hisat2_stringtie_cov_tran_dir, "carpenter_{sample}_hisat2_stringtie_covered_transcripts.gtf"), sample = config["28abs_ALL"]),
		expand(os.path.join(hisat2_stringtie_gene_dir, "carpenter_{sample}_hisat2_stringtie_genes.txt"), sample = config["28abs_ALL"]),
		os.path.join(hisat2_stringtie_final_count_dir, "carpenter_gene_count_matrix_28abs.csv"),
		os.path.join(hisat2_stringtie_final_count_dir, "carpenter_transcript_count_matrix_28abs.csv"),
		# COUNTING WITH FEATURECOUNTS AFTER HISAT2 ALIGNMENT
		expand(os.path.join(hisat2_featurecounts_gene_dir, "carpenter_{sample}_hisat2_featurecounts_genes_s2.txt"), sample = config["28abs_ALL"])

###############
## 05_counts ##
###############

######################################################
## Counting with StringTie (after HISAT2 alignment) ##
######################################################

rule hisat2_stringtie_count:
	priority: 1
	input:
		BAM = os.path.join(hisat2_sorted_bams_dir, "carpenter_{sample}_hisat2_pair_trim_sort.bam"),
		GUIDE_GTF = annotation_gtf
	output:
		TRANSCRIPT_GTF = os.path.join(hisat2_stringtie_count_tran_dir, "{sample}/carpenter_{sample}_hisat2_stringtie_count_transcripts.gtf"),
		COVERED_TRANSCRIPTS = os.path.join(hisat2_stringtie_cov_tran_dir, "carpenter_{sample}_hisat2_stringtie_covered_transcripts.gtf"),
		GENE_ABUNDANCES = os.path.join(hisat2_stringtie_gene_dir, "carpenter_{sample}_hisat2_stringtie_genes.txt")
	params:
		stringtie = stringtie_path,
		p = 8, # Processing threads; Default = 1
		M = 0.95, # Maximum fraction of multiple-location-mapped reads allowed at a locus; Default = 0.95
		g = 50 # Minimum locus gap separation value; Default = 50
	shell:
		"""
		{params.stringtie} {input.BAM} -e -B -p {params.p} -M {params.M} -g {params.g} -G {input.GUIDE_GTF} \
		-o {output.TRANSCRIPT_GTF} -C {output.COVERED_TRANSCRIPTS} -A {output.GENE_ABUNDANCES}
		"""

rule run_prepDE_28abs:
	priority: 10
	input: 
		PREP_TXT = os.path.join(misc_dir, "carpenter_prepDE_28abs.txt") # Contains sample IDs and assembled transcript GTF locations
	output:
		GENE_MATRIX = os.path.join(hisat2_stringtie_final_count_dir, "carpenter_gene_count_matrix_28abs.csv"),
		TRANSCRIPT_MATRIX = os.path.join(hisat2_stringtie_final_count_dir, "carpenter_transcript_count_matrix_28abs.csv")
	params:
		gtf_dir = hisat2_stringtie_count_tran_dir,
		prepDE3 = prepDE_path,
		length = 150 # Average read length
	shell:
		"""
		mkdir -p {params.gtf_dir}
		python3 {params.prepDE3} -i {params.gtf_dir} -l {params.length} -g {output.GENE_MATRIX} -t {output.TRANSCRIPT_MATRIX} -v
		"""

##########################################################
## Counting with featureCounts (after HISAT2 alignment) ##
##########################################################

rule hisat2_featurecounts_gene_count_per_sample_s2:
	priority: 1
	input:
		BAM = os.path.join(hisat2_sorted_bams_dir, "carpenter_{sample}_hisat2_pair_trim_sort.bam"),
		GTF = annotation_gtf
	output:
		GENE_COUNTS = os.path.join(hisat2_featurecounts_gene_dir, "carpenter_{sample}_hisat2_featurecounts_genes_s2.txt")
	message: "Counting {wildcards.sample} reads (paired, reverse-stranded) using featureCounts for primary alignments only."
	params:
		s = 2, # Strandness
		threads = 5
	shell:
		"""
		conda deactivate && conda activate salmon
		featureCounts -T {params.threads} --primary -p -s {params.s} -t gene -g gene_id -a {input.GTF} -o {output.GENE_COUNTS} {input.BAM}
		"""