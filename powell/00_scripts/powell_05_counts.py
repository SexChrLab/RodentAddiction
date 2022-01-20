# TO-DO
# Create rule for powell_prepDE.txt if still using StringTie

########################################################
## Workflow for counting reads after HISAT2 alignment ##
########################################################

## - Information about Project & Samples - ##
# Original fastq files come from our laboratory and are separated into 4 runs per sample.
# All sample names (indicated by {sample}) include both the SRA sample ID (SRS#) and the treatment group - e.g. SRS6085261_CI1cue
# Run # is also listed for applicable files.
# Concatenated files can be downloaded from SRA (though FASTQ headers have been stripped).
# SRA Project ID: SRP246331; Bioproject ID: PRJNA604189; Gene Expression Omnibus Series ID: GSE144606

## - Using This Script - ##
# See README file in GitHub before running these scripts.
# Remember to change variables under the heading "USERS SHOULD CHANGE THE FOLLOWING VARIABLES"


#####################
## Start of Script ##
#####################

import os

## - USERS SHOULD CHANGE THE FOLLOWING VARIABLES - ##
# Config file
configfile: "powell_config.json"

# Tool Paths
stringtie_path = "/home/avannan/miniconda3/envs/rodent_addiction/bin/stringtie" # Version 2.1.4
featurecounts_path = "/home/avannan/miniconda3/envs/rodent_addiction/bin/featureCounts/" # Subread, version 2.0.1
# References & Reference Directories
annotation_gtf = "/scratch/avannan/refs/rat/mRatBN7.2/Rattus_norvegicus.mRatBN7.2.105.gtf"
# Starting Directory
start_dir = "/data/CEM/wilsonlab/projects/rodent_addiction/"
# Script Paths
prepDE_path = "/data/CEM/wilsonlab/projects/rodent_addiction/prepDE.py3" # Last updated for StringTie version 2.0.4 release

## -- END SECTION FOR CHANGING VARIABLES -- ##


# Directory Variables
# Miscellaneous
misc_dir = os.path.join(start_dir, config["misc_dir"]) # Contains miscellaneous files, including text file for converting to raw counts
# FASTQs
cat_fastq_dir = os.path.join(start_dir, config["cat_fastq_dir"]) # Initial FASTQs, concatenated by sample
trimmed_fastq_dir = os.path.join(start_dir, config["trimmed_fastq_dir"]) # Trimmed FASTQs, concatenated by sample
# Processing
hisat2_sorted_bams_dir = os.path.join(start_dir, config["hisat2_sorted_bams_dir"]) # Sorted HISAT2 BAM files
# Counting with StringTie
hisat2_stringtie_dir = os.path.join(start_dir, config["hisat2_stringtie_dir"]) # Main HISAT2/StringTie directory
hisat2_stringtie_count_tran_dir = os.path.join(start_dir, config["hisat2_stringtie_count_tran_dir"]) # Transcript counts
hisat2_stringtie_cov_tran_dir = os.path.join(start_dir, config["hisat2_stringtie_cov_tran_dir"]) # Covered transcripts
hisat2_stringtie_gene_dir = os.path.join(start_dir, config["hisat2_stringtie_gene_dir"]) # Genes
hisat2_stringtie_final_count_dir = os.path.join(start_dir, config["hisat2_stringtie_final_count_dir"]) # Final raw counts from prepDE.py3
# Counting with featureCounts
hisat2_featurecounts_dir = os.path.join(start_dir, config["hisat2_featurecounts_dir"]) # Main HISAT2/featureCounts directory
hisat2_featurecounts_gene_dir = os.path.join(start_dir, config["hisat2_featurecounts_gene_dir"]) # Genes


######################
## All Output Files ##
######################

rule all:
	input:
		# COUNTING WITH STRINGTIE AFTER HISAT2 ALIGNMENT
		# Tab files are also created in the same directory as the counted transcripts, but not listed here
		expand(os.path.join(hisat2_stringtie_count_tran_dir, "{sample}/powell_{sample}_hisat2_stringtie_count_transcripts.gtf"), sample = config["Iso_ALL"]),
		expand(os.path.join(hisat2_stringtie_cov_tran_dir, "powell_{sample}_hisat2_stringtie_covered_transcripts.gtf"), sample = config["Iso_ALL"]),
		expand(os.path.join(hisat2_stringtie_gene_dir, "powell_{sample}_hisat2_stringtie_genes.txt"), sample = config["Iso_ALL"]),
		os.path.join(hisat2_stringtie_final_count_dir, "powell_gene_count_matrix.csv"),
		os.path.join(hisat2_stringtie_final_count_dir, "powell_transcript_count_matrix.csv"),

		# COUNTING WITH FEATURECOUNTS AFTER HISAT2 ALIGNMENT
		expand(os.path.join(hisat2_featurecounts_gene_dir, "powell_{sample}_hisat2_featurecounts_genes.txt"), sample = config["Iso_ALL"])

###############
## 05_counts ##
###############

######################################################
## Counting with StringTie (after HISAT2 alignment) ##
######################################################

rule hisat2_stringtie_count:
	priority: 1
	input:
		BAM = os.path.join(hisat2_sorted_bams_dir, "powell_{sample}_hisat2_trim_sort.bam"),
		GUIDE_GTF = annotation_gtf
	output:
		TRANSCRIPT_GTF = os.path.join(hisat2_stringtie_count_tran_dir, "{sample}/powell_{sample}_hisat2_stringtie_count_transcripts.gtf"),
		COVERED_TRANSCRIPTS = os.path.join(hisat2_stringtie_cov_tran_dir, "powell_{sample}_hisat2_stringtie_covered_transcripts.gtf"),
		GENE_ABUNDANCES = os.path.join(hisat2_stringtie_gene_dir, "powell_{sample}_hisat2_stringtie_genes.txt")
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

rule run_prepDE:
	priority: 10
	input: 
		PREP_TXT = os.path.join(misc_dir, "powell_prepDE_iso.txt") # Contains sample IDs and assembled transcript GTF locations
	output:
		GENE_MATRIX = os.path.join(hisat2_stringtie_final_count_dir, "powell_gene_count_matrix.csv"),
		TRANSCRIPT_MATRIX = os.path.join(hisat2_stringtie_final_count_dir, "powell_transcript_count_matrix.csv")
	params:
		gtf_dir = hisat2_stringtie_count_tran_dir,
		prepDE3 = prepDE_path,
		length = 75 # Average read length
	shell:
		"""
		mkdir -p {params.gtf_dir}
		python3 {params.prepDE3} -i {params.gtf_dir} -l {params.length} -g {output.GENE_MATRIX} -t {output.TRANSCRIPT_MATRIX} -v
		"""

##########################################################
## Counting with featureCounts (after HISAT2 alignment) ##
##########################################################

rule hisat2_featurecounts_gene_count_per_sample:
	priority: 1
	input:
		BAM = os.path.join(hisat2_sorted_bams_dir, "powell_{sample}_hisat2_trim_sort.bam"),
		GTF = annotation_gtf
	output:
		GENE_COUNTS = os.path.join(hisat2_featurecounts_gene_dir, "powell_{sample}_hisat2_featurecounts_genes.txt")
	message: "Counting {wildcards.sample} reads (single end, unstranded) using featureCounts for primary alignments only."
	params:
		s = 0, # Strandness, 0 for unstranded
		threads = 5
	shell:
		"""
		featureCounts -T {params.threads} --primary -s 0 -t gene -g gene_id -a {input.GTF} -o {output.GENE_COUNTS} {input.BAM}
		"""