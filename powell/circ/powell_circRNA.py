
# Need to add:
# miRanda, and cropping the rat files


#####################
## Start of Script ##
#####################

import os

## - USERS SHOULD CHANGE THE FOLLOWING VARIABLES - ##
# Tool Paths
CIRI2_path = "/home/avannan/circRNA_tools/CIRI2/CIRI_v2.0.6/CIRI2.pl"
samtools_path = "/home/avannan/miniconda3/envs/circ_analysis/lib/samtools" # Version 2.5.1
bamtools_path = "/home/avannan/miniconda3/envs/circ_analysis/lib/bamtools" # Version 1.7
# References
bwa_index =  "/scratch/avannan/refs/rat/mRatBN7.2/Rattus_norvegicus.mRatBN7.2.dna.toplevel.fa"
genome_fasta = "/scratch/avannan/refs/rat/mRatBN7.2/Rattus_norvegicus.mRatBN7.2.dna.toplevel.fa"
annotation_gtf = "/scratch/avannan/refs/rat/mRatBN7.2/Rattus_norvegicus.mRatBN7.2.105.gtf "
gene_pred = "/data/CEM/wilsonlab/projects/rodent_addiction/powell/mRatBN7.2_annotation.txt"
# Starting Directory
start_dir = "/data/CEM/wilsonlab/projects/rodent_addiction/"
## -- END -- ##

# Config file
configfile: "powell_config.json"

# Directory Variables
# FASTQs
trimmed_fastq_dir = os.path.join(start_dir, config["trimmed_fastq_dir"]) # Trimmed FASTQs, concatenated by sample
# SAMS, BAMS, and stats
circ_sam_dir = os.path.join(start_dir, config["circ_sam_dir"]) # BWA-aligned SAMs for CIRI2 analysis
circ_bam_dir = os.path.join(start_dir, config["circ_bam_dir"]) # BAMs for getting stats
circ_sort_bam_dir = os.path.join(start_dir, config["circ_sort_bam_dir"]) # BAMs for getting stats
circ_bam_stats_dir = os.path.join(start_dir, config["circ_bam_stats_dir"]) # Stats from BAMs
# Circ identification
ciri2_out_dir = os.path.join(start_dir, config["ciri2_out_dir"]) # CIRI2 identification & counting of circRNAs
ciri2_log_dir = os.path.join(start_dir, config["ciri2_out_dir"], "logs") # Log file directory for CIRI2
circexplorer2_out_dir = os.path.join(start_dir, config["circexplorer2_out_dir"]) # CIRCexplorer2 identification & counting of circRNAs

######################
## All Output Files ##
######################

rule all:
	input:
		# BWA-aligned SAMs
		expand(os.path.join(circ_sam_dir, "powell_{sample}_bwa_trim.sam"), sample = config["samples_ALL"]),
		# BWA-aligned BAMS, including sorted and stats
		expand(os.path.join(circ_bam_dir, "powell_{sample}_bwa_trim.bam"), sample = config["samples_ALL"]),
		expand(os.path.join(circ_sort_bam_dir, "powell_{sample}_bwa_trim_sort.bam"), sample = config["samples_ALL"]),
		expand(os.path.join(circ_bam_stats_dir, "powell_{sample}_bwa_trim_sort_stats.txt"), sample = config["samples_ALL"]),
		# CIRI2 outputs
		expand(os.path.join(ciri2_out_dir, "powell_{sample}_bwa.ciri2.out"), sample = config["samples_ALL"]),
		# CIRCexplorer2 outputs
		expand(os.path.join(circexplorer2_out_dir, "parse/powell_{sample}_back_spliced_junction.bed"), sample = config["samples_ALL"]),
		expand(os.path.join(circexplorer2_out_dir, "annotate/powell_{sample}_circularRNA_known.txt"), sample = config["samples_ALL"]),
		expand(os.path.join(circexplorer2_out_dir, "annotate/logs/CIRCexplorer2_annotate_powell_{sample}.log"), sample = config["samples_ALL"])
		

##############################################################
## Identify circRNAs and analyze for differental expression ##
##############################################################

#############################
## BWA Alignment and Stats ##
#############################

rule bwa_align:
	input:
		TRIMMED_FQ = os.path.join(trimmed_fastq_dir, "powell_{sample}_trim.fastq")
	output:
		TRIMMED_SAM = os.path.join(circ_sam_dir, "powell_{sample}_bwa_trim.sam")
	params:
		threads = 16,
		threshold = 19, # Recommended by CIRI2
		bwa_files = bwa_index
	shell:
		"""
		bwa mem -t {params.threads} -T {params.threshold} {params.bwa_files} {input.TRIMMED_FQ} > {output.TRIMMED_SAM}
		"""

rule bwa_sam_to_bam:
	input:
		BWA_SAM = os.path.join(circ_sam_dir, "powell_{sample}_bwa_trim.sam")
	output:
		BWA_BAM = os.path.join(circ_bam_dir, "powell_{sample}_bwa_trim.bam")
	params:
		samtools = samtools_path
	shell:
		"""
		samtools view -b {input.BWA_SAM} > {output.BWA_BAM}
		"""

rule bwa_sort:
	input:
		BWA_UNSORTED_BAM = os.path.join(circ_bam_dir, "powell_{sample}_bwa_trim.bam")
	output:
		BWA_SORTED_BAM = os.path.join(circ_sort_bam_dir, "powell_{sample}_bwa_trim_sort.bam")
	params:
		bamtools = bamtools_path
	shell:
		"""
		bamtools sort -in {input.BWA_UNSORTED_BAM} -out {output.BWA_SORTED_BAM}
		"""

rule bwa_sorted_stats:
	input:
		BWA_SORTED_BAM = os.path.join(circ_sort_bam_dir, "powell_{sample}_bwa_trim_sort.bam")
	output:
		BWA_SORTED_STATS = os.path.join(circ_bam_stats_dir, "powell_{sample}_bwa_trim_sort_stats.txt")
	params:
		bamtools = bamtools_path
	shell:
		"""
		bamtools stats -in {input.BWA_SORTED_BAM} > {output.BWA_SORTED_STATS}
		"""

##################################################
## CircRNA Identification & Counting with CIRI2 ##
##################################################

rule ciri2_bwa:
	input:
		SAM = os.path.join(circ_sam_dir, "powell_{sample}_bwa_trim.sam")
	output:
		CIRI2_OUT = os.path.join(ciri2_out_dir, "powell_{sample}_bwa.ciri2.out")
	params:
		ciri2 = CIRI2_path,
		genome = genome_fasta,
		log_dir = ciri2_log_dir,
		gtf = annotation_gtf
	shell:
		"""
		perl {params.ciri2} -I {input.SAM} -O {output.CIRI2_OUT} -F {params.genome} -A {params.gtf} -0
		"""

##########################################################
## CircRNA Identification & Counting with CIRCexplorer2 ##
##########################################################

rule circexplorer2_bwa_parse:
	input:
		SAM = os.path.join(circ_sam_dir, "powell_{sample}_bwa_trim.sam")
	output:
		CIRCEXPLORER2_BED = os.path.join(circexplorer2_out_dir, "parse/powell_{sample}_back_spliced_junction.bed"),
		PARSE_LOG = os.path.join(circexplorer2_out_dir, "parse/logs/CIRCexplorer2_parse_powell_{sample}.log")
	shell:
		"""
		source ~/miniconda3/etc/profile.d/conda.sh && conda deactivate && conda activate circexplorer2 &&
		CIRCexplorer2 parse -t BWA {input.SAM} -b {output.CIRCEXPLORER2_BED} > {output.PARSE_LOG}
		"""

rule circexplorer2_annotate:
	input:
		CIRCEXPLORER2_BED = os.path.join(circexplorer2_out_dir, "parse/powell_{sample}_back_spliced_junction.bed")
	output:
		KNOWN_CIRCS = os.path.join(circexplorer2_out_dir, "annotate/powell_{sample}_circularRNA_known.txt"),
		ANNOTATE_LOG = os.path.join(circexplorer2_out_dir, "annotate/logs/CIRCexplorer2_annotate_powell_{sample}.log")
	params:
		gene_pred = gene_pred,
		genome = genome_fasta
	shell:
		"""
		source ~/miniconda3/etc/profile.d/conda.sh && conda deactivate && conda activate circexplorer2 &&
		CIRCexplorer2 annotate -r {params.gene_pred} -g {params.genome} -b {input.CIRCEXPLORER2_BED} \
		-o {output.KNOWN_CIRCS} > {output.ANNOTATE_LOG}
		"""
