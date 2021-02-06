#!/bin/bash
#SBATCH --job-name=practice          # Job name
#SBATCH --mail-type=END,FAIL         # Send a notification when the job stops or fails
#SBATCH --mail-user=avannan@asu.edu  # Send-to address
#SBATCH --qos=normal                 # Quality of service (priority)
#SBATCH -N 1                         # Number of nodes
#SBATCH -n 2                         # Number of tasks per node
#SBATCH -c 4                         # Number of CPUs per task
#SBATCH -t 00:01:00                  # Wall time (DD:HH:MM)

# Activate environment
conda activate rodent_addiction

# Go to the directory where the snakefile is.
cd /scratch/avannan/rodent_addiction/carpenter/

# Run the snakefile.
snakemake --snakefile <snakemakefilename> -j <numberofparalleljobs> --cluster "sbatch -n <numberofcores> --mem=<memorypernode> -c 4 -t 24:00:00" --rerun-incomplete
# <numberofparalleljobs> == 43, <numberofcores> == 4, <memorypernode> == 12000, <t> == 24:00:00


From FastQC:
AGATCGGAAGAG

TruSeq2-PE.fa
>PrefixPE/1
AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT
>PrefixPE/2
CAAGCAGAAGACGGCATACGAGATCGGTCTCGGCATTCCTGCTGAACCGCTCTTCCGATCT
>PCR_Primer1
AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT
>PCR_Primer1_rc
AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT
>PCR_Primer2
CAAGCAGAAGACGGCATACGAGATCGGTCTCGGCATTCCTGCTGAACCGCTCTTCCGATCT
>PCR_Primer2_rc
AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGATCTCGTATGCCGTCTTCTGCTTG
>FlowCell1
TTTTTTTTTTAATGATACGGCGACCACCGAGATCTACAC
>FlowCell2
TTTTTTTTTTCAAGCAGAAGACGGCATACGA

TruSeq3-PE.fa
>PrefixPE/1
TACACTCTTTCCCTACACGACGCTCTTCCGATCT
>PrefixPE/2
GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT

TruSeq3-PE-2.fa
>PrefixPE/1
TACACTCTTTCCCTACACGACGCTCTTCCGATCT
>PrefixPE/2
GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT
>PE1
TACACTCTTTCCCTACACGACGCTCTTCCGATCT
>PE1_rc
AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTA
>PE2
GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT
>PE2_rc
AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC


# POWELL 2020
#RNA was prepared using the Nugen Ovation RNA-seq system 
#(#7102-32; Tecan Genomics, Inc., Redwood City, California USA). 
#cDNA was quantified using Nanodrop and sheared to approximately 300 base pair fragments. 
#Libraries were generated using the Kapa Biosystem library preparation kit. 
#Fragments were then end repaired and A-tailed. 
#Sequencing was performed on a 1 × 75 flow cell using an Illumina NextSeq500 instrument.

# Save for later
				# HISAT2 BAMs with duplicates removed
				expand(hisat2_rmdup_dir + "carpenter_{sample}_hisat2_rmdup.bam", sample = config["samples_ALL"]),
				# Indexes for HISAT2 BAMs with duplicates removed
				expand(hisat2_rmdup_dir + "carpenter_{sample}_hisat2_rmdup.bam.bai", sample = config["samples_ALL"]),
				# Metrics for HISAT2 BAMs with duplicates removed
				expand(hisat2_rmdup_metrics_dir + "carpenter_{sample}_hisat2_rmdup_marked_dup_metrics.txt", sample = config["samples_ALL"])

rule hisat2_picard_remove_duplicates:
		input:
				HISAT2_SORTED_BAM = hisat2_sorted_bams_dir + "carpenter_{sample}_hisat2_pair_trim_sort.bam"
		output:
				HISAT2_SORTED_RMDUP_BAM = hisat2_rmdup_dir + "carpenter_{sample}_hisat2_rmdup.bam",
				HISAT2_BAM_METRICS = hisat2_rmdup_metrics_dir + "carpenter_{sample}_hisat2_rmdup_marked_dup_metrics.txt"
		message: "Removing duplicates from {input.HISAT2_SORTED_BAM}."
		params:
				picard = picard_path
		shell:
				"""
				{params.picard} -Xmx14g MarkDuplicates -I {input.HISAT2_SORTED_BAM} -O {output.HISAT2_SORTED_RMDUP_BAM} -M {output.HISAT2_BAM_METRICS}
				"""

rule hisat2_index_rmdup_bams:
		input:
				HISAT2_RMDUP_BAM = hisat2_rmdup_dir + "carpenter_{sample}_hisat2_rmdup.bam"
		output: hisat2_rmdup_dir + "carpenter_{sample}_hisat2_rmdup.bam.bai"
		params:
				bamtools = bamtools_path
		shell:
				"""
				{params.bamtools} index -in {input.HISAT2_RMDUP_BAM}
				"""

rule hisat2_picard_readgroups:
		input:
				HISAT2_SORTED_RMDUP_BAM = hisat2_rmdup_dir + "carpenter_{sample}_hisat2_rmdup.bam"
		output:
				HISAT2_RMDUP_RG_BAM = hisat2_rmdup_rg_dir + "carpenter_{sample}_hisat2_rmdup_rg.bam"
		params:
				picard = picard_path,
				id = lambda wildcards: config[wildcards.sample]["ID"], # ID
				sm = lambda wildcards: config[wildcards.sample]["SM"], # Sample name; Required
				lb = lambda wildcards: config[wildcards.sample]["LB"], # Library; Required
				pu = lambda wildcards: config[wildcards.sample]["PU"], # Platform unit; Required
				pl = lambda wildcards: config[wildcards.sample]["PL"], # Platform; Required
		message: "Adding readgroups to {input.HISAT2_SORTED_RMDUP_BAM}."
		shell:
				"""
				{params.picard} -Xmx14g AddOrReplaceReadGroups I={input.HISAT2_SORTED_RMDUP_BAM} O={output.HISAT2_RMDUP_RG_BAM}
				RGID={params.id} RGPU={params.pu} RGSM={params.sm} RGPL={params.pl} RGLB={params.lb} VALIDATION_STRINGENCY=LENIENT
				"""

rule hisat2_sorted_bam_list:
		input:
		output:
				HISAT2_SORTED_BAM_LIST = hisat2_align_dir + "carpenter_hisat2_sorted_bam_list.txt"
		message: "Creating list of HISAT2 BAMs as a .txt file for use in StringTie merge."
		params:
				hisat2_align_dir = hisat2_align_dir
		shell:
				"""
				dir {param.hisat2_processed_dir} > {output.HISAT2_SORTED_BAM_LIST}
				"""


# CARPENTER INFO
#Libraries were prepared using the TruSeq Stranded mRNA HT Sample Prep Kit protocol
#(Illumina, San Diego, CA). Briefly, poly A selection and fragmentation of 300 ng of RNA was
#converted to cDNA with random hexamers. Adapters were ligated and samples were sizeselected with AMPur XP beads (Beckman Coulter, Brea, CA). 
#Barcode bases (6 bp) were introduced at one end of the adaptors during PCR amplification steps. Library size and 
#concentration was assessed using Tape Station (Life Technologies, Grand Island, NY) before
#sequencing. Libraries were pooled for multiplexing (4 pools of ~60 samples with each group and
#brain region equally represented across each pool) and sequenced on a HighSeq2500 System
#using V4 chemistry with 50 base pair single-end reads at GeneWiz LLC (South Plainfield, NJ).
#Each pool was sequenced 8 times with the goal of obtaining ~25 million reads per sample. Initial
#quality control assessments revealed 43 samples, which did not meet standards for read depth
#and were excluded from analysis. Therefore, the final number of samples included in the analysis
#were between 5 – 8 per group apart from the CS group in VTA (N = 3). 

# WALKER INFO
#RNA integrity number (RIN) and concentration were assessed using a Bioanalyzer (Agilent). 
#Libraries were constructed using the ScriptSeq Complete Gold Kit (Epicentre, Illumina) preceded by ribosomal RNA depletion.