# RodentAddiction
Annika Vannan, email: avannan@asu.edu \
https://github.com/SexChrLab/RodentAddiction \
Last modified: 06/11/2022

## Quick-Start Guide
1. Clone this repository.
2. Set up Conda environment.
3. Download FASTQ files for each project from SRA. Rename them to match naming conventions as described in this README.
4. Using the appropriate config file (`{dataset}_config.json`), run Snakemake scripts (use SLURM scheduler if available). For each dataset, run scripts in the following order: 
    a. `{dataset}_01_fastq_02_trimmed_fastq.snk`
    b. `{dataset}_03_alignment_04_processing.snk` and `{dataset}_03_alignment_SALMON.snk`
    c. `{dataset}_05_counts.snk`
5. For each dataset, use the sample manifests/phenotype file (`{dataset}_pheno.csv`) and final gene counts run the R script `{dataset}_differential_expression.R`. Note that all R scripts in this repository should be run line-by-line to check results.
6. Download GTEX data and process data according to `gtex_subset.slurm.sh`.
7. Run R scripts in **post_processing/** directory to compare DEGs from each dataset. Start with `individual_overlaps_homologs.R`, and then run the other 3 scripts in any order. These scripts require downloaded data (**post_processing/downloaded_data/**). The script `individual_gtex.R` requires files from the previous step.


## Project Description
The purpose of this project is to identify genes relevant to drug-motivated behavior in rodent models of cocaine addiction and develop a pipeline for narrowing candidate genes from RNA-seq data. RNA-seq datasets were obtained from the papers listed below. Data from each paper/experiment are referred to by their first author (e.g. Carpenter data) throughout this README and other files in the directory.

Carpenter, M.D., Hu, Q., Bond, A.M. *et al.* Nr4a1 suppresses cocaine-induced behavior via epigenetic regulation of homeostatic target genes. *Nat Commun* **11,** 504 (2020). https://doi.org/10.1038/s41467-020-14331-y

Powell, G.L., Vannan, A., Bastle, R.M. *et al.* Environmental enrichment during forced abstinence from cocaine self-administration opposes gene network expression changes associated with the incubation effect. *Sci Rep* **10,** 11291 (2020). https://doi.org/10.1038/s41598-020-67966-8

Walker, D.M., Cates, H.M, Loh, Y.E. *et al.* Cocaine self-administration alters transcriptome-wide responses in the brain's reward circuitry. *Biol Psychiatry* **15**, 84(12):867-880. (2018) https://doi.org/10.1016/j.biopsych.2018.04.009

## GitHub Repository Structure
This repository contains directories for each dataset (**carpenter/**, **powell/**, **walker/**). Within these, you will find two subdirectories: (1) **00_misc_and_scripts/**, which contains the Snakemake scripts and config files used to process the data (**scripts/**), the final gene counts after analysis with HISAT2 and featureCounts (**gene_counts/**), and miscellaneous files including sample manifests (**miscellaneous/**)and the initial sizes of each FASTQ file, and (2) **01_fastq/**, where you should place the FASTQ files downloaded from SRA.

The repository also has a directory scripts for the post-processing of all 3 datasets (**post_processing/**). Within this directory are several scripts that should be used to analyze all 3 datasets together after differential expression analysis, along with other publicly available data necessary for these scripts (**downloaded_data/**) and results from the post-processing analysis (**results/**).

## Data Overview & Download
Data descriptions below only include basic information relevant to the present analysis. Please see original papers for further detail.

**Carpenter** (Paired End, Reverse Stranded) \
Experiment Description: Mice self-administered either cocaine or saline for 10 days, and were then placed into either 1d or 28d abstinence. \
The C1abs and C28abs conditions could not be compared because of extreme batch effects - each condition was sequenced in a separate run. See `carpenter_config.json` for FASTQ headers that verify this. Because of this, the Carpenter snakemake scripts and the R scripts are set to run only the 28d samples for a comparison of S28abs vs. C28abs. \
Species: *Mus musculus* (mouse) \
Strain: C57BL/6J \
Sex: Male \
Tissue: Nucleus accumbens (whole) \
Treatment Groups:
- Cocaine, 1d abstinence (C1abs)
- Cocaine, 28d abstinence (C28abs)
- Saline, 1d abstinence (S1abs)
- Saline, 28d abstinence (S21abs)

**Walker** (Single End, Reverse Stranded) \
Experiment Description: Mice self-administered cocaine or saline for 10-15 days, and then experienced either 1d of withdrawal, or 30d withdrawal and a saline challenge before sacrifice. \
The C1no and C30sal conditions were not compared because the conditions don't just differ with length of abstinence, but also with no/saline injection and nothing/context re-exposure. Because of this, the Walker snakemake files and the R files are set to run only the 30d samples for a comparison of S30sal vs. C30sal. \
Species: *Mus musculus* (mouse) \
Strain: C57BL/6J \
Sex: Male \
Tissue: Nucleus accumbens (whole) \
Treatment Groups:
- Cocaine, 1d withdrawal, no challenge injection (C1no)
- Saline, 1d withdrawal, no challenge injection (S1no)
- Cocaine, 30d withdrawal, saline challenge (C30sal)
- Saline, 30d withdrawal, saline challenge (S30sal)

**Powell** (Single End, Unstranded) \
Experiment Description: Rats self-administered cocaine for >= 15 days, and were then placed into abstinence either in their original isolated housing or in an enriched environment with social interaction and novel toys. After either 1d or 21d abstinence, rats were measured for cocaine-seeking behavior in a 1h cue reactivity test. Enrichment conditions were not of interest for this project, so Powell snakemake and R scripts only run the isolated samples. \
Species: *Rattus norvegicus* (rat) \
Strain: Sprague-Dawley \
Sex: Male \
Tissue: Nucleus accumbens shell\
Treatment Groups:
- Cocaine, 1d abstinence in isolated housing (CI1cue)
- Cocaine, 21d abstinence in isolated housing (CI21cue)

FASTQ files for each experiment can be downloaded from NCBI's Sequence Read Archive (SRA) under the project IDs SRP234876 (Carpenter), SRP246331 (Powell), and SRP132477 (Walker). Not all of the original files are used in this project. See above for specific information on which samples were used.

*NOTE:* FASTQ files from SRA are stripped of their headers, which can contain valuable information. The original FASTQ headers have been obtained from the authors and are included in the manifests and JSON config files.

## Getting Reference Files
Reference genomes (DNA) transcriptomes (cDNA) and annotation files were downloaded from Ensembl's FTP site, version 105. \
This pipeline uses only the transcriptome files.

Mouse transcriptome: GRCm39 \
Rat transcriptome: mRatBN7.2

```bash
# Download mouse genome, transcriptome, and annotation file: GRCm39, database version 105
wget ftp://ftp.ensembl.org/pub/release-105/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna.toplevel.fa.gz
wget ftp://ftp.ensembl.org/pub/release-105/fasta/mus_musculus/cdna/Mus_musculus.GRCm39.cdna.all.fa.gz
wget ftp://ftp.ensembl.org/pub/release-105/gtf/mus_musculus/Mus_musculus.GRCm39.105.gtf.gz
gzip -d Mus_musculus.GRCm39.dna.toplevel.fa.gz
gzip -d Mus_musculus.GRCm39.cdna.all.fa.gz
gzip -d Mus_musculus.GRCm39.105.gtf.gz

# Download rat transcriptome and annotation file: mRatBN7.2, database version 105
wget ftp://ftp.ensembl.org/pub/release-105/fasta/rattus_norvegicus/dna/Rattus_norvegicus.mRatBN7.2.dna.toplevel.fa.gz
wget ftp://ftp.ensembl.org/pub/release-105/fasta/rattus_norvegicus/cdna/Rattus_norvegicus.mRatBN7.2.cdna.all.fa.gz
wget ftp://ftp.ensembl.org/pub/release-105/gtf/rattus_norvegicus/Rattus_norvegicus.mRatBN7.2.105.gtf.gz
gzip -d Rattus_norvegicus.mRatBN7.2.dna.toplevel.fa.gz
gzip -d Rattus_norvegicus.mRatBN7.2.cdna.all.fa.gz
gzip -d Rattus_norvegicus.mRatBN7.2.105.gtf.gz
```

### Create or Download Reference Indexes
This pipeline uses HISAT2 for alignment of transcriptomes. Additionally, it uses Salmon to determine library type.

#### Creating HISAT2 Indexes
The HISAT2 indexes used here were created directly from the transcriptome and annotation files downloaded from ENSEMBL. Successful index building will result in the creation of a folder *{transcriptome_name}_transcriptome_tran* containing 8 index files with the names *{transcriptome_build}_cdna_hisat2*.{1-8}.ht2*

```bash
hisat2-build -f GRCm39_transcriptome.fa GRCm39_transcriptome_tran # Mouse
hisat2-build -f mRatBN7.2_transcriptome.fa mRatBN7.2_transcriptome_tran # Rat
```

#### Download Salmon Indexes
Salmon was used to check library type, using transcriptome indices from refgenie for older versions of the Ensembl transcriptomes - mm10/GRCm38 (mouse) and Rn6 (rat). These were downloaded from the refegenie website (refgenomes.databio.org) for the mm10_cdna and rn6_cdna indexes. \

Digests: \
mm10_cdna = fa159612d40b1bedea9a279eb24999b3d27145f9dd70dcca \
rn6_cdna = d982f0b1888504cccc131d5fa2c11eb3522b9a8e95e0ea89

## Before Using These Scripts
### Conda Environment & Dependencies
All data were processed on a HPC cluster using SLURM and Snakemake. To create a conda environment with the tools used in data processing, install conda (e.g. Miniconda; I used Miniconda3), version 4.9.2. Then run the following to create the environment and install the necessary dependencies:

```bash
conda create --name rodent_addiction # Create environment
conda activate rodent_addiction # Activate environment
conda config --add channels bioconda # Add bioconda channel
conda config --add channels conda-forge # Add conda-forge channel

# Install packages
conda install perl==5.26.2 python==3.6.12 snakemake==5.31.1 fastqc==0.11.9 multiqc==1.9 \
bbmap==38.90 hisat2==2.2.1 samtools==1.7 bamtools==2.5.1 stringtie==2.1.4 subread==2.0.1
```

You will also need to create a second environment with just salmon, bedtools, and mashmap:

```bash
conda create --name salmon
conda activate salmon

conda install salmon==1.1.0 bedtools==2.28.0 mashmap==2.0
```

## Before Using These Scripts
### Troubleshooting Issues with Perl
You may run into problems when using tools that require Perl if the HPC cluster has a different version of Perl than the custom environment. Perl commands may attempt to access the HPC version of Perl instead of the local installation. To reset the Perl library, use the following code with appropriate modifications based on the location of your Perl installation  (only needs to be run once):

```bash
unset PERL5LIB
export PERL5LIB=$PERL5LIB:/home/avannan/miniconda3/envs/rodent_addiction/lib/site_perl/5.26.2
```

### Using SLURM on the HPC
All sbatch scripts utilizing tools in the rodent_addiction environment should start with:

```bash
# Activate environment
source ~/miniconda3/etc/profile.d/conda.sh
conda activate rodent_addiction
```

### Other Important Considerations
Before using these Snakemake or R scripts, see the beginning of the script for directory variables that should be changed by the user.

## Workflow and Script Order
The structure of the GitHub directory is described at the end of this document.

Snakemake scripts are used to process FASTQ files, align to a reference genome, and get raw gene counts. For each dataset, run Snakemake scripts in the following order:
1. `{dataset}_01_fastq_02_trimmed_fastq.snk` - Creates files in the 01_fastq and 02_trimmed_fastq directories
2. `{dataset}_03_alignment_04_processing.snk` - Creates files in the 03_alignment and 04_processing directories
3. `{dataset}_05_counts.snk` - Creates files in the 05_counts directory

R scripts are used to perform differential expression analysis for each dataset and then analyze the 3 datasets together. R scripts should be run line-by-line in the following order.
1. `final_{dataset}_DE.R` (Must be run before any other R scripts)
2. `individual_homologs.R` (Obtains Candidate gene lists for pipeline)
3. `candidates_overlaps_homologs.R` (OPTIONAL)
4. `individual_conservation.R`
5. `individual_gtex.R`
6. `power_analysis.R` (OPTIONAL; Can be run any time after the `final_{dataset}_DE.R` scripts)

Note that `individual_conservation.R` requires data downloads from Cardoso-Moreira et al. 2020 (PubMed ID: 33113372) and `individual_gtex.R` requires data from GTEx (https://gtexportal.org/home/), which are not included here - see "Directory Structure & Obtaining Additional Files".

### Naming Conventions for Initial FASTQ files
Raw FASTQ files should be named according to the following conventions: \
Samples with 1 run: **{dataset}\_{SRA_Sample_ID}\_{Treatment_Group}.fq** \
Samples with multiple runs: **{dataset}\_{SRA_Sample_ID}\_{Treatment_Group}_run{#}.fq**

Carpenter FASTQ files should be located in the **carpenter/01_fastq** directory. Powell and Walker have multiple FASTQ files per sample, and should be located in the **powell/01_fastq/by_run** and **walker/01_fastq/by_run** directories, respectively.

Note: One sample (SRS2926581_S1no) from the Walker experiment was produced in a pilot experiment and contains only one run. Its FASTQ file should be renamed `walker_SRS2926581_S1no_cat.fq` and placed in the **walker/01_fastq/concatenated** directory. \
If you download the Powell files from SRA, they will *already be concatenated* so this should be taken into account when running these scripts.

## Note on BBDuk Parameters
For all datasets, FASTQ files were trimmed with the following parameters: `qtrim=rl trimq=30 minlen={half of read length} maq=20`. Based on the FastQC output, there were some additional trimming parameters for Powell and Walker. For Powell, poly(A/T) and poly(G) were also trimmed. For Walker, additional adapters were trimmed by adding them to the adapters file included with BBDuk - this modified file is called `adapters_overrepresented.fa` and is located in the miscellaneous files for Walker. 

### Directory Structure & Obtaining Additional Files
The RNA-seq workflow will create directories for each dataset that follow the structures outlined in the 00_misc folder for each dataset. Snakemake will generate most of the directories for you, but to start you will need the following structure and files. You will also need to download FASTQ files from SRA and rename as described above. For each dataset, the final gene counts are located in this repository for convenience (**dataset/gene_counts/**), though they will also be created separately when the scripts are run. * indicates data files that need to be downloaded from external sources by the user. ^ indicates data files that are available upon request, but can also be obtained using custom scripts and publicly available GTEx data.

```bash
├── carpenter
|   ├── 00_misc_and_scripts
│   |   ├── gene_counts
|   |   |   └── carpenter_{sample}_hisat2_featurecounts_genes_s2.txt (12 files)
│   |   ├── miscellaneous
│   |   |   ├── carpenter_file_sizes.txt
│   |   |   └── carpenter_pheno.csv
|   |   └── scripts
│   |       ├── carpenter_01_fastq_02_trimmed_fastq.snk
│   |       ├── carpenter_03_alignment_04_processing.snk
|   |       ├── carpenter_03_alignment_SALMON.snk
│   |       ├── carpenter_05_counts.snk
│   |       ├── carpenter_config.json
|   |       └── final_carpenter_DE.R
|   └── 01_fastq
|       ├── carpenter_{sample}_fq1.fastq (12 files) *
|       └── carpenter_{sample}_fq2.fastq (12 files) *
├── powell
|   ├── 00_misc_and_scripts
│   |   ├── gene_counts
|   |   |   └── powell_{sample}_hisat2_featurecounts_genes.txt (6 files)
│   |   ├── miscellaneous
│   |   |   ├── powell_file_sizes.txt
│   |   |   └── powell_pheno.csv
|   |   └── scripts
│   |       ├── powell_01_fastq_02_trimmed_fastq.snk
│   |       ├── powell_03_alignment_04_processing.snk
|   |       ├── powell_03_alignment_SALMON.snk
│   |       ├── powell_05_counts.snk
│   |       ├── powell_config.json
|   |       └── final_powell_DE.R
|   └── 01_fastq
|       ├── by_run
|       |   └── powell_{sample}_run{1-4}.fastq (24 files; 6 samples) *
|       └── concatenated
├── post_processing
|   ├── results
|   |   ├── gene_info
|   |   |   ├── carpenter_DE_FINAL.txt
|   |   |   ├── powell_DE_FINAL_treat_lane.txt
|   |   |   └── walker_DE_FINAL_treat_batch.txt
|   |   └── other
|   |       └── all_genes_list.RDS
|   └── scripts
|       ├── candidates_overlaps_homologs.R
|       ├── gtex_subset_slurm.sh
|       ├── individual_conservation.R
|       ├── individual_gtex.R
|       ├── individual_homologs.R
|       └── power_analysis.R
└── walker
    ├── 00_misc_and_scripts
    |   ├── gene_counts
    |   |   └── walker_{sample}_hisat2_featurecounts_genes_s2.txt (11 files)
    |   ├── miscellaneous
    |   |   ├── adapters_overrepresented.fa
    |   |   ├── walker_file_sizes.txt
    |   |   └── walker_pheno.csv
    |   └── scripts
    |       ├── walker_01_fastq_02_trimmed_fastq.snk
    |       ├── walker_03_alignment_04_processing.snk
    |       ├── walker_03_alignment_SALMON.snk
    |       ├── walker_05_counts.snk
    |       ├── walker_config.json
    |       └── final_walker_DE.R
    └── 01_fastq
        ├── by_run
        |   └── walker_{sample}_run{1-8}.fastq (80 files; 10 samples) *
        └── concatenated
            └── walker_SRS2926581_S1no_cat.fq *
```

