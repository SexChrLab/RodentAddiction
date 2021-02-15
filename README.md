# RodentAddiction

Annika Vannan, email: avannan@asu.edu \
https://github.com/SexChrLab/RodentAddiction \
Last modified: 02/15/2021

## Project Description
The purpose of this project is to identify genes relevant to drug-motivated behavior in rodent models of cocaine addiction. RNA-seq datasets were obtained from the papers listed below. Data from each paper/experiment are referred to by their first author (e.g. Carpenter data) throughout this README and other files in the directory.

Carpenter, M.D., Hu, Q., Bond, A.M. *et al.* Nr4a1 suppresses cocaine-induced behavior via epigenetic regulation of homeostatic target genes. *Nat Commun* **11,** 504 (2020). https://doi.org/10.1038/s41467-020-14331-y

Powell, G.L., Vannan, A., Bastle, R.M. *et al.* Environmental enrichment during forced abstinence from cocaine self-administration opposes gene network expression changes associated with the incubation effect. *Sci Rep* **10,** 11291 (2020). https://doi.org/10.1038/s41598-020-67966-8

Walker, D.M., Cates, H.M, Loh, Y.E. *et al.* Cocaine self-administration alters transcriptome-wide responses in the brain's reward circuitry. *Biol Psychiatry* **15**, 84(12):867-880. (2018) https://doi.org/10.1016/j.biopsych.2018.04.009

## GitHub Directory Structure
The folder **all/** contains scripts and other files pertaining to the overall project. There are also separate folders for each the Carpenter, Powell, and Walker experiments that contain experiment-specific files and scripts. The folder **junk/** contains junk files and scripts that I'm not ready to part with yet.

```bash
├── all
│   ├── all_manifest.xlsx
│   └── RNASeqWorkflow02062021.pdf
├── carpenter
│   ├── carpenter_01_fastq_02_trimmed_fastq.py
│   ├── carpenter_03_alignment_04_processing.py
│   ├── carpenter_05_counts.py
│   ├── carpenter_config.json
|   └── carpenter_directory_structure.txt
├── junk
│   └── (various files)
├── powell
│   ├── powell_01_fastq_02_trimmed_fastq.py
│   ├── powell_03_alignment_04_processing.py
│   ├── powell_05_counts.py
│   ├── powell_config.json
|   └── powell_directory_structure.txt
├── walker
│   ├── walker_01_fastq_02_trimmed_fastq.py
│   ├── walker_03_alignment_04_processing.py
│   ├── walker_05_counts.py
│   ├── walker_config.json
|   ├── walker_directory_structure.txt
|   └── walker_file_sizes.txt
├── .gitattributes
├── LICENSE
└── README.md
```

## Data Overview & Download
Data descriptions below only include basic information relevant to the present analysis. Please see original papers for further detail.

**Carpenter** (Paired End) \
Experiment Description: Mice self-administered either cocaine or saline for 10 days, and were then placed into either 1d or 28d abstinence. \
Species: *Mus musculus* (mouse) \
Strain: C57BL/6J \
Sex: Male \
Tissue: Nucleus accumbens (whole) \
Treatment Groups:
- Cocaine, 1d abstinence (C1abs)
- Cocaine, 28d abstinence (C28abs)
- Saline, 1d abstinence (S1abs)
- Saline, 28d abstinence (S21abs)

**Powell** (Single End) \
Experiment Description: Rats self-administered cocaine for >= 15 days, and were then placed into abstinence either in their original isolated housing or in an enriched environment with social interaction and novel toys. After either 1d or 21d abstinence, rats were measured for cocaine-seeking behavior in a 1h cue reactivity test. \
Species: *Rattus norvegicus* (rat) \
Strain: Sprague-Dawley \
Sex: Male \
Tissue: Nucleus accumbens shell\
Treatment Groups:
- Cocaine, 1d abstinence in isolated housing (CI1cue)
- Cocaine, 21d abstinence in isolated housing (CI21cue)
- Cocaine, 1d abstinence in enriched housing (CE1cue)
- Cocaine, 21d abstinence in enriched housing (CE21cue)

**Walker** (Single End) \
Experiment Description: Mice self-administered cocaine or saline for 10-15 days, and then experienced either 1d of withdrawal, or 30d withdrawal and a saline challenge before sacrifice. \
Species: *Mus musculus* (mouse) \
Strain: C57BL/6J \
Sex: Male \
Tissue: Nucleus accumbens (whole) \
Treatment Groups:
- Cocaine, 1d withdrawal, no challenge injection (C1no)
- Saline, 1d withdrawal, no challenge injection (S1no)
- Cocaine, 30d withdrawal, saline challenge (C30sal)
- Saline, 30d withdrawal, saline challenge (S30sal)

FASTQ files for each experiment can be downloaded from NCBI's Sequence Read Archive (SRA) under the project IDs SRP234876 (Carpenter), SRP246331 (Powell), and SRP132477 (Walker). Not all of the original files are used in this project. See the JSON or manifest files for specific information about the samples that were used.

*NOTE:* FASTQ files from SRA are stripped of their headers, which can contain valuable information. The original FASTQ headers have been obtained from the authors and are included in the experiment manifests.

**Reference Files** \
Reference genomes and annotation files were downloaded from ENSEMBL's FTP site.

Mouse genome: GRCm38 (also called mm10) \
Rat genome: Rnor_6.0 (also called Rn6)

```bash
# Download mouse genome and annotation file: GRCm38.p6/mm10, database version 102.38
wget ftp://ftp.ensembl.org/pub/release-102/fasta/rattus_norvegicus/dna/Rattus_norvegicus.Rnor_6.0.dna.toplevel.fa.gz
gzip -d Rattus_norvegicus.Rnor_6.0.dna.toplevel.fa.gz
mv Rattus_norvegicus.Rnor_6.0.dna.toplevel.fa Rn6_genome.fa

wget ftp://ftp.ensembl.org/pub/release-102/gtf/rattus_norvegicus/Rattus_norvegicus.Rnor_6.0.102.gtf.gz
gzip -d Rattus_norvegicus.Rnor_6.0.102.gtf.gz
mv Rattus_norvegicus.Rnor_6.0.102.gtf Rn6_annotation.gtf

# Download rat genome and annotation file: Rnor_6.0, database version 102.6
wget ftp://ftp.ensembl.org/pub/release-102/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.toplevel.fa.gz
gzip -d Mus_musculus.GRCm38.dna.toplevel.fa.gz
mv Mus_musculus.GRCm38.dna.toplevel.fa mm10_genome.fa

wget ftp://ftp.ensembl.org/pub/release-102/gtf/mus_musculus/Mus_musculus.GRCm38.102.gtf.gz
gzip -d Mus_musculus.GRCm38.102.gtf.gz
mv Mus_musculus.GRCm38.102.gtf mm10_annotation.gtf
```

## Using These Scripts
### Conda Environment & Dependencies
All data were processed on a HPC cluster using SLURM and Snakemake. To create a conda environment with the tools used in data processing, install conda (e.g. Miniconda; I used Miniconda3) and run:

```bash
conda create -y rodent_addiction
```

Environment Name: rodent_addiction \
Conda Version: 4.9.2

Activate the environment using:
```bash
conda activate rodent_addiction
```
Note that in some cases it is preferable to use ***source activate*** rather than ***conda activate***.

#### Channels:
```bash
conda config --add channels bioconda
conda config --add channels conda-forge
```
#### Dependencies:
```bash
conda install perl==5.26.2
conda install python==3.6.12
conda install snakemake==5.31.1
conda install fastqc==0.11.9
conda install multiqc==1.9
conda install trimmomatic==0.39-1
conda install hisat2==2.2.1
conda install samtools==1.7
conda install bamtools==2.5.1
conda install stringtie==2.1.4
```

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

### Creating HISAT2 Indexes
The HISAT2 indexes used here were created directly from the genome and annotation files downloaded from ENSEMBL. Creating HISAT2 indexes requires high memory and should be run on a high memory node on the HPC cluster. Successful index building will result in the creation of a folder *{index_name}* containing 8 index files with the names *{index_name}.{1-8}.ht2*

```bash
# Create splice-aware genome indexes that incorporate known splice sites and exons from GTF annotations
# Mouse
hisat2_extract_splice_sites.py mm10_annotation.gtf > mm10_splice_sites.ss
hisat2_extract_exons.py mm10_annotation.gtf > mm10_exons.exon
hisat2-build -p 16 --exon mm10_exons.exon --ss mm10_splice_sites.ss /mm10_genome.fa mm10_genome_tran

# Rat
hisat2_extract_splice_sites.py Rn6_annotation.gtf > Rn6_splice_sites.ss
hisat2_extract_exons.py Rn6_annotation.gtf > Rn6_exons.exon
hisat2-build -p 16 --exon Rn6_exons.exon --ss Rn6_splice_sites.ss Rn6_genome.fa Rn6_genome_tran
```

### Workflow, Script Order, and Output Directory Structure
The RNA-seq workflow will create directories for each dataset that follow the structures outlined in their respective .txt files. Because Snakemake creates directories as it runs, the only directories that must be created in advance are the reference directories and the following:
- carpenter/01_fastq
- powell/01_fastq/by_run
- walker/01_fastq/by_run
- walker/01_fastq/concatenated

Before running the Snakemake scripts, place FASTQ files into the appropriate directory as follows:

Raw FASTQ files should be named according to the following conventions: \
Samples with 1 run: **{paper}\_{SRA_Sample_ID}\_{Treatment_Group}.fq** \
Samples with multiple runs: **{paper}\_{SRA_Sample_ID}\_{Treatment_Group}_run{#}.fq**

Carpenter FASTQ files should be located in the **carpenter/01_fastq** directory. Powell and Walker have multiple FASTQ files per sample, and should be located in the **powell/01_fastq/by_run** and **walker/01_fastq/by_run** directories, respectively.

Note: One sample (SRS2926581_S1no) from the Walker experiment was produced in a pilot experiment and contains only one run. Its FASTQ file should be renamed **walker_SRS2926581_S1no_*cat*.fq** and placed in the **walker/01_fastq/concatenated** directory.

For each dataset, run scripts in the following order:
1. **{paper}_01_fastq_02_trimmed_fastq.py** - Creates files in the 01_fastq and 02_trimmed_fastq directories
2. **{paper}_03_alignment_04_processing.py** - Creates files in the 03_alignment and 04_processing directories
3. **{paper}_05_counts.py** - Creates files in the 05_counts directory

## Activity Log
### 02/05/2021
Uploaded the following files:
- all/**all_manifest.xlsx**
- carpenter/**carpenter_01_fastq_02_trimmed_fastq.py**
- carpenter/**carpenter_03_alignment_04_processing.py**
- carpenter/**carpenter_05_counts.py**
- powell/**powell_01_fastq_02_trimmed_fastq.py**
- powell/**powell_03_alignment_04_processing.py**
- powell/**powell_05_counts.py**
- walker/**walker_01_fastq_02_trimmed_fastq.py**
- walker/**walker_03_alignment_04_processing.py**
- walker/**walker_05_counts.py**
- Several junk files

### 02/06/2021
Updated README:
- Added conda environment/dependencies
- Filled out dataset descriptions
- Added links for reference files
- Added info for using scripts on the HPC cluster and building HISAT2 indexes
- Added GitHub directory structure

Uploaded/created files:
- all/**RNASeqWorkflow02062021.pdf** (Workflow image)
- carpenter/**carpenter_directory_structure.txt**
- powell/**powell_directory_structure.txt**
- walker/**walker_directory_strucutre.txt**

### 02/12/2021
Uploaded files:
- walker/**walker_file_sizes.txt**
- powell/**powell_file_sizes.txt**

Includes the original file sizes for all Walker & Powell files, so they can be double-checked if there are errors in the workflow.

### 02/15/2021
Updated all Snakemake and JSON files:
- Included more information on how to run the scripts, including how to set up the initial directory structure
- File paths are now more customizable to the user and their starting directory
- File paths in JSON files changed to reflect alterations to Snakemake scripts

Uploaded files:
- carpenter/**carpenter_file_sizes.txt**

## To-Do List
### General
- Change trimmomatic parameters
- Change StringTie parameters
- Finish filling out all_manifest.xlsx, and then split into experiment-specific .csv files
- Change extension of .py files to be .snk or .snakefile
- Create differential expression scripts