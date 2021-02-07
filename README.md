# RodentAddiction

Annika Vannan, email: avannan@asu.edu \
https://github.com/SexChrLab/RodentAddiction \
Last modified: 02/05/2021

## Project Description
The purpose of this project is to identify genes relevant to drug-motivated behavior in rodent models of cocaine addiction. RNA-seq datasets were obtained from the papers listed below. Data from each paper/experiment are referred to by their first author (e.g. Carpenter data) throughout this README and other files in the directory.

Carpenter, M.D., Hu, Q., Bond, A.M. *et al.* Nr4a1 suppresses cocaine-induced behavior via epigenetic regulation of homeostatic target genes. *Nat Commun* **11,** 504 (2020). https://doi.org/10.1038/s41467-020-14331-y

Powell, G.L., Vannan, A., Bastle, R.M. *et al.* Environmental enrichment during forced abstinence from cocaine self-administration opposes gene network expression changes associated with the incubation effect. *Sci Rep* **10,** 11291 (2020). https://doi.org/10.1038/s41598-020-67966-8

Walker, D.M., Cates, H.M, Loh, Y.E. *et al.* Cocaine self-administration alters transcriptome-wide responses in the brain's reward circuitry. *Biol Psychiatry* **15**, 84(12):867-880. (2018) https://doi.org/10.1016/j.biopsych.2018.04.009

## Directory Structure
### GitHub Directory
The folder **all/** contains scripts and other files pertaining to the overall project. There are also separate folders for each the Carpenter, Powell, and Walker experiments that contain experiment-specific files and scripts. The folder **junk/** contains junk files and scripts that I'm not ready to part with yet.

### Workflow Directory
The RNA-seq workflow will create a directory that follows the structure outlined in the file **ADD FILE**

## Data Overview & Download

**Carpenter** \
Experiment Description: \
Species: *Mus musculus* (mouse) \
Strain: C57BL/6J \
Sex: Male \
Tissue: Nucleus accumbens (whole) \
Treatment Groups:
- Cocaine, 1d abstinence (C1abs)
- Cocaine, 28d abstinence (C28abs)
- Saline, 1d abstinence (S1abs)
- Saline, 28d abstinence (S21abs)

**Powell** \
Experiment Description: \
Species: *Rattus norvegicus* (rat) \
Strain: Sprague-Dawley \
Sex: Male \
Tissue: Nucleus accumbens shell\
Treatment Groups:
- Cocaine, 1d abstinence in isolated housing (CI1cue)
- Cocaine, 21d abstinence in isolated housing (CI21cue)
- Cocaine, 1d abstinence in enriched housing (CE1cue)
- Cocaine, 21d abstinence in enriched housing (CE21cue)

**Walker** \
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
Mouse genome: \
Mouse genome annotation: \
Rat genome: \
Rat genome annotation: \

## Using These Scripts
Conda environment and dependencies

## Activity Log
### 02/06/2021
Uploaded the following files:
- all/**all_manifest.xlsx**
- carpenter/scripts/**carpenter_01_fastq_02_trimmed_fastq.py**
- carpenter/scripts/**carpenter_03_alignment_04_processing.py**
- carpenter/scripts/**carpenter_05_counts.py**
- powell/scripts/**powell_01_fastq_02_trimmed_fastq.py**
- powell/scripts/**powell_03_alignment_04_processing.py**
- powell/scripts/**powell_05_counts.py**
- walker/scripts/**walker_01_fastq_02_trimmed_fastq.py**
- walker/scripts/**walker_03_alignment_04_processing.py**
- walker/scripts/**walker_05_counts.py**
- Several junk files

## To-Do List
### General
- Change trimmomatic parameters
- Change StringTie parameters
- Finish filling out all_manifest.xlsx, and then split into experiment-specific .csv files
- Change extension of .py files to be .snk or .snakefile
- Upload directory structure document
- Fill out README! :)
  - Add project description
  - Add rodent strain
  - Conda environment/dependencies
  - Add reference file links
  - Fill out data description
  
  
  
