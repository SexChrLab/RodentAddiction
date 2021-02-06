# GOALS
1. HISAT2 alignment
	a. Align to transcriptomes?
	b. Can StringTie deal with transcriptome output from HISAT2?
2. Use FastQC to detect adapters for Carpenter, Powell, & Walker experiments
	a. Get info on FASTQ files for Powell experiment - were they pre-trimmed by the Genomics core?
		i. Re-email Greg
3. Test Snakefiles for Carpenter project
	a. Test Snakefiles for Walker project (single-end reads)
4. Re-email authors for FASTQ headers, reported sex, & sample barcodes for NAc samples + adapter file
5. Update directory structure Google doc
6. Incorporate FASTQ headers from Carpenter et al. into config file


# QUESTIONS
1. FASTQ header/batch effect issue with Carpenter
2. HISAT2 transcriptome

# STRINGTIE ERROR
# Please make sure the -G annotation file uses the same naming convention for the genome sequences.

snakemake --snakefile carpenter_01_fastq_02_trimmed_fastq.py -j 8 --rerun-incomplete --latency-wait 60 --force --nolock --cluster "sbatch -n 1 --mem=60000 -c 4"

--nolock
--latency-wait 60
--rerun-incomplete
--verbose

# Get run info
# CARPENTER
esearch -db sra -query PRJNA593775 | efetch -format runinfo > carpenter_SRA_info.csv
cat carpenter_SRA_info.csv | cut -f 1 -d , | grep SRR > carpenter_SRR_IDs.txt
cat carpenter_SRA_info.csv | cut -f 1,25 -d , | sed 's/,/\t/' | grep SRR > carpenter_SRR_SRS_IDs.tsv

# WALKER
esearch -db sra -query PRJNA433508 | efetch -format runinfo > walker_SRA_info.csv
cat walker_SRA_info.csv | cut -f 1 -d , | grep SRR > walker_SRR_IDs.txt
cat walker_SRA_info.csv | cut -f 1,25 -d , | sed 's/,/\t/' | grep SRR > walker_SRR_SRS_IDs.tsv



# HISAT2 INDEXES
# Get index for Rat genome - one of the following commands
wget --content-disposition https://genome-idx.s3.amazonaws.com/hisat/rn6_genome.tar.gz



# IDEAS FOR ADAPTER DETECTION
fastp https://github.com/OpenGene/fastp




# GENERAL NOTES
# Create a bioinformatics enviroment called rodent_addiction.
conda create -y --name rodent_addiction

# Activate the bioinformatics enviroment.
conda activate rodent_addiction

# Set up bioconda and conda-forge channels for installing bioinformatics tools.
conda config --add channels bioconda
conda config --add channels conda-forge

# Install packages.
conda install hisat2==2.2.1
conda install samtools==1.7
conda install bamtools==2.5.1
conda install sra-tools==2.10.9
conda install parallel==20201122
conda install entrez-direct==13.9
conda install fastqc==0.11.9
conda install multiqc==1.9
conda install trimmomatic==0.39-1
conda install snakemake==5.31.1
conda install picard==2.24.0
conda install stringtie==2.1.4

# Fix Perl library problem, so that Perl commands correctly access the local installation of Perl.
unset PERL5LIB
export PERL5LIB=$PERL5LIB:/home/avannan/miniconda3/envs/rodent_addiction/lib/site_perl/5.26.2

# Download latest rat genome (DNA), transcriptome (cDNA), and annotation file: Rnor_6.0, database version 102.6
# Genome
wget ftp://ftp.ensembl.org/pub/release-102/fasta/rattus_norvegicus/dna/Rattus_norvegicus.Rnor_6.0.dna.toplevel.fa.gz
gzip -d Rattus_norvegicus.Rnor_6.0.dna.toplevel.fa.gz
mv Rattus_norvegicus.Rnor_6.0.dna.toplevel.fa Rn6_genome.fa
# Transcriptome
wget ftp://ftp.ensembl.org/pub/release-102/fasta/rattus_norvegicus/cdna/Rattus_norvegicus.Rnor_6.0.cdna.all.fa.gz
gzip -d Rattus_norvegicus.Rnor_6.0.cdna.all.fa.gz
mv Rattus_norvegicus.Rnor_6.0.cdna.all.fa Rn6_transcriptome.fa
# Annotation
wget ftp://ftp.ensembl.org/pub/release-102/gtf/rattus_norvegicus/Rattus_norvegicus.Rnor_6.0.102.gtf.gz
gzip -d Rattus_norvegicus.Rnor_6.0.102.gtf.gz
mv Rattus_norvegicus.Rnor_6.0.102.gtf Rn6_annotation.gtf
# Build HGFM index of genome with transcripts.
hisat2_extract_splice_sites.py Rn6_annotation.gtf > Rn6_splice_sites.ss
hisat2_extract_exons.py Rn6_annotation.gtf > Rn6_exons.exon
hisat2-build -p 16 --exon Rn6_exons.exon --ss Rn6_splice_sites.ss Rn6_genome.fa Rn6_genome_tran
# Build HFM index of transcriptome.
hisat2-build -p 16 Rn6_transcriptome.fa hisat2_Rn6_transcriptome

# Most recent mouse build (GRCm39/mm39) does not have an available transcriptome available on ENSEMBL.
# Download previous mouse genome (DNA), transcriptome (cDNA), and annotation file: GRCm38.p6/mm10, database version 102.38
# Genome
wget ftp://ftp.ensembl.org/pub/release-102/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.toplevel.fa.gz
gzip -d Mus_musculus.GRCm38.dna.toplevel.fa.gz
mv Mus_musculus.GRCm38.dna.toplevel.fa mm10_genome.fa
# Transcriptome
wget ftp://ftp.ensembl.org/pub/release-102/fasta/mus_musculus/cdna/Mus_musculus.GRCm38.cdna.all.fa.gz
gzip -d Mus_musculus.GRCm38.cdna.all.fa.gz
mv Mus_musculus.GRCm38.cdna.all.fa mm10_transcriptome.fa
# Annotation
wget ftp://ftp.ensembl.org/pub/release-102/gtf/mus_musculus/Mus_musculus.GRCm38.102.gtf.gz
gzip -d Mus_musculus.GRCm38.102.gtf.gz
mv Mus_musculus.GRCm38.102.gtf mm10_annotation.gtf
# Build HGFM index of genome with transcripts.
hisat2_extract_splice_sites.py mm10_annotation.gtf > mm10_splice_sites.ss
hisat2_extract_exons.py mm10_annotation.gtf > mm10_exons.exon
hisat2-build -p 16 --exon /home/avannan/rodent_addiction/refs/mouse/HISAT2/mm10_exons.exon --ss /home/avannan/rodent_addiction/refs/mouse/HISAT2/mm10_splice_sites.ss /home/avannan/rodent_addiction/refs/mouse/mm10_genome.fa mm10_genome_tran
# Build HFM index of transcriptome.
hisat2-build -p 16 /home/avannan/rodent_addiction/refs/mouse/mm10_transcriptome.fa hisat2_mm10_transcriptome


# Download mm11 genome/transcriptome/etc. from GenCode
# Genome
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M11/GRCm38.p4.genome.fa.gz
gzip -d GRCm38.p4.genome.fa.gz
mv GRCm38.p4.genome.fa mm11_genome.fa
# Transcriptome
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M11/gencode.vM11.transcripts.fa.gz
gzip -d gencode.vM11.transcripts.fa.gz
mv gencode.vM11.transcripts.fa mm11_transcriptome.fa
# Annotation
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M11/gencode.vM11.annotation.gtf.gz
gzip -d gencode.vM11.annotation.gtf.gz
mv gencode.vM11.annotation.gtf mm11_annotation.gtf
# Build HGFM index of genome with transcripts.
hisat2_extract_splice_sites.py mm11_annotation.gtf > mm11_splice_sites.ss
hisat2_extract_exons.py mm11_annotation.gtf > mm11_exons.exon
hisat2-build -p 16 --exon /home/avannan/rodent_addiction/refs/mouse/HISAT2/mm11_exons.exon --ss /home/avannan/rodent_addiction/refs/mouse/HISAT2/mm11_splice_sites.ss /home/avannan/rodent_addiction/refs/mouse/mm11_genome.fa mm11_genome_tran
# Build HFM index of transcriptome.

"hisat2_rmdup_dir": "/scratch/avannan/rodent_addiction/carpenter/04_processing/HISAT2/picard/rmdup_bams/",
"hisat2_rmdup_metrics_dir": "/scratch/avannan/rodent_addiction/carpenter/04_processing/HISAT2/picard/metrics/rmdup_bams/",
"hisat2_rmdup_rg_dir": "/scratch/avannan/rodent_addiction/carpenter/04_processing/HISAT2/picard/rmdup_rg_bams/",
"hisat2_rmdup_rg_stats_dir": "/scratch/avannan/rodent_addiction/carpenter/04_processing/HISAT2/bamtools_stats/rmdup_rg_bams/",
"hisat2_align_dir": "/scratch/avannan/rodent_addiction/carpenter/04_processing/HISAT2/processed_bams/",