#!/bin/bash
#SBATCH --job-name=gtex_subset          # Job name
#SBATCH --mail-user=avannan@asu.edu    # Send-to address
#SBATCH -n 1                           # Number of tasks per node (N)
#SBATCH -o slurm.%j.out                # STDOUT (%j = JobId)
#SBATCH -e slurm.%j.err                # STDERR (%j = JobId)
#SBATCH --qos=normal                   # Quality of service (priority)
#SBATCH -t 0-2:00                      # Maximum wall time (D-HH:MM)

###
### This script gets a subset of genes from the GTEX brain files
###

# Go to directory
cd /scratch/avannan/TTR/counts/

# Set brain tissue files
brain_tissues=(Brain-Amygdala Brain-Anteriorcingulatecortex_BA24 Brain-Caudate_basalganglia Brain-CerebellarHemisphere \
Brain-Cerebellum Brain-Cortex Brain-FrontalCortex_BA9 Brain-Hippocampus Brain-Hypothalamus Brain-Nucleusaccumbens_basalganglia \
Brain-Putamen_basalganglia Brain-Spinalcord_cervicalc-1 Brain-Substantianigra Pituitary)

# Subset brain tissue files by genes of interest
for tissue in ${brain_tissues[@]}; do
  echo "START Subsetting genes for tissue: ${tissue}"
  awk 'NR==1; \
  # Craving Genes
  $1 ~ "ENSG00000067715"; $1 ~ "ENSG00000077235"; $1 ~ "ENSG00000102189"; $1 ~ "ENSG00000103197"; \
  $1 ~ "ENSG00000105058"; $1 ~ "ENSG00000116251"; $1 ~ "ENSG00000116353"; $1 ~ "ENSG00000121858"; \
  $1 ~ "ENSG00000122176"; $1 ~ "ENSG00000138030"; $1 ~ "ENSG00000141627"; $1 ~ "ENSG00000141750"; \
  $1 ~ "ENSG00000145982"; $1 ~ "ENSG00000146072"; $1 ~ "ENSG00000154144"; $1 ~ "ENSG00000159579"; \
  $1 ~ "ENSG00000166261"; $1 ~ "ENSG00000173175"; $1 ~ "ENSG00000196935"; \
  # GWAS Genes
  $1 ~ "ENSG00000189319"; $1 ~ "ENSG00000112078"; $1 ~ "ENSG00000144119"; $1 ~ "ENSG00000112079"; \
  $1 ~ "ENSG00000152056"; $1 ~ "ENSG00000164344"; $1 ~ "ENSG00000169291"; $1 ~ "ENSG00000257941"; \
  $1 ~ "ENSG00000271179"; \
  # Other Genes of Interest
  $1 ~ "ENSG00000118260"; $1 ~ "ENSG00000162374"; $1 ~ "ENSG00000123358"; $1 ~ "ENSG00000198576";
  $1 ~ "ENSG00000176697"\
  ${tissue}_counts.txt > subset_${tissue}.txt \
  echo "FINISH subsetting genes for tissue: ${tissue}"
done
