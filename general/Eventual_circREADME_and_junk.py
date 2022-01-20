#!/bin/bash
#SBATCH --job-name=rat_hisat2_build  # Job name
#SBATCH --mail-type=ALL                # Send a notification when the job starts, stops, or fails
#SBATCH --mail-user=avannan@asu.edu    # Send-to address
#SBATCH -o slurm.%j.out                # STDOUT (%j = JobId 
#SBATCH -e slurm.%j.err                # STDERR (%j = JobId)
#SBATCH --qos=normal                   # Quality of service (priority)
#SBATCH -t 0-02:00                     # Maximum wall time (D-HH:MM)
#SBATCH -n 8

source activate rodent_addiction
#hisat2-build -f /scratch/avannan/refs/mouse/GRCm39/Mus_musculus.GRCm39.cdna.all.fa /scratch/avannan/refs/mouse/GRCm39/HISAT2/GRCm39_cdna_hisat2 # Mouse
#hisat2-build -f /scratch/avannan/refs/rat/mRatBN7.2/Rattus_norvegicus.mRatBN7.2.cdna.all.fa /scratch/avannan/refs/rat/mRatBN7.2/HISAT2/mRatBN7.2_cdna_hisat2 # Rat
#hisat2-build -f /scratch/avannan/refs/rat/Rn6/Rattus_norvegicus.Rnor_6.0.cdna.all.fa /scratch/avannan/refs/rat/Rn6/HISAT2/Rn6_cdna_hisat2 # Rat
hisat2-build -f /scratch/avannan/refs/rat/mRatBN7.2/Rattus_norvegicus.mRatBN7.2.dna.toplevel.fa /scratch/avannan/refs/rat/mRatBN7.2/HISAT2/mRatBN7.2_dna_hisat2 # Rat
