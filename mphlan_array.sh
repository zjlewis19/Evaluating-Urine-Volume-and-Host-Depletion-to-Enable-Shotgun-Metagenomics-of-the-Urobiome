#!/bin/bash
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --account=PAS1331
#SBATCH -J mphlan_array
#SBATCH --output=/fs/ess/PAS1331/ZJL/VHale001_091323_update/mphlan_logs/%x_%A_%a
#SBATCH --array=1-68%10
#SBATCH --cpus-per-task=16
#SBATCH --nodes=1
#SBATCH --mem 117000mb
#SBATCH --time=04:00:00

cd /fs/ess/PAS1331/ZJL/VHale001_091323_update

module load bowtie2
PATH=$PATH:~/.local/bin

sample=$(sed "${SLURM_ARRAY_TASK_ID}q;d" ./samples.config)

R1=${sample}_R1_001.fastq.gz
R2=${sample}_R2_001.fastq.gz
metaphlan --bowtie2out mphlan_array_bowtie/${sample}_mphlan_bowtie_output.bz2 --input_type fastq -t rel_ab_w_read_stats reads/${R1},reads/${R2} mphlan_array_outputs/${sample}_metaphlan_output.txt 