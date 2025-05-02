#!/bin/bash
#======================================================
#
# Job script for running a parallel job on a single gpu
#
# Submit as follows:
#
#======================================================

#======================================================
# Propagate environment variables to the compute node
#SBATCH --export=ALL
#
# Run in the gpu partition (queue)
#SBATCH --partition=gpu
#
# Total number GPUs for the job
#SBATCH --gpus=1
#
# Number of GPUs to use per node (max 2)
#SBATCH --gpus-per-node=1
#
# Number of CPUs per GPU
#SBATCH --cpus-per-gpu=1
#
# Specify (hard) runtime (HH:MM:SS)
#SBATCH --time=00:30:00
#
# Job name
#SBATCH --job-name=score-Dmel-gene
#
# Output file
#SBATCH --output=/home/tuo90294/slurm/out/score-Dmel-gene-%j.out
#SBATCH --error=/home/tuo90294/slurm/err/score-Dmel-gene-%j.err
#======================================================

# Load CUDA always
module load cuda

# Activate conda environment always
export PATH=/home/tuo90294/miniconda3/envs/evo1/bin:$PATH
# source activate torch-2.5.1

# change to directory where 'sbatch' was called
cd /home/tuo90294/evo

# runs on GPU 0
srun --gpus 1 \
    python3 -m scripts.score \
    --input-fasta /home/tuo90294/FBgn0037797_NT_033777.3_Dmel_allSynonymousMutants.fasta \
    --output-tsv /home/tuo90294/FBgn0037797_NT_033777.3_Dmel_allSynonymousMutants_scores.tsv \
    --model-name 'evo-1-131k-base' \
    --device cuda:0

