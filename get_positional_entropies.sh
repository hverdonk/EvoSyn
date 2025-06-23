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
#SBATCH --time=00:15:00
#
# Job name
#SBATCH --job-name=getPositionalEntropies
#
# Output file
#SBATCH --output=/home/tuo90294/slurm/out/posEntropies-flanked-conserved-seqs-%j.out
#SBATCH --error=/home/tuo90294/slurm/err/posEntropies-flanked-conserved-seqs-%j.err
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
    python3 -m scripts.get_positional_entropies \
    --input-fasta /home/tuo90294/genes/20250606_conserved_seqs_flanked.fasta \
    --output-tsv /home/tuo90294/20250606_conserved_seqs_flanked_positionalEntropies_evo-1.5-8k-base.tsv \
    --model-name 'evo-1.5-8k-base' \
    --device cuda:0

