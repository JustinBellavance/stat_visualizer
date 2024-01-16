#!/bin/bash
#SBATCH --account=def-gsarah
#SBATCH --time=06:00:00
#SBATCH --mem=100G
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1

module load python/3.11.5
python3 tests.py