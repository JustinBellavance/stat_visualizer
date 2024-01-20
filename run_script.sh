#!/bin/bash
#SBATCH --account=def-gsarah
#SBATCH --time=23:58:00
#SBATCH --mem=20G
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1
#SBATCH --nodes=1

module load python/3.11.5
python3 tests.py