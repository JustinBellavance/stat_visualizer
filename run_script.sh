#!/bin/bash
#SBATCH --account=def-gsarah
#SBATCH --time=04:00:00
#SBATCH --mem=200G

module load python/3.11.5
python3 tests.py