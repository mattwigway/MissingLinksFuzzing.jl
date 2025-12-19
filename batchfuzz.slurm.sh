#!/bin/bash

#SBATCH --time=10-00:00:00
#SBATCH --time-min=12:00:00
#SBATCH --mem-per-cpu=2G

# Julia LTS
module load julia/1.10.10
./fuzz --log fuzzlog%i.log --csv fuzzlog%i.csv --headless