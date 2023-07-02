#!/bin/bash
#
#SBATCH --job-name=sobolSTYd
#SBATCH --ntasks=256
#SBATCH --nodes=1
#SBATCH --partition=cpu256
#SBATCH --time=20-00:00:00
#SBATCH --mem=450G
#SBATCH --mail-type=BEGIN,TIME_LIMIT_50,END,FAIL,ALL
#SBATCH --mail-user=m.giraud@fz-juelich.de


cd ..

python3 sobol_phloem.py 11 dry ${SLURM_JOB_ID}
#zip AllAuxC1 AllAuxC1/*.txt

#sbatch --nodelist=node03 All_C1.sh
