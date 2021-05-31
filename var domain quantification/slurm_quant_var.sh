#!/bin/tcsh
#
##SBATCH -p common                # Partition to submit to (comma separated)
#SBATCH -J align		         # Job name
#SBATCH -n 8                     # Number of cores
#SBATCH -N 1                     # Ensure that all cores are on one machine
#SBATCH -t 2-00:00                # Runtime in D-HH:MM (or use minutes)
#SBATCH --mem 3000               # Memory in MB
#SBATCH -o _QuantVAR_k15_%j.out # File for STDOUT (with jobid = %j) 
#SBATCH -e QuantVAR_k15_%j.err       # File for STDERR (with jobid = %j) 
#SBATCH --mail-type=ALL          # Type of email notification: BEGIN,END,FAIL,ALL
#SBATCH --mail-user=jws48@duke.edu  # Email where notifications will be sent
#Your actual work goes after this line

./quant_var.pl /gpfs/fs1/data/taylorlab/Genomes/varDomains/PfEMP1domains_index_15 \
/gpfs/fs1/data/taylorlab/HbAS/maliHbAS/VAR/var_quant_k15 \
/gpfs/fs1/data/taylorlab/HbAS/maliHbAS/VAR/out/varReads/singleton
