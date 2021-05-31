#!/bin/tcsh
#
##SBATCH -p common                # Partition to submit to (comma separated)
#SBATCH -J align		         # Job name
#SBATCH -n 8                     # Number of cores
#SBATCH -N 1                     # Ensure that all cores are on one machine
#SBATCH -t 2-00:00                # Runtime in D-HH:MM (or use minutes)
#SBATCH --mem 30000               # Memory in MB
#SBATCH -o _ExtractVAR_%j.out # File for STDOUT (with jobid = %j) 
#SBATCH -e ExtractVAR_%j.err       # File for STDERR (with jobid = %j) 
#SBATCH --mail-type=ALL          # Type of email notification: BEGIN,END,FAIL,ALL
#SBATCH --mail-user=jws48@duke.edu  # Email where notifications will be sent
#Your actual work goes after this line

module load samtools/1.10-gcb01
module load STAR/2.7.2b-gcb01


./extractVAR.pl /gpfs/fs1/data/taylorlab/Genomes/3D7_ASM276v2/FASTA \
/gpfs/fs1/data/taylorlab/HbAS/maliHbAS/VAR/out \
/gpfs/fs1/data/taylorlab/HbAS/maliHbAS/align/out/HsUnmappedReads/1 \
/gpfs/fs1/data/taylorlab/HbAS/maliHbAS/align/out/HsUnmappedReads/2 \
/gpfs/fs1/data/taylorlab/Genomes/3D7_ASM276v2/GTF/Plasmodium_falciparum.ASM276v2.45.gtf
