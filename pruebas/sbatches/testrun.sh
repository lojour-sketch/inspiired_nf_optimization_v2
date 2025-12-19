#!/bin/bash
#SBATCH --job-name=testrun
#SBATCH --output=sbatches/logs/slurm-%j.out  # Archivo que slurm generara con el output del trabajo
#SBATCH --error=sbatches/logs/slurm-%j.err   # Archivo que slurm generara con los errores
#SBATCH --ntasks=1             # Tareas independientes que se van a ejecutar
#SBATCH --cpus-per-task=1     # Numero de CPUs por tarea
#SBATCH --mem=1G              # Memoria RAM necesaria
#SBATCH --time=4:00:00        # Tiempo que tu trabajo podra estar ejecutandose
#SBATCH --partition=short      # Particion del cluster
#SBATCH --mail-type=END,FAIL   # Slurm puede mandar un mail cuando el trabajo acabe
#SBATCH --mail-user=lrenteria@external.unav.es



cd /home/lrenteria/inspiired_nf


module load Nextflow/24.04.2
module load UMI-tools

nextflow run testmain.nf --linkerdata "/home/lrenteria/inspiired_nf/linkerdata.csv"