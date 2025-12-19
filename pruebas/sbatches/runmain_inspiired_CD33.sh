#!/bin/bash
#SBATCH --job-name=main_running
#SBATCH --output=sbatches/logs/slurm-%j.out  # Archivo que slurm generara con el output del trabajo
#SBATCH --error=sbatches/logs/slurm-%j.err   # Archivo que slurm generara con los errores
#SBATCH --ntasks=1             # Tareas independientes que se van a ejecutar
#SBATCH --cpus-per-task=15      # Numero de CPUs por tarea
#SBATCH --mem=100G              # Memoria RAM necesaria
#SBATCH --time=10:00:00        # Tiempo que tu trabajo podra estar ejecutandose
#SBATCH --partition=short      # Particion del cluster
#SBATCH --mail-type=END,FAIL   # Slurm puede mandar un mail cuando el trabajo acabe
#SBATCH --mail-user=lrenteria@external.unav.es

cd /beegfs/home/lrenteria/inspiired_nf


module load Nextflow/24.04.2
module load Python


nextflow run main_inspiired_CD33.nf \
    --runfolderDir /home/lrenteria/inspiired_nf/221122_M05935_0264_000000000-KJVKL/ \
    --outputDir /home/lrenteria/inspiired_nf/output_221122/ \
    --samplesheet /home/lrenteria/inspiired_nf/221122_M05935_0264_000000000-KJVKL/SampleSheet.csv \
    --projectName CD33 \
    -resume
