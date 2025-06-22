#!/bin/bash

#SBATCH --job-name=tcga_analysis
#SBATCH --partition=cpu
#SBATCH --ntasks=1
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4
#SBATCH --time=24:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=k24058218@kcl.ac.uk

# Define paths
CONTAINER="/scratch/prj/bmb_phyl_chron_disease/esta_proj/bioconductor_latest.sif"
WORKDIR="/scratch/prj/bmb_phyl_chron_disease/esta_proj/"
SCRIPT="TCGA_DE_masterscript.R"

# Create output directories (if they don't exist)
mkdir -p ${WORKDIR}/DE_results
mkdir -p ${WORKDIR}/plots
mkdir -p ${WORKDIR}/stats
mkdir -p ${WORKDIR}/logs

# Run the analysis in the container
singularity exec \
  --bind ${WORKDIR}:/workdir \
  ${CONTAINER} \
  Rscript /workdir/${SCRIPT}
