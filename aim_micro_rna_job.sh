#!/bin/bash
#SBATCH --job-name=Micro_rna_job
#SBATCH --time=18:00:00
#SBATCH --mem=32000

export APPTAINER_CACHEDIR=/scratch/p309221/apptainer

apptainer exec --overlay /scratch/p309221/micro_rna:/mnt -B /scratch/p309221/micro_rna:/mnt srnatoolbox_latest.sif /bin/bash /mnt/inside_container_v1.sh





