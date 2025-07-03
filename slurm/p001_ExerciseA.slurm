#!/bin/bash
#SBATCH --job-name=summer_school
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32           
#SBATCH --time=17-00:00:00
#SBATCH --mem=20G
#SBATCH --partition=gpu
#SBATCH --output=/dev/null

EXEC=/share/home/lsy_luzhen/summer_school_gwas_exercise_2025/src/summer_school_gwas_exercise_2025/ExerciseA.sh
SCRIPT_NAME=$(basename $EXEC .sh)
TIMESTAMP=$(date +"%Y%m%d_%H%M%S")
OUTPUT_FILE="${SCRIPT_NAME}_${TIMESTAMP}.out"
exec > $OUTPUT_FILE 2>&1

echo 'start at:'
date
bash $EXEC
echo 'end at:'
date
