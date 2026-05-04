#!/bin/bash -l
#SBATCH --job-name=benchH100-freq
#SBATCH --nodes=1
#SBATCH --gres=gpu:h100:1
#SBATCH --time=24:00:00
#SBATCH --partition=h100

set -eu

SCRIPT_DIR="FIRESTARTER-scripts"

exec "${SCRIPT_DIR}/benchX-freq.sh" "${SCRIPT_DIR}/../FIRESTARTER-H100-frequency" \
    345 1980 15
