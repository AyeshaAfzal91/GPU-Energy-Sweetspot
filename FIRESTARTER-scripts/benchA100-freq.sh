#!/bin/bash -l
#SBATCH --job-name=benchA100-freq
#SBATCH --nodes=1
#SBATCH --gres=gpu:a100:1
#SBATCH --constraint=a100_40
#SBATCH --time=24:00:00

set -eu

SCRIPT_DIR="FIRESTARTER-scripts"

exec "${SCRIPT_DIR}/benchX-freq.sh" "${SCRIPT_DIR}/../FIRESTARTER-A100-frequency" \
    210 1410 15
