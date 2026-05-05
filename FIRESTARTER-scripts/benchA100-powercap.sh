#!/bin/bash -l
#SBATCH --job-name=benchA100-powercap
#SBATCH --nodes=1
#SBATCH --gres=gpu:a100:1
#SBATCH --constraint=a100_40
#SBATCH --time=24:00:00

set -eu

SCRIPT_DIR="FIRESTARTER-scripts"

exec "${SCRIPT_DIR}/benchX-powercap.sh" "${SCRIPT_DIR}/../FIRESTARTER-A100-powercap" \
    100 400 400 10
