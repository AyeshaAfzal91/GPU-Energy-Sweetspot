#!/bin/bash -l
#SBATCH --job-name=benchH100-powercap
#SBATCH --nodes=1
#SBATCH --gres=gpu:h100:1
#SBATCH --time=24:00:00
#SBATCH --partition=h100

set -eu

SCRIPT_DIR="FIRESTARTER-scripts"

exec "${SCRIPT_DIR}/benchX-powercap.sh" "${SCRIPT_DIR}/../FIRESTARTER-H100-powercap" \
    200 700 500 10
