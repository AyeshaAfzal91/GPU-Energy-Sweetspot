#!/bin/bash -l
#SBATCH --job-name=benchH200-powercap
#SBATCH --nodes=1
#SBATCH --gres=gpu:h200:1
#SBATCH --time=24:00:00
#SBATCH --partition=h200

set -eu

SCRIPT_DIR="FIRESTARTER-scripts"

exec "${SCRIPT_DIR}/benchX-powercap.sh" "${SCRIPT_DIR}/../FIRESTARTER-H200-powercap" \
    200 700 500 10
