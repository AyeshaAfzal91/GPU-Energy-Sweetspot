#!/bin/bash -l
#SBATCH --job-name=benchH200-freq
#SBATCH --nodes=1
#SBATCH --gres=gpu:h200:1
#SBATCH --time=24:00:00
#SBATCH --partition=h200

set -eu

SCRIPT_DIR="FIRESTARTER-scripts"

exec "${SCRIPT_DIR}/benchX-freq.sh" "${SCRIPT_DIR}/../FIRESTARTER-H200-frequency" \
    345 1980 15
