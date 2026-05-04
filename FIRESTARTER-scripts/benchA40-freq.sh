#!/bin/bash -l
#SBATCH --job-name=benchA40-freq
#SBATCH --nodes=1
#SBATCH --gres=gpu:a40:1
#SBATCH --time=24:00:00

set -eu

SCRIPT_DIR="FIRESTARTER-scripts"

exec "${SCRIPT_DIR}/benchX-freq.sh" "${SCRIPT_DIR}/../FIRESTARTER-A40-frequency" \
    210 1740 15
