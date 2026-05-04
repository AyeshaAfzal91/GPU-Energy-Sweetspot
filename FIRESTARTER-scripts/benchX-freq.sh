#!/bin/bash

set -eu

if [ $# -ne 4 ]; then
    echo "Usage: ${0} [result-dir] [freq-start] [freq-end] [freq-step]"
    exit 1
fi

function reset_clocks() {
    echo "Resetting clocks..."
    sudo nvidia-smi -rgc
    echo "Done!"
}

trap 'reset_clocks' 15

RESULT_DIR="$(readlink -f "${1}")"
rm -fr "${RESULT_DIR}"
mkdir -p "${RESULT_DIR}"

FREQ_START="${2}"
FREQ_END="${3}"
FREQ_STEP="${4}"

LOG_FILE="${RESULT_DIR}/firestarter-bench-`date +%Y%m%d`-`hostname -s`.$$.log"
echo "Saving bench log also to log file ..."
echo "Please be patient ..."
exec > >( tee "${LOG_FILE}" ) 2>&1

# save script for reference :)
echo "### This is ${0}"
cat "${0}"

echo "Start: `date`  / `date +%s`"

echo "########################################"

BINARY=/home/hpc/unrz/unrz104h/Projects/FIRESTARTER/sources/build/src/FIRESTARTER_CUDA

#cat /proc/cpuinfo

module load cuda/12.9.0

export CUDA_VISIBLE_DEVICES=0
nvidia-smi -i ${CUDA_VISIBLE_DEVICES}

for ((FREQ=FREQ_START; FREQ<=FREQ_END; FREQ+=FREQ_STEP)); do
    if (( FREQ > FREQ_END )); then
        FREQ=${FREQ_END}
    fi
    # === PARSE GPU CLOCK SETTINGS ===
    GPU_GRAPHICS_CLOCK="${FREQ}"

    # === SET GPU CLOCKS ===  
    echo "SETTING GPU CLOCKS TO GRAPHICS=${GPU_GRAPHICS_CLOCK}"
    sudo /usr/bin/nvidia-smi --lock-gpu-clocks="${GPU_GRAPHICS_CLOCK},${GPU_GRAPHICS_CLOCK}"

    SUB_RESULT_DIR="${RESULT_DIR}/${GPU_GRAPHICS_CLOCK}"
    mkdir -p "${SUB_RESULT_DIR}"

    NVLOG_FILE="${SUB_RESULT_DIR}/${GPU_GRAPHICS_CLOCK}_powerlog.csv"
    FSLOG_FILE="${SUB_RESULT_DIR}/${GPU_GRAPHICS_CLOCK}_pid-$$.fsout"

    echo "Nvidia Log: ${NVLOG_FILE}"
    echo "FIRESTARTER Log: ${FSLOG_FILE}"

    nvidia-smi -i "${CUDA_VISIBLE_DEVICES}" --loop-ms=100 --query-gpu=timestamp,power.draw,utilization.gpu --format=csv > "${NVLOG_FILE}" & POWER_PID=$!
    "${BINARY}" -n 1 -l 1 -t 60 --period 1000000 > "${FSLOG_FILE}"
    kill "${POWER_PID}"

    echo "End: `date`  / `date +%s`"
done

reset_clocks
