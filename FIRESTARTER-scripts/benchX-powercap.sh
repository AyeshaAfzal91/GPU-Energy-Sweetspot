#!/bin/bash

set -eu

if [ $# -ne 5 ]; then
    echo "Usage: ${0} [result-dir] [freq-start] [freq-end] [freq-reset] [freq-step]"
    exit 1
fi

# The reason POWERCAP_RESET differs from POWERCAP_END is because of H100 & H200.
# On our cluster they don't run on the default power limit of 700W, but 500W instead.
# So allow the caller to reset the powerlimit to a custom one after benching.
POWERCAP_START="${2}"
POWERCAP_END="${3}"
POWERCAP_RESET="${4}"
POWERCAP_STEP="${5}"

function reset_powercap() {
    echo "Resetting clocks..."
    sudo nvidia-smi --power-limit="${POWERCAP_RESET}"
    echo "Done!"
}

trap 'reset_powercap' 15

RESULT_DIR="$(readlink -f "${1}")"
rm -fr "${RESULT_DIR}"
mkdir -p "${RESULT_DIR}"

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

for ((POWERCAP=POWERCAP_START; POWERCAP<=POWERCAP_END; POWERCAP+=POWERCAP_STEP)); do
    if (( POWERCAP > POWERCAP_END )); then
        POWERCAP=${POWERCAP_END}
    fi

    echo "SETTING GPU POWERCAP TO GRAPHICS=${POWERCAP}"
    sudo /usr/bin/nvidia-smi --power-limit="${POWERCAP}"

    SUB_RESULT_DIR="${RESULT_DIR}/${POWERCAP}"
    mkdir -p "${SUB_RESULT_DIR}"

    NVLOG_FILE="${SUB_RESULT_DIR}/${POWERCAP}_powerlog.csv"
    FSLOG_FILE="${SUB_RESULT_DIR}/${POWERCAP}_pid-$$.fsout"

    echo "Nvidia Log: ${NVLOG_FILE}"
    echo "FIRESTARTER Log: ${FSLOG_FILE}"

    nvidia-smi -i "${CUDA_VISIBLE_DEVICES}" --loop-ms=100 --query-gpu=timestamp,power.draw,utilization.gpu --format=csv > "${NVLOG_FILE}" & POWER_PID=$!
    "${BINARY}" -n 1 -l 1 -t 60 --period 1000000 > "${FSLOG_FILE}"
    kill "${POWER_PID}"

    echo "End: `date`  / `date +%s`"
done

reset_powercap
