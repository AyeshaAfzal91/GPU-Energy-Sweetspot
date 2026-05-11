#!/bin/sh

set -eu

batch() {
    ssh "${1}" << EOF
    cd "$(readlink -f "$(dirname "${0}")/..")"
    sbatch "./FIRESTARTER-scripts/${2}"
EOF
}

batch alex benchA40-freq.sh
batch alex benchA100-freq.sh
batch helma benchH100-freq.sh
batch helma benchH200-freq.sh
batch alex benchA40-powercap.sh
batch alex benchA100-powercap.sh
batch helma benchH100-powercap.sh
batch helma benchH200-powercap.sh
