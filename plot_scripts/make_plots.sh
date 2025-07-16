#!/bin/bash

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
echo ${SCRIPT_DIR}
#python3 ${SCRIPT_DIR}/plot_intersect.py --input_file $1 --group Tracking4D --postfix $2
#python3 ${SCRIPT_DIR}/plot_intersect.py --input_file $1 --group Clustering4D --postfix $2
python3 ${SCRIPT_DIR}/plot_residual.py --input_file $1 --postfix $2
python3 ${SCRIPT_DIR}/plot_tracks.py --input_file $1 --postfix $2
