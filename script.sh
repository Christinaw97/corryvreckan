#!/bin/bash
source etc/setup_lxplus.sh


BASE_PATH=2025_08_SNSPD
config=mimosa_v1.conf
run_number=$(printf "%06d" "$1")

cp /eos/uscms/store/group/cmstestbeam/${BASE_PATH}/Tracks/RawData/run${run_number}.raw input.raw

bin/corry -c ${config}
cp output/TreeWriter/output.root /eos/uscms/store/user/cmstestbeam/${BASE_PATH}/Tracks/RecoData/v1/Run${1}_CMSTiming_FastTriggerStream_converted.root
rm -rf output
