#!/bin/bash
source etc/setup_lxplus.sh


BASE_PATH=2025_08_SNSPD
config=mimosa_v1.conf
run_number=$(printf "%06d" "$1")

#cp /eos/uscms/store/group/cmstestbeam/${BASE_PATH}/Tracks/RawData/run${run_number}.raw input.raw
#files="-o TreeWriter.file_name=output.root -o EventLoaderEUDAQ2.file_name=input.raw "
#bin/corry -c ${config} ${files} 

files="-o TreeWriter.file_name=output.root -o EventLoaderEUDAQ2.file_name=/eos/uscms/store/group/cmstestbeam/${BASE_PATH}/Tracks/RawData/run${run_number}.raw"
bin/corry -c ${config} ${files}

cp output/TreeWriter/output.root /eos/uscms/store/user/cmstestbeam/${BASE_PATH}/Tracks/RecoData/v1/Run${1}_CMSTiming_FastTriggerStream_converted.root

#make plots
mkdir -p plots_${BASE_PATH}
cd plot_scripts
./make_plots.sh ../output/out_MIMOSA.root Run${1} plots_${BASE_PATH}
cd ../

#rm ROOT files
rm -rf output
