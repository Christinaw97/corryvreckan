# Corryvreckan

This is a fork of the main Corryvreckan software: https://gitlab.cern.ch/corryvreckan/corryvreckan/commits/master

## Installation (should work on lpc and lxplus el9 systems)

```bash
git clone https://github.com/eudaq/eudaq.git
git clone git@github.com:Christinaw97/corryvreckan.git
source corryvreckan/etc/setup_lxplus.sh

# install eudaq
cd eudaq
mkdir build && cd build
cmake .. -DEUDAQ_BUILD_EXECUTABLE=OFF -DEUDAQ_BUILD_GUI=OFF -DUSER_TLU_BUILD=ON 
cmake ..
make install -j8

## install corryvreckan
cd ../../corryvreckan
mkdir build && cd build
cmake .. -DBUILD_EventLoaderEUDAQ2=ON
export eudaq_DIR=/path/to/eudaq #use absolute path here
cmake ..
make install -j8
```

Run `./testbeam_script.sh [run number]` for test beam reconstruction

### Uninstallation
To uninstall and remove compilation files (e.g. when you want to do a clean recompilation/reinstallation):

```
$ rm -r build/ lib/ bin/
```
