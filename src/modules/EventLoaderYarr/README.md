---
# SPDX-FileCopyrightText: 2017-2024 CERN and the Corryvreckan authors
# SPDX-License-Identifier: CC-BY-4.0 OR MIT
---
# EventLoaderYarr
**Maintainer**: pmanoj (pigwhale.dhcp.lbl.gov)
**Module Type**: *DETECTOR* 
**Status**: Work in progress

### Description
This module loads data from YARR format binary files adding it to the clipboard and plotting the pixel hit positions. The module can read data from multiple detectors at the same time. Each detector must have a separate data file with all files located in a common directory. Each datafile name must contain the the detector name (as provided in the detector configuration file) as well as the extension ".raw". For example, if a detector is named "0x20c52", a valid YARR datafile name would be "datafile_0x20c52.raw". Trigger multiplicity should be provided in config file. 

### Parameters
* `input_directory`: Path of the directory above the data files.
* `trigger_multiplicity`: Trigger multiplicity or event length in number of triggers. Defaults to 17.
* `log_level`: Specifies the lowest log level that should be reported. Possible values are FATAL, STATUS, ERROR, WARNING, INFO, and DEBUG, where all options are case-insensitive. Defaults to the INFO level.

### Plots produced
For each detector the following plots are produced:

* 2D histogram of pixel hit positions
* 1D histogram of number of hits vs timestamp

### Usage
```toml
[EventLoaderYarr]
input_directory = "data/000699_bella_exttrigger"
trigger_multiplicity = 16
log_level = "INFO"

```
