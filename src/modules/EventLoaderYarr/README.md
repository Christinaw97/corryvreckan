---
# SPDX-FileCopyrightText: 2017-2024 CERN and the Corryvreckan authors
# SPDX-License-Identifier: CC-BY-4.0 OR MIT
---
# EventLoaderYarr
**Maintainer**: Pranav Manoj (mpranav579@gmail.com), Luc Le Pottier (luclepot@lbl.gov), Timon Heim (theim@lbl.gov)
**Module Type**: *DETECTOR* 
**Status**: Functional

### Description
This module loads data from YARR format binary files adding it to the clipboard and plotting the pixel hit positions. The module can read data from multiple detectors at the same time. Each detector must have a separate data file with all files located in a common directory. Each datafile name must contain the the detector name (as provided in the detector configuration file) as well as the extension ".raw". For example, if a detector is named "0x20c52", a valid YARR datafile name would be "datafile_0x20c52.raw". The event loader supports two timestamp settings: Trigger tag timing makes use of the timestamp in the data for each event and represents the number of milliseconds since 12AM on the day the run began. If Trigger tag timing is not used, an arbitrary timestamp based on the event number is assigned (25ns time resolution of each trigger window is still preserved). 

### Parameters
* `input_directory`: Path of the directory above the data files.
* `log_level`: Specifies the lowest log level that should be reported. Possible values are FATAL, STATUS, ERROR, WARNING, INFO, and DEBUG, where all options are case-insensitive. Defaults to the INFO level.
* `trigger_tag_timing`: Toggle to specify timing setting. Defaults to False.

### Plots produced
For each detector the following plots are produced:

* 2D histogram of pixel hit positions
* 1D histogram of number of hits vs timestamp (only made if trigger_tag_timing is on)
* 1D histogram of number of events vs timestamp (only made if trigger_tag_timing is on)

### Usage
```toml
[EventLoaderYarr]
input_directory = "path_to_Yarr_data_directory/" (ex. "home/usr/path/to/data/directory")
log_level = "INFO"
trigger_tag_timing = True
```
