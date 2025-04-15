---
# SPDX-FileCopyrightText: 2023-2024 CERN and the Corryvreckan authors
# SPDX-License-Identifier: CC-BY-4.0 OR MIT
---
# EventLoaderHDF5
**Maintainer**: Christian Bespin (<cbespin@uni-bonn.de>)
**Module Type**: *DETECTOR*  
**Status**: Work in progress

### Description
This module loads data from hdf5 files and adds it to the clipboard. The input file must have extension `.h5` and follow the structure below:

| column |  row  |  raw  |  charge  |   timestamp \[ns\]   | trigger_number |
|:------:|:-----:|:-----:|:--------:|:--------------------:|:--------------:|
| `int`  | `int` | `int` | `double` |       `double`       | `unsigned int` |

Compressed hdf5 files (filters) are supported via [hdf5_plugins]( https://github.com/HDFGroup/hdf5_plugins) and must be found by the library at runtime (e.g., by setting the environment variable `HDF5_PLUGIN_PATH`).
The module is capable of defining an event as well as adding records based on timestamp or trigger. In case of the latter, trigger information must be present in the events.

### Parameters
* `filename`: Input file name.
* `dataset_name`: Name of the node in the hdf5 file.
* `buffer_depth`: Buffer size (entries) for chunking. Default is 100,000.
* `event_length`: Duration of the event if this module is the first event loader and defines the event. Defaults to `1 us`.
* `sync_by_trigger`: Add records to the clipboard based on its trigger instead of timestamp. This requires an event definition with trigger information and can therefore not be used as first event loader.
* `timestamp_shift`: Shift the timestamp of the record by the defined value in nanoseconds.
* `trigger_shift`: Shift the trigger of the record by the defined value.

### Plots produced

The following plots are produced:

* Histogram with pixel raw value
* Histogram with pixel charge
* 2D map of hit positions
* 2D map of raw values per pixel

### Usage
```toml
[EventLoaderHDF5]
filename = "path/to/file.h5"
dataset_name = "Hits"
```
