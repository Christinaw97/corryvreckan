---
# SPDX-FileCopyrightText: 2023-2024 CERN and the Corryvreckan authors
# SPDX-License-Identifier: CC-BY-4.0 OR MIT
---
# EventLoaderHDF5
**Maintainer**: Christian Bespin (<cbespin@uni-bonn.de>)
**Module Type**: *DETECTOR*  
**Status**: Work in progress

### Description
This module loads data from hdf5 files and adds it to the clipboard. The input file must have extension `.h5` without any compression filters and follow the structure below:

| column |  row  | charge |       timestamp      | trigger_number |
|:------:|:-----:|:------:|:--------------------:|:--------------:|
| `int`  | `int` | `int`  | `unsigned long long` | `unsigned int` |

Decimal values for the charge are not supported yet and both raw and charge of the `Pixel` object are populated with the same value.

The module is capable of defining an event as well as adding records based on timestamp or trigger. In case of the latter, trigger information must be present in the events.

### Parameters
* `filename`: Input file name.
* `dataset_name`: Name of the node in the hdf5 file.
* `buffer_size`: Buffer size for chunking.
* `event_length`: Duration of the event if this module is the first event loader and defines the event. Ignored otherwise.
* `sync_by_trigger`: Add records to the clipboard based on its trigger instead of timestamp. This requires an event definition with trigger information and can therefore not be used as first event loader.
* `timestamp_shift`: Shift the timestamp of the record by the defined value in nanoseconds.
* `trigger_shift`: Shift the trigger of the record by the defined value.

### Plots produced

The following plots are produced:

* 2D map of pixel positions
* Histogram with pixel charge

### Usage
```toml
[EventLoaderHDF5]
input_directory = "path/to/file"
dataset_name = "Hits"
```
