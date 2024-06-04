---
# SPDX-FileCopyrightText: 2017-2023 CERN and the Corryvreckan authors
# SPDX-License-Identifier: CC-BY-4.0 OR MIT
---
# EventLoaderTimepix4
**Maintainer**: Uwe Kraemer (uwe.kraemer@nikhef.nl)  
**Module Type**: *DETECTOR*  
**Detector Type**: *Timepix4*  
**Status**: Experimental

### Description
This module loads raw data from a Timepix4 device and adds it to the clipboard. The input files must have extension `.dat` and are sorted into time order into for both halves. 


The hit timestamps are derived from the 40 MHz ToA counter, the fast on-pixel 640 MHz Voltage Controlled Oscillator (VCO) which measures the number of cyclles completed until end of 40 MHz cycle called fine ToA (fToA), and the 4 phase shifted fToA values resulting in an ultra fine time of arrival (ufToA)

This module requires either another event loader of another detector type before which defines the event start and end times (Event object on the clipboard) or an instance of the Metronome module which provides this information.

The calibration is performed as described in and requires a Timepix4 plane to be set as `role = DUT`.

Currently calibration such as ToT->Charge, VCO calibration and Timewalk calibration are not supported and to be added later.

### Parameters
* `input_directory`: Path to the directory above the data directory for each device. The device name is added to the path during the module.

### Plots produced

For all detectors, the following plots are produced:

* 2D map of pixel positions
* raw ToA, extended ToA and hittime for each pixel hit
* Histogram with pixel ToT in different states (raw, full, corrected)

### Usage
```toml
[Timepix4EventLoader]
input_directory = "path/to/directory"
```
