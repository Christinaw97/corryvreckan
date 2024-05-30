---
# SPDX-FileCopyrightText: 2017-2023 CERN and the Corryvreckan authors
# SPDX-License-Identifier: CC-BY-4.0 OR MIT
---
# EventDisplay
**Maintainer**: ruizdaza (fhl-ripper02.desy.de)
**Module Type**: *DETECTOR* **Detector Type**: *<add types here>*  
**Status**: Immature

### Description
This is a demonstrator module only, taking data every detector on the clipboard and plots the pixel hit positions.
It serves as template to create new modules.

### Parameters
No parameters are used from the configuration file.

### Plots produced
* Histogram of event numbers

For each detector the following plots are produced:

* 2D histogram of pixel hit positions

### Usage
```toml
[EventDisplay]

```
