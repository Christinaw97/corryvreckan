---
# SPDX-FileCopyrightText: 2017-2023 CERN and the Corryvreckan authors
# SPDX-License-Identifier: CC-BY-4.0 OR MIT
---
# EventDisplay
**Maintainer**: Simon Spannagel (simon.spannagel@desy.de), Sara Ruiz Daza (sara.ruiz.daza@desy.de)
**Module Type**: *DETECTOR*
**Status**: Immature

### Description
This module creates hitmaps of individual events.

It should be noted that this is only useful for a small number of events. It is suggested to use this module only after filtering interesting events using the `FilterEvents` module.

### Parameters
No parameters are used from the configuration file.

### Plots produced
For each detector the following plots are produced:

* 2D histogram of pixel raw data.

### Usage
```toml
[EventDisplay]

```
