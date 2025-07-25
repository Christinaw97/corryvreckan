---
# SPDX-FileCopyrightText: 2018-2024 CERN and the Corryvreckan authors
# SPDX-License-Identifier: CC-BY-4.0 OR MIT
---
# MaskCreator

**Maintainer**: Simon Spannagel (<simon.spannagel@cern.ch>)  
**Module Type**: *DETECTOR*  
**Detector Type**: *all*  
**Status**: Work in progress

### Description

This module reads in `pixel` objects for each device from the clipboard, and masks pixels considered noisy or dead.

Currently, three methods are available. The `localdensity` noise estimation method is taken from the Proteus framework [@proteus-repo] developed by Université de Genève.
It uses a local estimate of the expected hit rate to find pixels that are a certain number of standard deviations away from this estimate.
The second method, `frequency`, is a simple cut on a global pixel firing frequency which masks pixels with a hit rate larger than `frequency_cut` times the mean global hit rate.
The bare occupancy is adjusted to account for pixel area in detectors with non-identical pixels (e.g. `PixelModuleDetector` and `RadialStripDetector`).
The third method, `dead`, simply masks all pixels without any recorded hits and adds them to the mask file similar to the `mask_dead_pixels` flag. Unlike this flag it however ignores any sort of noise masking method. This can be useful when pixels in the chip were masked before data acquisition and appear as unresponsive in the data set.

The module appends the pixels to be masked to the mask files provided in the geometry file for each device.
If no mask file is specified there, a new file `mask_<detector_name>.txt` is created in the globally configured output directory.
Already existing masked pixels are maintained.
No masks are applied in this module as this is done by the respective event loader modules when reading input data.

### Parameters

* `method`: Select the method to evaluate dead or noisy pixels. Can be either `dead`, `localdensity` or `frequency`, where the latter is chosen by default.
* `frequency_cut`: Frequency threshold to declare a pixel as noisy, defaults to 50. This means, if a pixel exhibits 50 times more hits than the average pixel on the sensor, it is considered noisy and is masked. Only used in `frequency` mode.
* `bins_occupancy`: Number of bins for occupancy distribution histograms, defaults to 128.
* `density_bandwidth`: Bandwidth for local density estimator, defaults to `2` and is only used in `localdensity` mode.
* `sigma_above_avg_max`: Cut for noisy pixels, number of standard deviations above average, defaults to `5`. Only used in `localdensity` mode.
* `rate_max`: Maximum rate, defaults to `1`. Only used in `localdensity` mode.
* `mask_dead_pixels`: If `true`, the module will search for pixels without any recorded hits and add them to the mask file. Default is `false`.
* `write_new_config`: If `true` and the detector config did not previously hold a mask file, then the new mask file is added to the outgoing config. Default is `false`.
* `new_config_suffix`: If `write_new_config=true` and the detector config did not previously hold a mask file, add a suffix to the new mask file. Default is `""` (no suffix).
* `square_big_pixel_weight`: The occupancy for masking is adjusted for pixel area - however, in some detectors, cross-talk between neighbours can also scale linearly with pixel size, giving localised occupancies that increase by the square of the pixel area. Setting this option to `true` divides the bare occupancy by the area squared (rather than just the area) to account for this. Defaults to `false`.

### Plots produced

For each detector the following plots are produced:

* Map of masked pixels
* 2D histogram of occupancy

### Usage

```toml
[MaskCreator]
frequency_cut = 10
```

[@proteus-repo]: https://gitlab.cern.ch/proteus/proteus
