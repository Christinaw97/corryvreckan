AnalysisTrackIntercept
**Maintainer**: Philipp Windischhofer (<philipp.windischhofer@cern.ch>)
**Module Type**: *GLOBAL*
**Status**: Functional

### Description
This module produces 2D histograms of the `x` and `y` coordinates at which reconstructed tracks intercept a given `z`-plane (or `z`-planes). 

### Parameters
* `n_bins_x`: Number of bins to be used along the `x` direction.
* `xmin`, `xmax`: Domain of histogram along the `x` direction.
* `n_bins_y`: Number of bins to be used along the `y` direction.
* `ymin`, `ymax`: Domain of histogram along the `y` direction.
* `plane_z`: Array of `z` coordinates specifying the planes for which track intercepts should be histogrammed.

### Usage
```
[AnalysisTrackIntercept]
n_bins_x = 100
xmin = 0 
xmax = 30
n_bins_y = 100
ymin = 0 
ymax = 15 
plane_z = [30, 60]
```
