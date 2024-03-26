# AnalysisElectronCT
**Maintainer**: Paul Schuetze (paul.schuetze@desy.de)
**Module Type**: *DETECTOR* **Detector Type**: *all*  
**Status**: Functional

### Description
This is a module analyzing data for electronCT measurements. It interprets individual events of individual detector planes as beam profiles and analyses these profiles with respect to their center position, widths and the total charge.

As a result, it produces a single cluster per event representing the beam properties. Note that these clusters do not correspond well to clusters created in other clustering algorithms as this algorithm does not consider the proximity of pixel hits but assigns all pixel hits of an event to a single cluster.

### Parameters
* `charge_weighting`: Use pixel charge as weighting parameter for profile measurements. Defaults to `true`.
* `fitted_profile`: Use fit to gaussian distribution instead of statistical evaluation for profile quantification. Defaults to `false`.
* `plot_frames`: Initial axis scaling, number of events to be plotted on event-dependent graphs. Defaults to `1000`.

### Plots produced
For each detector the following plots are produced:

* Hit maps and projections onto x- and y-axis
* Number of hits and charge per frame
* Histogram of center positions and widths of the beam spot
* All of the above as a function of the frame number

### Usage
```toml
[AnalysisElectronCT]

```
