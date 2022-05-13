# AnalysisElectronCT
**Maintainer**: Paul Schuetze (paul.schuetze@desy.de)
**Module Type**: *DETECTOR* **Detector Type**: *all*  
**Status**: Immature

### Description
This is a module analyzing data for electron CT measurements. It interprets individual events (or frames) of individual detector planes as beam profiles and analyses these profiles with respect to the center position, widths and the number of hit pixels.

### Parameters
* `charge_weighting`: Use charge as weighting for profile measurements. Defaults to `true`.
* `fitted_profile`: Use fit to gaussian distribution instead of statistical evaluation for profile quantification. Defaults to `false`.
* `plot_frames`: Initial axis scaling, number of events to be plotted on event-dependent graphs. Defaults to `1000`.

### Plots produced
For each detector the following plots are produced:

* Hit maps and projections onto x- and y-axis
* Number of hits per frame
* Center positions and widths of the beam spot
* All above - as a function of the frame number

### Usage
```toml
[AnalysisElectronCT]

```
