# compositeR

A package to create composite timeseries. 

[Insert text about papers that use composites]

The composite R package is designed to apply 3 different methods:
- Standard Composite Calibration (SCC)
  - Mean value within a specified window is subtracted from each data point 
- Dynamic Composite Calibration (DCC)
  - Mean of each record iteratively adjusted to minimize differences among records.
- Composite Plus Scale (CPS)
  - DCC plus a scaling of the variance using prior information  

We do not have the functionality for PAI or GAM yet.

TODO:
- Add PAI, GAM
- Finish documenting standardize.R
- Look at scale.R
- Finalize DCC vignettes. Add additional vignettes?
- Finalize README

New:
- Plotting and printing function
- Added documentation
- Simplified binning.R

## Package family

compositeR is part of a family of interoperable paleogeoscience packages. It
builds on **[ens](https://github.com/nickmckay/ens)** for ensemble binning,
gaussianization, and age-ensemble simulation, and on
**[lipdViz](https://github.com/nickmckay/lipdViz)** for visualization.

| Package | Role |
|---|---|
| [ens](https://github.com/nickmckay/ens) | Ensemble methods: binning, uncertainty propagation, correlation, regression, PCA, spectra |
| [lipdViz](https://github.com/nickmckay/lipdViz) | Visualization of LiPD data and ensemble analyses |
| [geoChronR](https://github.com/nickmckay/geoChronR-chronOnly) | Age modeling (Bacon, Bchron, OxCal, BAM) |
| [actR](https://github.com/LinkedEarth/actR) | Abrupt-change detection |
| **compositeR** | Paleoclimate record compositing (this package) |

Compositing methods plug in through the `binFun` and `stanFun` arguments of
`compositeEnsembles2()`, so custom binning or standardization strategies can be
supplied without modifying the package.
