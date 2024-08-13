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
