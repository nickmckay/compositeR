# compositeR 0.3.0

* compositeR now builds on ens and lipdViz instead of the geoChronR monolith, so age modeling software is no longer installed with it.
* Fixed `scaleComposite()`, which passed the old `binvec` argument name to `bin()`.
* Added a `summary()` method for composites and deprecated the old API.
* lipdR moved from Suggests to Imports, since `compositeEnsembles()` calls it unconditionally. Requires lipdR >= 0.7.0.
