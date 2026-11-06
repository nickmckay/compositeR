# compositeR 0.5.0

* `compositeEnsembles()` now composites on each record's `time` by default (`ageVar = "time"`), converted to a common scale with the new `timeUnits` argument (default `"yr BP"`). It uses the time metadata that lipdR 0.7 adds in `extractTs()`, so records in year CE (tree rings, corals, documents), yr BP and ka can be composited together without converting them first. `binvec` is in `timeUnits`. A record's age ensemble (`ageEnsemble`) is used when it has one, except with `binFun = simpleBinTs`. Requires lipdR >= 0.7.3.
* Any other `ageVar` (for example `"age"` or `"ageEnsemble"`) is used as it is, as before. Results for Temp12k-style inputs, whose `age` is already in yr BP, are unchanged.
* Composites record their `timeUnits`. `plotComposite()` and `printComposite()` label the time axis with them, and year CE axes run forward instead of reversed.
* `simpleBinTs()` and `sampleEnsembleThenBinTs()` accept a one-row timeseries tibble as well as a TS element.

# compositeR 0.4.0

* `compositeEnsembles2()` is now `compositeEnsembles()`. The old name still works but warns that it is deprecated. The old single-ensemble `compositeEnsembles()`, deprecated in 0.3.0, has been removed; scripts that need it (for example the Temperature12k analysis code) should use compositeR 0.3.0.
* `compositeEnsembles()` accepts a timeseries tibble (from `lipdR::ts2tibble()`, or any data frame with one row per timeseries) as well as a TS list.
* `alignInterpDirection` works with current LiPD files. It used to look only for `*_interpDirection` fields, so records with `interpretation1_direction` were never flipped. The new `interpDirectionSign()` reads `interpretationN_direction` and `interpretationN_scope` (and the legacy fields), accepts positive/negative and increase/decrease in any case, and treats ambiguous values as unknown.
* New `uncUnits` argument to `compositeEnsembles()` and `sampleEnsembleThenBinTs()`. With `uncUnits = "sd"`, a numeric uncertainty is a multiple of each record's standard deviation (e.g. `uncVar = 0.5, uncUnits = "sd"`), so records in different units no longer need rescaling first. The default, `"native"`, keeps the old behavior. Uncertainties read from the timeseries are always in native units.
* `sampleEnsembleThenBinTs()` no longer simulates BAM age uncertainty for precisely dated archives. A new `exactAgeArchives` argument (default `c("Wood", "Documents")`) lists the archive types whose single age vector is used as is.

# compositeR 0.3.0

* compositeR now builds on ens and lipdViz instead of the geoChronR monolith, so age modeling software is no longer installed with it.
* Fixed `scaleComposite()`, which passed the old `binvec` argument name to `bin()`.
* Added a `summary()` method for composites and deprecated the old API.
* lipdR moved from Suggests to Imports, since `compositeEnsembles2()` calls it unconditionally. Requires lipdR >= 0.7.0.
