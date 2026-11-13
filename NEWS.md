# compositeR 0.6.2

* Requires lipdR 0.8.4 or later, which includes the `pullTsVariable()` fix from 0.8.3 and `filterByInterpretation()`.

# compositeR 0.6.1

* Compositing with spreading (`spread = TRUE`, the default) now gives the same composite whatever the time scale: the same records composited on yr BP, ka or year CE agree after converting the axes. Before, results on a year CE axis were close to, but not the same as, those on yr BP or ka. Results on a yr BP axis are unchanged.
  * `spreadPaleoData()` now always spreads on a yr BP axis, converting other scales to yr BP and back, so ties between two equally near values break the same way on every scale. `spreadAge` is still returned in `time.units`, in increasing order.
  * `simpleBinTs()` (and so `sampleEnsembleThenBinTs()` and `compositeEnsembles()`) puts a value that falls exactly on a bin edge in the younger bin on every scale. On a year CE axis it used to go in the older bin. `scaleComposite()` does the same when it knows the composite's `time.units`.
* `spreadPaleoData()` now honors `maxAge`, which it used to ignore. `minAge` and `maxAge` are the youngest and oldest times kept, in `time.units`: on a year CE axis `maxAge` is the earliest year kept and `minAge` the latest. (In 0.6.0, `minAge` on a year CE axis was the earliest year kept.) The output still never extends past the oldest end of the record.

# compositeR 0.6.0

Composites can be made on any time scale and any time variable, using the time scales lipdR 0.8 stores on every time column. Requires lipdR >= 0.8.1.

* `compositeEnsembles()` arguments are now `time.variable` and `time.units`, matching the other LiPD packages. `ageVar` and `timeUnits` still work, with a deprecation warning.
* `time.variable` can be any time variable of the records (e.g. `"age"`, `"year"`, or a column such as `"calAge"`; see `lipdR::tsTimeVariables()`), and it is always converted to `time.units`. Before, any variable other than `"time"` was used without conversion; for records whose variable is already on `time.units` (e.g. Temp12k `age` in yr BP), results are unchanged.
* `time.units` now defaults to the records' own time scale when they all share one, and to yr BP otherwise (it was always yr BP). For example, a set of tree-ring records in year CE is now composited in year CE, with a message saying so; `binvec` is read in that scale. Set `time.units` to choose. A `binvec` that doesn't overlap the records (usually one on the wrong scale) is now an error.
* Records whose time scale was assumed (units that don't say) give one warning per call, listing them.
* The output records `timeVariableName` as well as `timeUnits`, and `parameters` records `time.variable` and `time.units` (instead of `ageVar` and `timeUnits`). Normalized labels are used, e.g. `"yr AD"` for `"CE"`.
* `plotComposite()` gains `time.units`, which converts the axis to another scale (`ageUnits`, which only relabeled it, is deprecated). The axis direction comes from the scale.
* `simpleBinTs()` and `sampleEnsembleThenBinTs()`: `ageVar` is now `time.variable` (old name deprecated), and a new `time.units` gives the scale of that variable.
* `spreadPaleoData()` works on any time scale, with a new `time.units` argument (default "yr BP"). The default `minAge` of -75 (2025 CE in yr BP) dropped every year before 75 BCE on a year CE axis; the default is now no lower limit, but nothing after 2025 CE. The spreading grid is rounded to whole years on the scale, so ka axes work (it used to round to whole ka).
* `standardizeOverRandomInterval()` picks its random start from steps of a year, or of the bin spacing if that is finer, so it works on ka axes.
* `scaleComposite()` accepts a composite object, and gains `time.units` and `scale.time.units` to put the scaling data on the composite's scale.

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
