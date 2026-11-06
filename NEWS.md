# compositeR 0.4.0

* `compositeEnsembles2()` is now `compositeEnsembles()`. The old name still works but warns that it is deprecated. The old single-ensemble `compositeEnsembles()`, deprecated in 0.3.0, has been removed; scripts that need it (for example the Temperature12k analysis code) should use compositeR 0.3.0.
* `compositeEnsembles()` accepts a timeseries tibble (from `lipdR::ts2tibble()`, or any data frame with one row per timeseries) as well as a TS list.
* `alignInterpDirection` works with current LiPD files. It used to look only for `*_interpDirection` fields, so records with `interpretation1_direction` were never flipped. The new `interpDirectionSign()` reads `interpretationN_direction` and `interpretationN_scope` (and the legacy fields), accepts positive/negative and increase/decrease in any case, and treats ambiguous values as unknown.
* `sampleEnsembleThenBinTs()` no longer simulates BAM age uncertainty for precisely dated archives. A new `exactAgeArchives` argument (default `c("Wood", "Documents")`) lists the archive types whose single age vector is used as is.

# compositeR 0.3.0

* compositeR now builds on ens and lipdViz instead of the geoChronR monolith, so age modeling software is no longer installed with it.
* Fixed `scaleComposite()`, which passed the old `binvec` argument name to `bin()`.
* Added a `summary()` method for composites and deprecated the old API.
* lipdR moved from Suggests to Imports, since `compositeEnsembles2()` calls it unconditionally. Requires lipdR >= 0.7.0.
