# data-raw

`Temp12k.R` rebuilds the package datasets `Temp12k` and `Temp12k_ageens`.

`Temp12k_ageens` is built from 698 LiPD files with age ensembles. They used to be kept here, but they made the repository too large to install with `remotes::install_github()`. They are identical to the files in the Temperature12k repository at commit `ea85457`:

https://github.com/nickmckay/Temperature12k/tree/ea854570be681680dc724fd6b2f31c0f4667483c/ScientificDataAnalysis/lipdFilesWithEnsembles

To rebuild `Temp12k_ageens`, put them in `data-raw/lipdFilesWithEnsembles/` (git ignores that folder), for example with a sparse checkout, which downloads only that folder (about 610 MB):

```sh
git clone --filter=blob:none --no-checkout https://github.com/nickmckay/Temperature12k.git /tmp/Temperature12k
git -C /tmp/Temperature12k sparse-checkout set --no-cone /ScientificDataAnalysis/lipdFilesWithEnsembles/
git -C /tmp/Temperature12k checkout ea854570be681680dc724fd6b2f31c0f4667483c
mv /tmp/Temperature12k/ScientificDataAnalysis/lipdFilesWithEnsembles data-raw/
```

The same files are also in this repository's history, up to commit `65d2fc0`.
